# quiver.sage

import unicodedata as ud

""" Module to generate exchange graphs and sub-exchange graphs of cluster 
    algebras for total positivity and k-positivity """

def minor(mat,rows,cols):
    """Returns minor of mat, where rows from rows and cols from cols  """
    if len(rows)==0:
        return 1
    sgn = 1
    total = 0
    for i in range(len(rows)):
        minr = minor(mat,rows[:i]+rows[i+1:], cols[1:])
        total += sgn*mat[rows[i]][cols[0]]*minr
        sgn *= -1
    return total

# order is: 1,1  1,2  2,1  12,12  1,3  3,1  12,23  23,12  123,123

def exchange_poly(q,j,var_lst):
    """Computes the exchange polynomial which results from mutating quiver
        q at vertex j, with vertex variables accessed via var_lst"""
    m = 1
    n = 1
    mat = q.digraph().adjacency_matrix()
    for t in range(mat.dimensions()[0]):
        m *= (var_lst[t] if mat[j][t]==1 else 1)
        n *= (var_lst[t] if mat[t][j]==1 else 1)
    return m+n

def get_name(var_lab, ind_mat, n,dic):
    """Given an iterable of variable names, returns an alphabetized string"""
    def cmp(p,fn):
        return p==fn  #return maxima.ratsimp(p-fn)==0
    lst = []
    if n==3:
        det = minor(ind_mat, range(3), range(3))
        for i in range(len(var_lab)):
            p = var_lab[i]
            name = "!!"
            if cmp(p,minor(ind_mat, [1,2], [1,2])):
                name = "A"
            elif cmp(p,minor(ind_mat, [1,2], [0,2])):
                name = "B"
            elif cmp(p, minor(ind_mat, [0,2], [1,2])):
                name = "D"
            elif cmp(p, minor(ind_mat, [0,2], [0,2])):
                name = "E"
            elif cmp(p, minor(ind_mat, [0,2], [0,1])):
                name = "F"
            elif cmp(p, minor(ind_mat, [0,1], [0,2])):
                name = "H"
            elif cmp(p, minor(ind_mat, [0,1], [0,1])):
                name = "J"
            elif cmp(p, ind_mat[0][0]):
                name = "a"
            elif cmp(p, ind_mat[0][1]):
                name = "b"
            elif cmp(p, ind_mat[1][0]):
                name = "d"
            elif cmp(p, ind_mat[1][1]):
                name = "e"
            elif cmp(p, ind_mat[1][2]):
                name = "f"
            elif cmp(p, ind_mat[2][1]):
                name = "h"
            elif cmp(p, ind_mat[2][2]):
                name = "j"
            elif cmp(p,(ind_mat[2][2]*minor(ind_mat, [0,1],[0,1]) - det)):
                name ="K"
            elif cmp(p, (ind_mat[0][0]*minor(ind_mat, [1,2],[1,2]) - det)):
                name = "L"
            lst.append(name)
    else:        
        for p in var_lab:
            try:
                name = dic[p]
            except:
                name = dic["extra"][0]
                dic[p] = name
                dic["extra"]=dic["extra"][1:]
                dic[name] = p
            lst.append(name) #lst.append(str(maxima.ratsimp(p)))
    lst.sort()
    name = ""
    for c in lst:
        name = name + c
    return name

def has_k(q,j, num_mat, n, k):
    """Checks whether the exchange polynomial of quiver q at vertex j
        'goes through' a minor of size greater than k."""
    d = q.digraph()
    for i in range(k,n):
        for l in range(k,n):
            if d.has_edge(j, num_mat[i][l]) or d.has_edge(num_mat[i][l],j):
                return True
    return False 

def highest_k(q,j,num_mat,n):
    d = q.digraph()
    h = 0
    for i in range(n):
        for l in range(n):
            if d.has_edge(j,num_mat[i][l]) or d.has_edge(num_mat[i][l],j):
                h = max(h,min(i,l)+1)
    return h

def get_clust(labels, unfrozen_inds):
    return tuple([labels[i] for i in unfrozen_inds])

def exchange_graph(n=3, all_comps=False, k=2):
    """ Returns tuple (G,d), where G is a graph corresponding to the nxn 
        cluster algebra, and d is a dictionary giving the polynomial for 
        each letter in the vertex label.
        all_comps is True if all connected components should be found,
        False if only the component with the lex minimal one should be
        if all_comps == False, k should also be provided.
        Note that for n>3, the exchange graph is infinite.
        """
    R = PolynomialRing(ZZ,n,n,var_array="x")
    ind_mat = [[R("x"+str(i)+str(j)) for j in range(n)] for i in range(n)]
    num_mat = [ range(i*n,i*n+n) for i in range(n)]
    all_unicode = ''.join(unichr(i) for i in xrange(65536))
    unicode_letters = ''.join(c for c in all_unicode
                          if ud.category(c)=='Lu' or ud.category(c)=='Ll')
    dic = {"extra":unicode_letters}

    d = []
    for i in range(n-1,0,-1):
        d.extend([ [num_mat[i][j],num_mat[i-1][j]] for j in range(n)]) # up ->
        d.extend([ [num_mat[j][i],num_mat[j][i-1]] for j in range(n)]) #left ->
    for i in range(n-1): # diag -> 
        d.extend([ [num_mat[i][j],num_mat[i+1][j+1]] for j in range(n-1)])
    quiv_lst = [ClusterQuiver(d)]
    lab = []
    frozen_inds = []
    unfrozen_inds = []
    test_lab = []
    for i in range(n):
        for j in range(n):
            m = min(i,j)
            lab.append(minor(ind_mat, range(i-m,i+1), range(j-m,j+1)))
            flag = highest_k(quiv_lst[0],num_mat[i][j],num_mat,n) > k
            if i==n-1 or j==n-1 or (not all_comps and flag):
                frozen_inds.append(num_mat[i][j]) 
            else:
                unfrozen_inds.append(num_mat[i][j])
    all_labels = [lab]
    lab = get_clust(all_labels[0], unfrozen_inds)
    var_names = [get_name(lab, ind_mat, n, dic)]

    G=Graph()
    i=0
    while i < len(quiv_lst):
        G.add_vertex(name=var_names[i])
        for j in unfrozen_inds:
            q = quiv_lst[i].mutate(j,False)
            curr = all_labels[i]
            div = exchange_poly(q,j,curr).quo_rem(curr[j])[0]
            lab = [curr[l] if l!= j else div for l in range(n**2)]  
            name = get_clust(lab, unfrozen_inds)
            name = get_name(name, ind_mat, n, dic)
            if name not in var_names:
                quiv_lst.append(q)
                all_labels.append(lab)
                var_names.append(name)
            G.add_edge(var_names[i],name,highest_k(q,j,num_mat,n))
        i += 1
    return (G,dic)

def sub_exchange(G,k=2):
    """ Returns the subgraph of G avoiding edges which go through minors
        of size >k. Requires: G is the result of exchange_graph(n) for
        some n>=k.
        Don't actually try to use this for n>3 since then the exchange graph
        is infinite.
    """
    H = Graph()
    H.add_vertices(G.vertices())
    for e in G.edge_iterator():
        if e[2]<=k:
            H.add_edge(e)
    return H
