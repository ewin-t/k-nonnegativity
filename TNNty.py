import numpy as np
import itertools as itr

def e(n, i, c):
    '''RHS: add c copies of column $i$ to column $i+1$
    LHS: add c copies of row $i+1$ to row $i$'''
    x = np.matrix(np.eye(n))
    x[i,i+1] = c
    return x

def f(n, i, c):
    '''RHS: add c copies of column $i+1$ to column $i$
    LHS: add c copies of row $i$ to row $i+1$'''
    x = np.matrix(np.eye(n))
    x[i+1,i] = c
    return x

def k_minors(X, k, callzeros = False, checkrows = False, prec=0):
    '''Used by other functions to compute value of minors
    if callzeros = true then it prints out all of the nonobvious 0 minors
    if checkrows = true then it checks the row-solid minors
    if prec != 0 then the determinant is rounded to prec places before checking'''
    n = X.shape[0]
    zero_minors = []
    good = True
    if checkrows:
        sets = itr.product((tuple(range(i, i+k)) for i in range(n - k + 1)), itr.combinations(range(n), k))
    else:
        sets = itr.product(itr.combinations(range(n), k), (tuple(range(i, i+k)) for i in range(n - k + 1)))

    for subset in sets:
        submatrix = X[np.array(subset[0])][:,np.array(subset[1])]
        minor = np.linalg.slogdet(submatrix)[0] if prec == 0 else np.round(np.linalg.det(submatrix), prec)
        if minor < 0:
            good = False
            break
        if callzeros and minor == 0 and (0 not in (submatrix != 0).sum(0)) and (0 not in (submatrix != 0).sum(1)):
            zero_minors.append(subset)

    if callzeros:
        return zero_minors
    else:
        return good

def is_TNNk(X, k, prec=0):
    '''Check if a matrix is k-nn
    X, k, callzeros, checkrows, prec'''
    for i in range(1, k+1):
        if not k_minors(X, i, prec=prec):
            return False
    return True

def max_TNNk(X, prec=0):
    '''Return the largest k for which X is k-nn
    X, k, callzeros=False, checkrows=False, prec=10'''
    for i in range(1, X.shape[0]+1):
        if not k_minors(X, i, prec=prec):
            return i-1
    return X.shape[0]

def zero_minors(X, k, prec=0, checkrows=False):
    n = X.shape[0]
    final = []
    for i in range(1, k+1):
        final.extend(k_minors(X, i, callzeros=True, checkrows=checkrows, prec=prec))
    return final

def generate_matrix(n, k, test_maximality = True, intsonly = True, high=20, zeros=0.1, triangular=False, max_check = 1e6):
    '''Generate k-nonnegative n by n matrices
    test_maximality: returns maximally k-nonnegative (default True)
    intsonly: matrix only has integer entries (default True)
    high: upper bound on values in matrix (default 20)
    zeros: more zeros = more zero blocks. fraction indicates probability, integer indicates that number. default 0.1
    triangular: returns upper triangular matrix. default false
    '''
    indices = np.triu_indices(n) if triangular else np.mask_indices(n, lambda x,k: np.ones(x.shape))
    total_count = 0
    while total_count < max_check:
        mat = np.zeros((n,n))
        if intsonly:
            mat[indices] = np.random.randint(1,high,len(indices[0]))
        else:
            mat[indices] = np.random.random_sample(len(indices[0]))*high
        if zeros > 0 and zeros < 1:
            for i,j in zip(indices[0], indices[1]):
                if i != j and np.random.rand() < zeros:
                    mat[i,j] = 0
        if zeros >= 1:
            nonzeros = np.nonzero(mat * (np.ones((n,n)) - np.eye(n)))
            zinds = np.random.choice(len(nonzeros[0]), zeros)
            for ind in zinds:
                mat[nonzeros[0][ind], nonzeros[1][ind]] = 0
        for i,j in itr.product(range(n), range(n)):
            if mat[i,j] == 0:
                irange = range(i, n) if i > j else range(0, i+1)
                jrange = range(0, j+1) if i > j else range(j, n)
                for i2, j2 in itr.product(irange, jrange):
                    mat[i2, j2] = 0
        mat = np.matrix(mat)
        good = (np.linalg.slogdet(mat)[0] != 0)
        knn = max_TNNk(mat)
        good = good and ((knn == k) if test_maximality else (knn >= k))
        if good:
            return mat
    print "total count exceeded max. returning..."
    return

def LDU(X):
    ''' Compute the LDU factorization '''
    n = X.shape[0]
    L = np.matrix(np.eye(n))
    D = np.matrix(np.eye(n))
    U = np.matrix(np.eye(n))
    princs = []
    for i in range(0,n):
        val = np.linalg.det(X[:i+1,:i+1])
        if val == 0:
            print "zero principal minor, aborting"
            return L,D,U
        princs.append(val)
    princs.append(1)
    # D
    for i in range(0, n):
        D[i,i] = princs[i] / princs[i-1]
    # L
    for j in range(0, n):
        jrange = range(0, j)
        for i in range(j+1, n):
            L[i,j] = np.linalg.det(X[np.array(jrange + [i])][:,:j+1]) / princs[j]
    # U
    for i in range(0, n):
        irange = range(0, i)
        for j in range(i+1, n):
            U[i,j] = np.linalg.det(X[:i+1][:,np.array(irange + [j])]) / princs[i]
    return L,D,U

def modwhile(X, k, mod, prec):
    while is_TNNk(mod(X), k):
        X = mod(X)
    return X

def irred_helper(X, k, c, prec):
    n = X.shape[0]
    Y = np.zeros(X.shape)
    while not np.array_equal(Y, X):
        Y = np.copy(X)
        for i in range(0, n-1):
            X = modwhile(X, k, (lambda x: x * e(n, i, -c)), prec)
            X = modwhile(X, k, (lambda x: e(n, i, -c) * x), prec)
            X = modwhile(X, k, (lambda x: x * f(n, i, -c)), prec)
            X = modwhile(X, k, (lambda x: f(n, i, -c) * x), prec)
    return X

def irred(X, k, prec):
    if not is_TNNk(X, k):
        print 'error; matrix is not k-nonnegative'
        return
    c_tmp = 1.0
    tmp = X
    while c_tmp > 1e-50:
        print tmp, c_tmp
        tmp = irred_helper(tmp, k, c_tmp, prec)
        c_tmp = c_tmp / 10
    return np.round(tmp, 40)


'''
def sixthree(a, b, c, d):
    mat = np.matrix(np.eye(6))
    for i in range(5):
        mat[i,i+1] = 1
        mat[i+1,i] = 1
    mat[0,0] = a
    mat[1,1] = b + 1/a
    mat[2,2] = c + 1/b
    mat[3,3] = 1/(mat[2,2]-1/mat[1,1])
    mat[4,4] = 1/(mat[3,3]-1/mat[2,2])
    mat[5,5] = d + 1/(mat[4,4] - 1/mat[3,3])
    return mat
    
def red(X, k):
    n = X.shape[0]
    for i in range(0, k):
        for j in range(0, k-i):
            if X[n-i-2,j] != 0:
                for entry in range(0, n):
                    X[n-j-1, entry] -= X[n-j-1, i]/X[n-j-2,i]*X[n-j-2,entry]
            if X[j,n-i-2] != 0:
                for entry in range(0, n):
                    X[entry, n-j-1] -= X[i, n-j-1]/X[i, n-j-2]*X[entry,n-j-2]
    return X

def generate_tri_matrices(n, k, num, test_maximality=True, intsonly = True, high=20):
    output_list = []
    count = 0
    total_count = 0
    while count < num:
        mat = np.zeros((n,n))
        mat[np.triu_indices(n)] = np.random.randint(1,high,(n*(n+1))/2) if intsonly else np.random.random_sample(n*(n+1)/2)*high
        total_count += 1
        for i in range(0, n):
            mat[i,i] = 1
        for i in range(0, n):
            for j in range(0, n):
                if mat[i,j] == 0:
                    irange = range(i, n) if i > j else range(0, i+1)
                    jrange = range(0, j+1) if i > j else range(j, n)
                    for i2, j2 in itr.product(irange, jrange):
                        mat[i2, j2] = 0
        if np.linalg.slogdet(mat)[0] != 0:
            knn = max_TNNk(mat)
            good = (knn == k) if test_maximality else (knn >= k)
            if good:
                count += 1
                print count
                output_list.append(mat)
    return output_list

'''
