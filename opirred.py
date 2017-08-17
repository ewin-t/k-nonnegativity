import itertools as itr
import numpy as np
import collections

def combs(n):
    ''' list of n choose k combinations, indexed by k'''
    return [map(lambda x: set(x), list(itr.combinations(range(n), i))) for i in range(n)]

def partitions_raw(n):
    ''' iterable of all partitions of n. helper method. '''
    # base case of recursion: zero is the sum of the empty list
    if n == 0:
        yield []
        return
    # modify partitions of n-1 to form partitions of n
    for p in partitions_raw(n-1):
        yield [1] + p
        if p and (len(p) < 2 or p[1] > p[0]):
            yield [p[0] + 1] + p[1:]

def partitions(k, n):
    ''' partitions of k using the numbers 1 through n-1, having only n parts'''
    raw = list(partitions_raw(k))
    final = []
    for r in filter(lambda x: (x[-1] < n) and len(x) == n, raw):
        c = collections.Counter(r)
        final.append(zip(c.keys(),c.values()))
    return final

def dual(s):
    ''' dual of some list of n subsets of [n]'''
    d = []
    for tmp in range(0, len(s)):
        d.append(set())
    for elt in range(0, len(s)):
        for ind in s[elt]:
            d[ind].add(elt)
    return d

def permute(s, sigma):
    ''' permute a list of subsets by sigma
    that is, send [s_1,\ldots,s_n] to [sigma(s_1),\ldots,sigma(s_n)]''' 
    return map(lambda i: set(map(lambda j: sigma[j], i)), s)

def eq(s, t):
    ''' helper method to check for equality
    Dependent on s,t having same number of elements and s having unique elements '''
    for sub in s:
        if sub not in t:
            return False
    return True

def opirreds(n, k):
    ''' enumerate all op-irreducible matrices
    n is the size, k is the number of nonzero entries '''
    final = []
    c = combs(n)
    
    parts = partitions(k,n)
    count = len(parts)
    for lam in parts:
        count -= 1
        subfinal = []
        subfinal_perms = []
        lam_combs = [itr.combinations(c[p[0]], p[1]) for p in lam]
        for tmp in itr.product(*lam_combs):
            m = []
            for k in tmp:
                m.extend(k)
            d = dual(m)
            good = True
            for i,j in itr.combinations(range(0, n), 2):
                if m[i].issubset(m[j]) or m[j].issubset(m[i]):
                    good = False
                    break
                if d[i].issubset(d[j]) or d[j].issubset(d[i]):
                    good = False
                    break
            if good:
                for t in subfinal_perms:
                    if eq(m, t):
                            good = False
                            break
                    if not good:
                        break
            if good:
                subfinal.append(m)
                subfinal_perms.extend(permute(m, sigma) for sigma in itr.permutations(range(n), n))
        final = final + subfinal
    print "length", len(final)
    return final


