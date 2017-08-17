import itertools as itr

# n = 4, k = 2 unitriangular
def orderr(x, y):
    if x[1] == y[1]:
        return x[0].bruhat_lequal(y[0])
    elif x[1] == 1:
        return False
    return x[0].bruhat_lequal(y[0]*Permutation([2,1,4,3])) or x[0].bruhat_lequal(y[0]*Permutation([1,4,2,3])) or x[0].bruhat_lequal(y[0]*Permutation([3,1,2,4]))

elms = list(itr.product(Permutations(4), range(2)))

P = Poset((elms, orderr))
P.is_eulerian()

heights = dict()
for perm in Permutations(4):
    l = float(perm.length())
    if l not in heights:
        heights[l] = []
    if l+2.5 not in heights:
        heights[l+2.5] = []
    heights[l].append(tuple([perm, 0]))
    heights[l+2.5].append(tuple([perm, 1]))

edges = {'black': [], 'red': []}
for x in P.cover_relations():
    if x[0][1] != x[1][1]:
        edges['red'].append(x)
    else:
        edges['black'].append(x)

P.show(figsize=[20.0, 15.0], partition = [[tuple([k, 0]) for k in Permutations(4)], [tuple([k,1]) for k in Permutations(4)]], fontsize = 20, heights=heights)

mod = DiGraph([elms, edges['red']])

####

P = Posets.SymmetricGroupBruhatIntervalPoset([1,2,3,4],[4,3,2,1])
relations = []
elements = []
not_elements = ((4,3,2,1), (3,4,1,2), (4,3,1,2))
cor_gees = ((2,3,4,1),(1,4,2,3),(1,4,2,3))
extra_elements = [tuple([k, 1]) for k in ((1,2,3,4),(2,1,3,4),(1,3,2,4),(1,2,4,3),(2,3,1,4),(1,4,2,3),(1,3,4,2),(2,1,4,3),(2,3,4,1))]
for k in Permutations(4):
    if k not in not_elements:
        elements.append(tuple([tuple(k), 0]))
elements.extend(extra_elements)
for rel in P.cover_relations():
    new_rel = []
    for k in rel:
        if k in not_elements:
            new_rel.append(tuple([cor_gees[not_elements.index(k)],1]))
        else:
            new_rel.append(tuple([k, 0]))
    relations.append(new_rel)
for gee in extra_elements:
    for g_zero in [Permutation([2,3,1,4]), Permutation([2,1,4,3]), Permutation([1,3,4,2])]:
        relations.append([tuple([Permutation(gee[0]).left_action_product(g_zero), 0]), gee])
