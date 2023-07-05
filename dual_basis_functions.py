from sympy import *
import itertools
from collections import deque

n = 4
var('q')
var('r')

f = open(f"so{2*n}_eBases.txt", 'w')

H = symbols('H(1:5)', commutative = False) # These are 0-indexed (so H[0] is H1)
E = symbols('E(1:5)', commutative = False)
F = symbols('F(1:5)', commutative = False)

print(f"q = {q} \t n = {n}", file = f)

def a(i, j): # Cartan matrix lookup
    if (i == j):
        return 2
    if (abs(i-j) == 1 and max(i, j) <= n-1):
        return -1
    if (i == n-2 and j == n or i == n and j == n-2):
        return -1
    return 0

def pair(list1, list2):
    ret = []
    if(len(list1) != len(list2)):
        return []
    if(len(list1) == 1):
        if(list1[0] == list2[0]):
            return [0]
        else:
            return []
    first = list2[0]
    for i in range(len(list1)):
        if(list1[i] == first):
            inductive_list1 = list1.copy()
            inductive_list1.pop(i)
            inductive_list2 = list2.copy()
            inductive_list2.pop(0)
            ih = pair(inductive_list1, inductive_list2)
            prefactors = 0
            for j in range(i):
                prefactors += a(list1[j], list1[i])
            for j in range(len(ih)):
                ih[j] += prefactors
            ret.extend(ih)
    return ret

def mat(listofLists):
    M = zeros(len(listofLists))
    for i in range(len(listofLists)):
        for j in range(len(listofLists)):
            lst = pair(listofLists[i], listofLists[j])
            for k in lst:
                M[i, j] += q**k
    return M

def perm(tentlist): # given list of lists, removes linear dependence
    fin = []
    for i in range(len(tentlist)):
        fin1 = list(fin)
        fin1.append(tentlist[i])
        M = mat(fin1)
        if(M.det() != 0):
            fin.append(tentlist[i])
    return fin

def result(setofindices): # gives a basis of elements (the first element is the basis is setofindices)
    tentlist = list(dict.fromkeys(itertools.permutations(setofindices)))
    for i in range(len(tentlist)):
        tentlist[i] = list(tentlist[i])
    tentlist = reduce(tentlist)
    return perm(tentlist)

def reduce(tentlist): # remove duplicates
    visited=[0 for i in range(len(tentlist))]
    adj = [[0 for i in range(len(tentlist))] for j in range(len(tentlist))]
    finlist = []
    for n in range(len(tentlist)):
        l = tentlist[n]
        for ind in range(len(l) - 1):
            if(a(l[ind], l[ind+1]) == 0):
                ledit = l[:]
                ledit[ind], ledit[ind+1] = ledit[ind+1], ledit[ind]
                m = tentlist.index(ledit)
                adj[m][n] = 1
                adj[n][m] = 1
    adjacents = [[] for i in range(len(tentlist))]
    for i in range(len(tentlist)):
        temp = []
        for j in range(len(tentlist)):
            if(adj[i][j] == 1):
                temp.append(j)
        adjacents[i] = temp
    for n in range(len(tentlist)):
        if(visited[n] == 0):
            finlist.append(tentlist[n])
            qu = deque([n])
            while(qu):
                curr = qu.popleft()
                for neigh in adjacents[curr]:
                    if(visited[neigh] == 0):
                        visited[neigh] = 1
                        qu.append(neigh)
    return finlist 

# [1, 2, 3] -> E1 E2 E3
def indicesToE(setofindices):
    solution = 1
    for index in setofindices:
        solution = solution * E[index-1]
    return solution

# [1, 2, 3] -> F1 F2 F3
def indicesToF(setofindices):
    solution = 1
    for index in setofindices:
        solution = solution * F[index-1]
    return solution


def printEBasis(case, i, j):
    if (case == 1):
        pathSet = [x for x in range(i, j)] # e_{μλ} = E_{pathSet} = E_{i, ..., j - 1}
        eBasis = result(pathSet)
        print("+i: ", i, ", +j: ", j)
        print(f"+i: {i} \t +j: {j} \t eBasis: {eBasis}", file = f)
    if (case == 2):
        pathSet = [] # e_{μλ} = +- E_{pathSet} = +- E_{i, ..., n-2, n, n-1, ..., j} (if i < n) or +- E_{n, n-2, ..., j} (if i = n) or +- E_{i, ..., n-2, n} (if j = n)
        for x in range(i,n-1):
            pathSet.append(x)
        pathSet.append(n)
        if (i != n and j != n):
            pathSet.append(n-1)
        for x in reversed(range(j,n-1)):
            pathSet.append(x) 
        eBasis = result(pathSet)
        print("+i: ", i, ", -j: ", j)
        print(f"+i: {i} \t -j: {j} \t eBasis: {eBasis}", file = f)
    if (case == 3):
        pathSet = [x for x in reversed(range(j, i))] # e_{μλ} = E_{pathSet} = E_{i-1, ..., j}
        eBasis = result(pathSet)
        print("-i: ", i, ", -j: ", j,)
        print(f"-i: {i} \t -j: {j} \t eBasis: {eBasis}", file = f)


def eBases():
    print("@@@@@@@@ eBases @@@@@@@@@@", file = f)
  
    for i in range(1,n+1):   # Computing dual elements e* and f* CASE 1
        for j in range(i+1,n+1): # μ = L_i > λ = L_j (i < j)
            printEBasis(1, i, j)
    
    for j in range(1,n+1): # CASE 3
        for i in range(j+1,n+1): # μ = -L_i > λ = -L_j (i > j)
            printEBasis(3, i, j)


    for i in reversed(range(1,n+1)): # CASE 2
        for j in reversed(range(1,n+1)): 
            if (i != n or j != n): # μ = L_i > λ = -L_j
                printEBasis(2, i, j)


eBases()


