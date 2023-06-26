from sympy import *
#from central_element_functions import *
import itertools
from collections import deque
import json

'''
with open('eDualSymbolic.json', 'r') as fSymbolic:
    eDualSymbolicData = json.load(fSymbolic) 
with open('eDualNumer.json', 'r') as fNumerous:
    eDualNumerousData = json.load(fNumerous)
with open('eBasis.json', 'r') as fEBasis:
    eBasisData = json.load(fEBasis) 
'''

eDualSymbolicData = {}
eDualNumerousData = {}
eBasisData = {}

q = Symbol('q')




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

def mat(listofLists, q):
    M = zeros(len(listofLists))
    for i in range(len(listofLists)):
        for j in range(len(listofLists)):
            lst = pair(listofLists[i], listofLists[j])
            for k in lst:
                M[i, j] += q**k
    return M.applyfunc(simplify)

def dual(listofLists, q): # for computing specific dual elements
    M = zeros(len(listofLists))
    for i in range(len(listofLists)):
        for j in range(len(listofLists)):
            lst = pair(listofLists[i], listofLists[j])
            for k in lst:
                M[i, j] += q**k
    N= M.inv()
    return N.applyfunc(simplify)

def perm(tentlist, q): # given list of lists, removes linear dependence
    fin = []
    for i in range(len(tentlist)):
        fin1 = list(fin)
        fin1.append(tentlist[i])
        M = mat(fin1, q)
        if(M.det() != 0):
            fin.append(tentlist[i])
    return fin

def result(setofindices, q): # gives a basis of elements (the first element is the basis is setofindices)
    tentlist = list(dict.fromkeys(itertools.permutations(setofindices)))
    for i in range(len(tentlist)):
        tentlist[i] = list(tentlist[i])
    tentlist = reduce(tentlist)
    return perm(tentlist, q)

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

def getPathSet(n, case, i, j):
    if (case == 1):
        pathSet = [x for x in range(i, j)] # e_{μλ} = E_{pathSet} = E_{i, ..., j - 1}
    elif (case == 2):
        pathSet = [] # e_{μλ} = +- E_{pathSet} = +- E_{i, ..., n-2, n, n-1, ..., j} (if i < n) or +- E_{n, n-2, ..., j} (if i = n) or +- E_{i, ..., n-2, n} (if j = n)
        for x in range(i,n-1):
            pathSet.append(x)
        pathSet.append(n)
        if (i != n and j != n):
            pathSet.append(n-1)
        for x in reversed(range(j,n-1)):
            pathSet.append(x)
    else:
        pathSet = [x for x in reversed(range(j, i))] # e_{μλ} = E_{pathSet} = E_{i-1, ..., j}
    return pathSet
         
    

# Takes in e and returns a list of the coefficients in e* (numerical)
# Case I: e_{L_i, L_j} (i < j)
# Case II: e_{L_i, -L_j}
# Case III: e_{-L_i, -L_j} (j > i)
def eDualNumerical(n, case, i, j, qNum):
    if ((n, case, i, j, qNum) in eDualNumerousData): 
        return eDualNumerousData[(n, case, i, j, qNum)]
    else:
        pathSet = getPathSet(n, case, i, j)
        if ((n, case, i, j) in eBasisData):
            eBasis = eBasisData[(n, case, i, j)]
        else:
            eBasis = result(pathSet, qNum) 
            eBasisData[(n, case, i, j)] = eBasis
        eDualCoeffs = dual(eBasis, qNum).row(0).tolist()[0] 
        eDualNumerousData[(n, case, i, j, qNum)] = eDualCoeffs
        eBasisData[(n, case, i, j)] = eBasis
        return eDualCoeffs  

# Takes in e and returns a list of the coefficients in e* (symbolic)
# Case I: e_{L_i, L_j} (i < j)
# Case II: e_{L_i, -L_j}
# Case III: e_{-L_i, -L_j} (j > i)
def eDualSymbolic(n, case, i, j):
    if ((n, case, i, j) in eDualSymbolicData): 
        return eDualSymbolicData[(n, case, i, j)]
    else: 
        currCoeffs = eDualNumerical(n, case, i, j, 10)
        pts = [[] for dual in currCoeffs]
        for qNum in range(10, 10 + 2*n):
            if qNum != 10:
                currCoeffs = eDualNumerical(n, case, i, j, qNum)
            for dual in range(len(pts)):
                pts[dual].append((qNum, currCoeffs[dual]))
        
        symbolicCoeffs = []
        for pt in pts:
            symbolicCoeffs.append(interpolate(pt, q)) # interpolate to find the coefficients
        
        eDualSymbolicData[(n, case, i, j)] = symbolicCoeffs # add the coefficients to the dictionary
        return symbolicCoeffs


eDualSymbolic(3, 1, 1, 2)

with open('eDualSymbolic.json', 'w') as fSymbolic:
    json.dump(eDualSymbolicData, fSymbolic)
with open('eDualNumer.json', 'rw') as fNumerous:
    json.dump(eDualNumerousData, fNumerous)
with open('eBasis.json', 'r') as fEBasis:
    json.dump(eBasisData, fEBasis)
