from sympy import symbols
from os import environ
import numpy as np
import numpy.linalg as la
import itertools
from collections import deque


n = 5
q = 10
perm_epsilon = 0.00000001
coeff_epsilon = q**(-2*n)

fname = f"so{2*n}_q_{q}_epsilon_{perm_epsilon}_np_test.txt"
f = open(fname, 'w')
p = open('perm.txt', 'w')
environ['OMP_NUM_THREADS'] = '6'


######## Generators ###################
H = symbols('H(1:6)', commutative = False) # These are 0-indexed (so H[0] is H1)
E = symbols('E(1:6)', commutative = False)
F = symbols('F(1:6)', commutative = False)
#eBasesCase2 = [[0 for j in range(n)] for i in range(n)] # Stores eBasis for different values of i and j (eBasis for i, j is the same as fBasis for j, i)
#eDualMatricesCase2 = [[0 for j in range(n)] for i in range(n)]
#eBasesCase1 = [[0 for j in range(n)] for i in range(n)] 
#eDualMatricesCase1 = [[0 for j in range(n)] for i in range(n)]

print(f"q = {q} \t n = {n}", file = f)


def inverse(M):
    detM = la.det(M)
    size = len(M)
    coM = np.zeros((size, size))
    signRow = 1
    for i in range(size):
        signCol = signRow       #
        for j in range(size):
            minor = np.delete(M, i, 0)
            minor = np.delete(minor, j, 1)
            coM[i,j] = signCol * la.det(minor)
            signCol *= -1
        signRow *= -1
    adjM = coM.T
    return (1/detM)*adjM    
  

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

def mat(listOfLists):
    l = len(listOfLists)
    M = np.zeros((l,l))
    for i in range(l):
        for j in range(l):
            lst = pair(listOfLists[i], listOfLists[j])
            for k in lst:
                M[i][j] += q**k
    return M

def dual(listOfLists): # for computing specific dual elements
    return la.inv(mat(listOfLists))

# Bareiss algorithm for determinant
def bareissDet(M):
    M = [row[:] for row in M] # make a copy to keep original M unmodified
    N, sign, prev = len(M), 1, 1
    for i in range(N-1):
        if M[i][i] == 0: # swap with another row having nonzero i's elem
            swapto = next( (j for j in range(i+1,N) if M[j][i] != 0), None )
            if swapto is None:
                return 0 # all M[*][i] are zero => zero determinant
            M[i], M[swapto], sign = M[swapto], M[i], -sign
        for j in range(i+1,N):
            for k in range(i+1,N):
                #assert ( M[j][k] * M[i][i] - M[j][i] * M[i][k] ) % prev == 0
                M[j][k] = ( M[j][k] * M[i][i] - M[j][i] * M[i][k] ) // prev
        prev = M[i][i]
    return sign * M[-1][-1]

def perm(tentlist): # given list of lists, removes linear dependence
    fin = []
    for i in range(len(tentlist)):
        fin1 = list(fin)
        fin1.append(tentlist[i])
        
        '''l = len(fin1)
        M = np.zeros((l,l))
        for i in range(l):
            mult = 0
            for j in range(l):
                lst = pair(fin1[i], fin1[j])
                for k in lst:
                    M[i][j] += q**k
                if min(lst) < mult:
                    mult = min(lst)
            for j in range(l):
                M[i][j] = int(round(M[i][j] * q**(-1 * mult)))'''
        M = mat(fin1)
        

        if (np.abs(la.det(M)) > perm_epsilon):
        #if (not np.allclose(la.det(M), 0)):
        #if (la.det(M) != 0):
            fin.append(tentlist[i])   
    return fin

def permTemp(tentlist): # given list of lists, removes linear dependence
    fin = []
    count = 1
    prevDet = 1
    for i in range(len(tentlist)):
        fin1 = list(fin)
        fin1.append(tentlist[i])
        M = mat(fin1)


        '''l = len(fin1)
        M = np.zeros((l,l))
        for i in range(l):
            mult = 0
            for j in range(l):
                lst = pair(fin1[i], fin1[j])
                for k in lst:
                    M[i][j] += q**k
                if min(lst) < mult:
                    mult = min(lst)
            for j in range(l):
                M[i][j] = int(round(M[i][j] * q**(-1 * mult)))'''
        
        currDet = la.det(M)
        if (np.abs(currDet / prevDet) > perm_epsilon):
        #if (not np.allclose(la.det(M), 0)):
        #if (la.det(M) != 0):
            #print(f"{count} M = {M}")
            fin.append(tentlist[i])
            print(f"{count} added {tentlist[i]}, det M = {currDet}\n", file=f)
            prevDet = currDet
        else:
            print(f"{count} rejected {tentlist[i]}, det M = {currDet}\n", file=f)
        count += 1
        print(count)
    print(f"final basis: {fin}", file = f)
    print(f"number of accepted terms: {len(fin)}", file = f)
    return fin

def removeperm(tentlist):
    count = 0
    for testlist in itertools.combinations(tentlist, 20):
        count += 1
        M = mat(testlist)

        if (np.abs(la.det(M)) > perm_epsilon):
        #if (not np.allclose(la.det(M), 0)):
        #if (la.det(M) != 0):
            #print(f"{count} M = {M}")
            return testlist
        
        if (count % 10000 == 0):
            print(count)
    print(f"failed", file = f)
    return tentlist


def getPathSet(n, i, j, c):
    match(c):
        case 1: pathSet = [x for x in range(i, j)] # e_{μλ} = E_{pathSet} = E_{i, ..., j - 1}
        
        case 2: 
            pathSet = [] # e_{μλ} = +- E_{pathSet} = +- E_{i, ..., n-2, n, n-1, ..., j} (if i < n) or +- E_{n, n-2, ..., j} (if i = n) or +- E_{i, ..., n-2, n} (if j = n)
            for x in range(i,n-1):
                pathSet.append(x)
            pathSet.append(n)
            if (i != n and j != n):
                pathSet.append(n-1)
            for x in reversed(range(j,n-1)):
                pathSet.append(x) 
        case 3: pathSet = [x for x in reversed(range(j, i))] # e_{μλ} = E_{pathSet} = E_{i-1, ..., j}
    
    return pathSet


def result(setofindices): # gives a basis of elements (the first element is the basis is setofindices)
    tentlist = list(dict.fromkeys(itertools.permutations(setofindices)))
    for i in range(len(tentlist)):
        tentlist[i] = list(tentlist[i])
    tentlist = reduce(tentlist)
    return perm(tentlist)

def resultTemp(setofindices): # gives a basis of elements (the first element is the basis is setofindices)
    tentlist = list(dict.fromkeys(itertools.permutations(setofindices)))
    for i in range(len(tentlist)):
        tentlist[i] = list(tentlist[i])
    tentlist = reduce(tentlist)
    return permTemp(tentlist)

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

# Takes in a path set for E and returns (e*, f*)
def dualElements(pathSet, i, j, case):
    '''if (case == 1):
        eBasis = result(pathSet)
        eDualMatrix = dual(eBasis)
        eBasesCase1[i-1][j-1] = eBasis
        eDualMatricesCase1[i-1][j-1] = eDualMatrix
        fBasis = result(list(reversed(pathSet)))
        fDualMatrix = dual(fBasis)
        eBasesCase1[j-1][i-1] = fBasis
        eDualMatricesCase1[j-1][i-1] = fDualMatrix
    if (case == 2):
        if (eBasesCase2[i-1][j-1] == 0):
            eBasis = result(pathSet)
            eDualMatrix = dual(eBasis)
            eBasesCase2[i-1][j-1] = eBasis
            eDualMatricesCase2[i-1][j-1] = eDualMatrix
        else:
            eBasis = eBasesCase2[i-1][j-1]
            eDualMatrix = eDualMatricesCase2[i-1][j-1]
        if (eBasesCase2[j-1][i-1] == 0):
            fBasis = result(list(reversed(pathSet)))
            fDualMatrix = dual(fBasis)
            eBasesCase2[j-1][i-1] = fBasis
            eDualMatricesCase2[j-1][i-1] = fDualMatrix
        else:
            fBasis = eBasesCase2[j-1][i-1]
            fDualMatrix = eDualMatricesCase2[j-1][i-1]
    if (case == 3):
        eBasis = eBasesCase1[i-1][j-1]
        eDualMatrix = eDualMatricesCase1[i-1][j-1]
        fBasis = eBasesCase1[j-1][i-1]
        fDualMatrix = eDualMatricesCase1[j-1][i-1]'''
    
    eBasis = result(pathSet)
    fBasis = result(list(reversed(pathSet)))
    eDualMatrix = dual(eBasis)
    fDualMatrix = dual(fBasis)

    eDual = 0 # e*
    fDual = 0 # f*

    for k in range(len(eBasis)):
        eDual += eDualMatrix[0][k] * indicesToF(eBasis[k]) if (np.abs(eDualMatrix[0][k]) > coeff_epsilon) else 0
    for k in range(len(fBasis)):
        fDual += fDualMatrix[0][k] * indicesToE(fBasis[k]) if (np.abs(fDualMatrix[0][k]) > coeff_epsilon) else 0
    return (eDual, fDual)


def leftsum():
    print("@@@@@@@@ Left sum: @@@@@@@@@@", file = f)
    sum = 0
    for i in range(1,n+1): # L_i
        secondterm = H[n-2] - H[n-1] # H_{n-1} - H_n
        for j in range(i,n):
            secondterm -= 2 * H[j-1]
        summand = q**(-2*n + 2*i) * (q**secondterm)
        sum += summand
        #print("+i: ", i, ", summand: ", summand)
        print("+i: ", i)
        print(f"+i: {i}\t summand: {summand}", file = f)

    for i in range(1,n+1): # -L_i
        secondterm = - H[n-2] + H[n-1] # H_{n-1} - H_n
        for j in range(i,n):
            secondterm += 2 * H[j-1]
        summand = q**(2*n - 2*i) * (q**secondterm)
        sum += summand
        #print("-i: ", i, ", summand: ", summand)
        print("-i: ", i)
        print(f"-i: {i} \t summand: {summand}", file = f)
    return sum

def rightsum():
    print("@@@@@@@@ Right sum: @@@@@@@@", file = f)
    sum = 0
  
    for i in range(1,n+1):   # Computing dual elements e* and f* CASE 1
        for j in range(i+1,n+1): # μ = L_i > λ = L_j (i < j)
            pathSet = getPathSet(n,i,j,1) # e_{μλ} = E_{pathSet} = E_{i, ..., j - 1}
            coeff = q**(1 - 2*n + 2*i) * (q - q**(-1))**(2 * len(pathSet))
            eDual, fDual = dualElements(pathSet, i, j, 1)
            qExponent = H[n-2] - H[n-1] # H_{-L_i - L_j} = H_{-μ - λ}
            for k in range(i,n):
                qExponent -= 2 * H[k-1]
            for k in range(i,j):
                qExponent += H[k-1]
            summand = ((coeff * eDual * (q**qExponent) * fDual))
            #print("i: ", i, ", j: ", j, ", summand: ", summand)
            print("i: ", i, ", j: ", j)
            print(f"i: {i} \t j: {j} \t summand: {summand}", file = f)
            sum += summand
    

    for i in range(1,n+1): # CASE 2
        for j in range(1,n+1): 
            if (i != n or j != n): # μ = L_i > λ = -L_j
                pathSet = getPathSet(n,i,j,2) # e_{μλ} = +- E_{pathSet} = +- E_{i, ..., n-2, n, n-1, ..., j} (if i < n) or +- E_{n, n-2, ..., j} (if i = n) or +- E_{i, ..., n-2, n} (if j = n)
                if (i == 1 and j == 1):
                    print(resultTemp(pathSet), file = p)

                #coeff = q**(2 + 2*n - 2*i) * (q - q**(-1))**(2 * len(pathSet)) if i == j else q**(1 + 2*n - 2*i) * (q - q**(-1))**(2 * len(pathSet))
                coeff = q**(2 - 2*n + 2*i) * (q - q**(-1))**(2 * len(pathSet)) if i == j else q**(1 - 2*n + 2*i) * (q - q**(-1))**(2 * len(pathSet))
                eDual, fDual = dualElements(pathSet, i, j, 2)
                qExponent = 0 # H_{-L_i + L_j} = H_{-μ - λ}
                for k in range(min(i,j), max(i,j)):
                    qExponent += H[k-1]
                if i < j:
                    qExponent *= -1
                #summand = simplify(((-1) * coeff * eDual * (q**qExponent) * fDual).subs(q**2 - 1, q*r).subs(q - 1/q, r))
                summand = (((-1) * coeff * eDual * (q**qExponent) * fDual))
                #print("i: ", i, ", -j: ", j, ", summand: ", summand)
                print("i: ", i, ", -j: ", j)
                print(f"i: {i} \t -j: {j} \t summand: {summand}", file = f)
                sum += summand

    for j in range(1,n+1): # CASE 3
        for i in range(j+1,n+1): # μ = -L_i > λ = -L_j (i > j)
            pathSet = getPathSet(n, i, j, 3) # e_{μλ} = E_{pathSet} = E_{i-1, ..., j}
            #coeff = q**(1 + 2*i - 2*n) * (q - q**(-1))**(2 * len(pathSet))
            coeff = q**(1 - 2*i + 2*n) * (q - q**(-1))**(2 * len(pathSet))
            eDual, fDual = dualElements(pathSet, i, j, 3)
            qExponent = - H[n-2] + H[n-1] # H_{L_i + L_j} = H_{-μ - λ}
            for k in range(i,n):
                qExponent += 2 * H[k-1]
            for k in range(j,i):
                qExponent += H[k-1]
            #summand = simplify((coeff * eDual * (q**qExponent) * fDual).subs(q**2 - 1, q*r).subs(q - 1/q, r))
            summand = ((coeff * eDual * (q**qExponent) * fDual))
            #print("-i: ", i, ", -j: ", j, ", summand: ", summand)
            print("-i: ", i, ", -j: ", j,)
            print(f"-i: {i} \t -j: {j} \t summand: {summand}", file = f)
            sum += summand 

    return sum

testBasis = resultTemp(getPathSet(5, 1, 5, 1))
#print((leftsum() + rightsum()))
#print(f"\n @@@@@@@@ Total Sum @@@@@@@@{(leftsum() + rightsum())}", file = f)
'''print(f"n = {n}, q = {q}, cutoff = {perm_epsilon}", file=f)
resultTemp([1,2,4,3,2,1])'''

#print(result([1,3,2]))

'''
testBasis11 = [[1, 2, 3, 5, 4, 3, 2, 1], [1, 2, 3, 5, 4, 3, 1, 2], [1, 2, 3, 5, 4, 2, 3, 1], [1, 2, 3, 5, 4, 1, 2, 3], [1, 2, 3, 5, 3, 4, 2, 1], [1, 2, 3, 5, 3, 4, 1, 2], [1, 2, 3, 5, 2, 3, 4, 1], [1, 2, 3, 5, 1, 2, 3, 4], [1, 2, 3, 4, 3, 5, 2, 1], [1, 2, 3, 4, 3, 5, 1, 2], [1, 2, 3, 4, 2, 3, 5, 1], [1, 2, 3, 4, 1, 2, 3, 5], [1, 2, 3, 3, 5, 4, 2, 1], [1, 2, 3, 3, 5, 4, 1, 2], [1, 2, 3, 2, 3, 5, 4, 1], [1, 2, 3, 1, 2, 3, 5, 4], [1, 2, 5, 3, 4, 3, 2, 1], [1, 2, 5, 3, 4, 3, 1, 2], [1, 2, 5, 3, 4, 2, 3, 1], [1, 2, 5, 3, 4, 1, 2, 3], [1, 2, 5, 3, 2, 3, 4, 1], [1, 2, 5, 3, 1, 2, 3, 4], [1, 2, 5, 4, 3, 2, 3, 1], [1, 2, 5, 4, 3, 1, 2, 3], [1, 2, 4, 3, 5, 2, 3, 1], [1, 2, 4, 3, 5, 1, 2, 3], [1, 2, 4, 3, 2, 3, 5, 1], [1, 2, 4, 3, 1, 2, 3, 5], [1, 2, 2, 3, 5, 4, 3, 1], [1, 2, 1, 2, 3, 5, 4, 3], [1, 3, 2, 5, 4, 3, 2, 1], [1, 3, 2, 5, 4, 3, 1, 2], [1, 3, 2, 5, 4, 1, 2, 3], [1, 3, 2, 5, 3, 4, 2, 1], [1, 3, 2, 5, 3, 4, 1, 2], [1, 3, 2, 5, 1, 2, 3, 4], [1, 3, 2, 4, 3, 5, 2, 1], [1, 3, 2, 4, 3, 5, 1, 2], [1, 3, 2, 4, 1, 2, 3, 5], [1, 3, 2, 3, 5, 4, 2, 1], [1, 3, 2, 3, 5, 4, 1, 2], [1, 3, 2, 1, 2, 3, 5, 4], [1, 3, 5, 4, 3, 2, 1, 2], [1, 3, 5, 3, 2, 4, 1, 2], [1, 3, 4, 3, 2, 5, 1, 2], [1, 5, 3, 2, 4, 3, 2, 1], [1, 5, 3, 2, 4, 3, 1, 2], [1, 5, 3, 2, 4, 1, 2, 3], [1, 5, 3, 2, 1, 2, 3, 4], [1, 5, 3, 4, 3, 2, 1, 2], [1, 5, 3, 3, 4, 1, 2, 2], [1, 5, 4, 3, 2, 1, 2, 3], [1, 4, 3, 2, 5, 3, 1, 2], [1, 4, 3, 2, 5, 1, 2, 3], [1, 4, 3, 2, 1, 2, 3, 5], [1, 1, 2, 3, 5, 4, 3, 2], [1, 1, 2, 2, 3, 3, 5, 4], [2, 1, 3, 5, 4, 3, 2, 1], [2, 1, 3, 5, 4, 2, 3, 1], [2, 1, 3, 5, 3, 4, 2, 1], [2, 1, 3, 5, 2, 3, 4, 1], [2, 1, 3, 4, 3, 5, 2, 1], [2, 1, 3, 4, 2, 3, 5, 1], [2, 1, 3, 3, 5, 4, 2, 1], [2, 1, 3, 2, 3, 5, 4, 1], [2, 1, 5, 3, 4, 3, 2, 1], [2, 1, 5, 3, 4, 2, 3, 1], [2, 1, 5, 3, 2, 3, 4, 1], [2, 1, 5, 4, 3, 2, 3, 1], [2, 1, 4, 3, 5, 2, 3, 1], [2, 1, 4, 3, 2, 3, 5, 1], [2, 1, 2, 3, 5, 4, 3, 1], [2, 5, 4, 2, 1, 3, 3, 1], [3, 2, 1, 5, 4, 3, 2, 1], [3, 2, 1, 5, 3, 4, 2, 1], [3, 2, 1, 4, 3, 5, 2, 1], [3, 2, 1, 3, 5, 4, 2, 1], [5, 3, 2, 1, 4, 3, 2, 1], [5, 3, 3, 2, 4, 2, 1, 1], [5, 4, 3, 3, 2, 2, 1, 1]]
testBasis12 = [[1, 2, 3, 5, 4, 3, 2], [1, 2, 3, 5, 4, 2, 3], [1, 2, 3, 5, 3, 4, 2], [1, 2, 3, 5, 2, 3, 4], [1, 2, 3, 4, 3, 5, 2], [1, 2, 3, 4, 2, 3, 5], [1, 2, 3, 3, 5, 4, 2], [1, 2, 3, 2, 3, 5, 4], [1, 2, 5, 3, 4, 3, 2], [1, 2, 5, 3, 4, 2, 3], [1, 2, 5, 3, 2, 3, 4], [1, 2, 5, 4, 3, 2, 3], [1, 2, 4, 3, 5, 2, 3], [1, 2, 4, 3, 2, 3, 5], [1, 2, 2, 3, 5, 4, 3], [1, 3, 2, 5, 4, 3, 2], [1, 3, 2, 5, 3, 4, 2], [1, 3, 2, 4, 3, 5, 2], [1, 3, 2, 3, 5, 4, 2], [1, 5, 3, 2, 4, 3, 2], [2, 1, 3, 5, 4, 3, 2], [2, 1, 3, 5, 4, 2, 3], [2, 1, 3, 5, 3, 4, 2], [2, 1, 3, 5, 2, 3, 4], [2, 1, 3, 4, 3, 5, 2], [2, 1, 3, 4, 2, 3, 5], [2, 1, 3, 3, 5, 4, 2], [2, 1, 3, 2, 3, 5, 4], [2, 1, 5, 3, 4, 3, 2], [2, 1, 5, 3, 4, 2, 3], [2, 1, 5, 3, 2, 3, 4], [2, 1, 5, 4, 3, 2, 3], [2, 1, 4, 3, 5, 2, 3], [2, 1, 4, 3, 2, 3, 5], [2, 1, 2, 3, 5, 4, 3], [2, 3, 5, 4, 3, 2, 1], [2, 3, 5, 4, 2, 1, 3], [2, 3, 5, 3, 4, 2, 1], [2, 3, 5, 2, 1, 3, 4], [2, 3, 4, 3, 5, 2, 1], [2, 3, 4, 2, 1, 3, 5], [2, 3, 3, 5, 4, 2, 1], [2, 5, 3, 4, 3, 2, 1], [2, 5, 3, 4, 2, 1, 3], [2, 4, 3, 5, 2, 1, 3], [3, 2, 1, 5, 4, 3, 2], [3, 2, 1, 5, 3, 4, 2], [3, 2, 1, 4, 3, 5, 2], [3, 2, 1, 3, 5, 4, 2], [3, 2, 5, 4, 3, 2, 1], [3, 2, 5, 3, 4, 2, 1], [3, 2, 4, 3, 5, 2, 1], [5, 3, 2, 1, 4, 3, 2], [5, 3, 2, 4, 3, 2, 1], [4, 3, 2, 1, 5, 3, 2]]
testBasis13 = [[1, 2, 3, 5, 4, 3], [1, 2, 3, 5, 3, 4], [1, 2, 3, 4, 3, 5], [1, 2, 3, 3, 5, 4], [1, 2, 5, 3, 4, 3], [1, 3, 2, 5, 4, 3], [1, 3, 2, 5, 3, 4], [1, 3, 2, 4, 3, 5], [1, 3, 2, 3, 5, 4], [1, 3, 5, 4, 3, 2], [1, 3, 5, 3, 2, 4], [1, 3, 4, 3, 2, 5], [1, 5, 3, 2, 4, 3], [1, 5, 3, 4, 3, 2], [1, 4, 3, 2, 5, 3], [2, 1, 3, 5, 4, 3], [2, 1, 3, 5, 3, 4], [2, 1, 3, 4, 3, 5], [2, 1, 3, 3, 5, 4], [2, 1, 5, 3, 4, 3], [3, 2, 1, 5, 4, 3], [3, 2, 1, 5, 3, 4], [3, 2, 1, 4, 3, 5], [3, 2, 1, 3, 5, 4], [3, 5, 4, 3, 2, 1], [3, 5, 3, 2, 1, 4], [3, 4, 3, 2, 1, 5], [5, 3, 2, 1, 4, 3], [5, 3, 4, 3, 2, 1], [4, 3, 2, 1, 5, 3]]
testBasis14 = [[1, 2, 3, 5, 4], [1, 2, 5, 3, 4], [1, 2, 5, 4, 3], [1, 2, 4, 3, 5], [1, 3, 2, 5, 4], [1, 5, 3, 2, 4], [1, 5, 4, 3, 2], [1, 4, 3, 2, 5], [2, 1, 3, 5, 4], [2, 1, 5, 3, 4], [2, 1, 5, 4, 3], [2, 1, 4, 3, 5], [3, 2, 1, 5, 4], [5, 3, 2, 1, 4], [5, 4, 3, 2, 1], [4, 3, 2, 1, 5]]
testBasis15 = [[1, 2, 3, 4], [1, 2, 4, 3], [1, 3, 2, 4], [1, 4, 3, 2], [2, 1, 3, 4], [2, 1, 4, 3], [3, 2, 1, 4], [4, 3, 2, 1]]
'''

print(len(testBasis))

fin = []
count = 1
for i in range(len(testBasis)):
    fin.append(testBasis[i])
    M = mat(fin)
    print(f"{count} det M = {la.det(M)}\n", file=f)
    count += 1

eDualMatrix = dual(testBasis)

eDual = 0 # e*

for k in range(len(testBasis)):
    eDual += eDualMatrix[0][k] * indicesToF(testBasis[k]) * (q - 1/q)**3 if (np.abs(eDualMatrix[0][k]) > coeff_epsilon) else 0
    print(eDualMatrix[0][k] * indicesToF(testBasis[k]), file=f)

print(eDual, file=f)

#print(result(getPathSet(4, 1, 1, 2)))


