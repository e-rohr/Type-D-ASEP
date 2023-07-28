from sympy import *
import itertools
from collections import deque





n = 3
q = Symbol('q')
r = Symbol('r')
K = []
KNeg = []
H = []
F = []
E = []

for i in range(2*n):
    K.append(MatrixSymbol(f'K{i+1}',2*n, 2*n))
    KNeg.append(MatrixSymbol(f'KNeg{i+1}',2*n, 2*n))
    H.append(MatrixSymbol(f'H{i+1}',2*n, 2*n))
    E.append(MatrixSymbol(f'E{i+1}',2*n, 2*n))
    F.append(MatrixSymbol(f'F{i+1}',2*n, 2*n))

eBasesCase2 = [[0 for j in range(n)] for i in range(n)] # Stores eBasis for different values of i and j (eBasis for i, j is the same as fBasis for j, i)
eDualMatricesCase2 = [[0 for j in range(n)] for i in range(n)]
eBasesCase1 = [[0 for j in range(n)] for i in range(n)] 
eDualMatricesCase1 = [[0 for j in range(n)] for i in range(n)]





########### Representations of generators ###############

# 1-indexed, returns K_index in the fundamental representation of so(2n)
def k(index):
    K = zeros(2*n)
    H = h(index)
    for i in range(2*n):
        K[i, i] = q**(H[i,i])
    return K


# 1-indexed, returns 1/K_index in the fundamental representation of so(2n)
def kNeg(index):
    K = zeros(2*n)
    H = h(index)
    for i in range(2*n):
        K[i, i] = q**(-H[i,i])
    return K


# 1-indexed, returns K_index in the fundamental representation of so(2n)
def h(index):
    H = zeros(2*n)
    if (index < n):
        H[index - 1, index - 1] = 1
        H[index, index] = -1
        H[n + index - 1, n + index - 1] = -1
        H[n + index, n + index] = 1
    else:
        H[n - 2, n - 2] = 1
        H[n - 1, n - 1] = 1
        H[2*n - 2, 2*n - 2] = -1
        H[2*n - 1, 2*n - 1] = -1
        
    return H


# 1-indexed, returns E_index in the fundamental representation of so(2n)
def e(index):
    E = zeros(2*n)
    if (index < n):
        E[index - 1, index] = 1
        E[n + index, n + index - 1] = -1
    else:
        E[n - 2, 2*n - 1] = 1
        E[n - 1, 2*n - 2] = -1
    return E


# 1-indexed, returns F_index in the fundamental representation of so(2n)
def f(index):
    F = zeros(2*n)
    if (index < n):
        F[index, index - 1] = 1
        F[n - 1 + index, n + index] = -1
    else:
        F[2*n - 2, n - 1] = 1
        F[2*n - 1, n - 2] = - 1
    return F





######## Matrix Inverse Using Adjugates ########

def inverse(M):
    detM = M.det(method = 'domain-ge')
    size = len(M.row(0))
    coM = zeros(size)
    signRow = 1
    for i in range(size):
        signCol = signRow       
        for j in range(size):
            minor = M[:,:]
            minor.col_del(j)
            minor.row_del(i)
            coM[i,j] = signCol * minor.det(method = 'domain-ge')
            signCol *= -1
        signRow *= -1
    adjM = coM.T
    return (1/detM)*adjM      





######## Andrew Lin Code ########

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
    return M.applyfunc(simplify)


def dual(listofLists): # for computing specific dual elements
    M = zeros(len(listofLists))
    for i in range(len(listofLists)):
        for j in range(len(listofLists)):
            lst = pair(listofLists[i], listofLists[j])
            for k in lst:
                M[i, j] += q**k
    N = inverse(M)
    return N.applyfunc(simplify)


def perm(tentlist): # given list of lists, removes linear dependence
    fin = []
    for i in range(len(tentlist)):
        fin1 = list(fin)
        fin1.append(tentlist[i])
        M = mat(fin1)
        d = M.det(method = 'domain-ge') # Modified from original code to use Gaussian Elimination for better performance
        if(d != 0):
            fin.append(tentlist[i])
    return fin


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


'''
The first line of result was changed from the original code so that 
the first row of the matrix returned by dual() would correspond to the desired dual element.
'''
def result(setofindices): # gives a basis of elements (the first element is the basis of setofindices)
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





########## Dual Element Code ########## 

# Takes in a set of indices and returns the corresponding product of E_i's
# Ex: [1, 2, 3] -> E1 E2 E3
def indicesToE(setofindices):
    solution = 1
    for index in setofindices:
        solution = solution * E[index-1]
    return solution


# Takes in a set of indices and returns the corresponding product of F_i's
# Ex: [1, 2, 3] -> F1 F2 F3
def indicesToF(setofindices):
    solution = 1
    for index in setofindices:
        solution = solution * F[index-1]
    return solution


# Takes in a path set for E and returns (e*, f*)
def dualElements(pathSet, i, j, case):
    # Calculate the basis and q-pairings matrix for e = E_{pathSet} and f = F_{pathSet}
    # Store these in eBasesCase1, eBasesCase2, eDualMatricesCase1, and eDualMatricesCase2 to avoid redundant computations
    if (case == 1):
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
        fDualMatrix = eDualMatricesCase1[j-1][i-1]

    # Calculate e*, f*
    eDual = ZeroMatrix(2*n, 2*n) # e*
    fDual = ZeroMatrix(2*n, 2*n) # f*
    for k in range(len(eBasis)):
        eDual += eDualMatrix[0,k] * indicesToF(eBasis[k])
        fDual += fDualMatrix[0,k] * indicesToE(fBasis[k])
    return (eDual, fDual)





############ Central Element ##############

# \sum_{\mu} q^{(-2\rho, \mu) q^{H_{-2\mu}}}

def sum():
    sum = ZeroMatrix(2*n, 2*n)
    
    # CASE 0: μ = λ
    for i in range(1,n+1): # μ = L_i
        secondterm = K[n-2] * KNeg[n-1]
        for j in range(i,n):
            secondterm *= KNeg[j-1]**2
        summand = q**(-2*n + 2*i) * secondterm
        sum += summand

    for i in range(1,n+1): # μ = -L_i
        secondterm = KNeg[n-2] * K[n-1]
        for j in range(i,n):
            secondterm *= K[j-1]**2
        summand = q**(2*n - 2*i) * secondterm
        sum += summand

    # CASE 1: μ = L_i > λ = L_j (i < j)
    for i in range(1,n+1):   
        for j in range(i+1,n+1): 
            pathSet = [x for x in range(i, j)] # e_{μλ} = E_{pathSet} = E_{i, ..., j - 1}
            coeff = q**(1 - 2*n + 2*i) * (q - q**(-1))**(2 * len(pathSet))
            eDual, fDual = dualElements(pathSet, i, j, 1)
            kTerm = K[n-2] * KNeg[n-1] # q**(H_{-L_i - L_j}) = q**H_{-μ - λ}
            for k in range(i,n):
                kTerm *= KNeg[k-1]**2
            for k in range(i,j):
                kTerm *= K[k-1]
            summand = ((coeff * eDual * kTerm * fDual))
            sum += summand
    
    # CASE 2: μ = L_i > λ = -L_j
    for i in range(1,n+1): 
        for j in range(1,n+1): 
            if (i != n or j != n):  
                pathSet = [] # e_{μλ} = +- E_{pathSet} = +- E_{i, ..., n-2, n, n-1, ..., j} (if i < n) or +- E_{n, n-2, ..., j} (if i = n) or +- E_{i, ..., n-2, n} (if j = n)
                for x in range(i,n-1):
                    pathSet.append(x)
                pathSet.append(n)
                if (i != n and j != n):
                    pathSet.append(n-1)
                for x in reversed(range(j,n-1)):
                    pathSet.append(x) 
                coeff = q**(2 - 2*n + 2*i) * (q - q**(-1))**(2 * len(pathSet)) if i == j else q**(1 - 2*n + 2*i) * (q - q**(-1))**(2 * len(pathSet))
                eDual, fDual = dualElements(pathSet, i, j, 2)
                kTerm = 1 # q**(H_{-L_i + L_j}) = q**H_{-μ - λ}
                if i < j:
                    for k in range(i, j):
                        kTerm *= KNeg[k-1]
                else:
                    for k in range(j, i):
                        kTerm *= K[k-1]
                summand = (((-1) * coeff * eDual * kTerm * fDual))
                sum += summand

    # CASE 3: μ = -L_i > λ = -L_j (i > j)
    for j in range(1,n+1): 
        for i in range(j+1,n+1): 
            pathSet = [x for x in reversed(range(j, i))] # e_{μλ} = E_{pathSet} = E_{i-1, ..., j}
            coeff = q**(1 - 2*i + 2*n) * (q - q**(-1))**(2 * len(pathSet))
            eDual, fDual = dualElements(pathSet, i, j, 3)
            kTerm = KNeg[n-2] * K[n-1] # q**(H_{L_i + L_j}) = q**H_{-μ - λ}
            for k in range(i,n):
                kTerm *= K[k-1]**2
            for k in range(j,i):
                kTerm *= K[k-1]
            summand = ((coeff * eDual * kTerm * fDual))
            sum += summand 

    return sum

print(sum())