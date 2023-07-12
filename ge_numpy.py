from sympy import *
import itertools
from collections import deque
import numpy as np
import numpy.linalg as la



n = 5
fout = open(f'so{2*n}_ge_num.txt', 'w')
#q = Symbol('q')
#r = Symbol('r')
q = 10
perm_epsilon = 0.00000001
coeff_epsilon = q**(-2*n)

'''K = []
KNeg = []
H = []
F = []
E = []
for i in range(2*n):
    K.append(MatrixSymbol(f'K{i+1}',2*n, 2*n))
    KNeg.append(MatrixSymbol(f'KNeg{i+1}',2*n, 2*n))
    H.append(MatrixSymbol(f'H{i+1}',2*n, 2*n))
    E.append(MatrixSymbol(f'E{i+1}',2*n, 2*n))
    F.append(MatrixSymbol(f'F{i+1}',2*n, 2*n))'''

H = symbols('H(1:6)', commutative = False) # These are 0-indexed (so H[0] is H1)
E = symbols('E(1:6)', commutative = False)
F = symbols('F(1:6)', commutative = False)
K = symbols('K(1:6)', commutative = False)
KNeg = symbols('KNeg(1:6)', commutative = False)
'''eBasesCase2 = [[0 for j in range(n)] for i in range(n)] # Stores eBasis for different values of i and j (eBasis for i, j is the same as fBasis for j, i)
eDualMatricesCase2 = [[0 for j in range(n)] for i in range(n)]
eBasesCase1 = [[0 for j in range(n)] for i in range(n)] 
eDualMatricesCase1 = [[0 for j in range(n)] for i in range(n)]'''
print(f"q = {q} \t n = {n}", file = fout)















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
    count = 1
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
            #print(f'{count} \n')
            count += 1
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

def perm(tentlist): # given list of lists, removes linear dependence
    fin = []
    count = 1
    for i in range(len(tentlist)):
        fin1 = list(fin)
        fin1.append(tentlist[i])
        M = mat(fin1)
        
        if (np.abs(la.det(M)) > perm_epsilon):
            fin.append(tentlist[i])  

        print(count)
        count += 1     
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




########## Dual Element Code ########## 

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

    #eDual = ZeroMatrix(2*n, 2*n) # e*
    #fDual = ZeroMatrix(2*n, 2*n) # f*
    eDual = 0
    fDual = 0
    for k in range(len(eBasis)):
        eDual += eDualMatrix[0][k] * indicesToF(eBasis[k]) if (np.abs(eDualMatrix[0][k]) > coeff_epsilon) else 0
    for k in range(len(fBasis)):
        fDual += fDualMatrix[0][k] * indicesToE(fBasis[k]) if (np.abs(fDualMatrix[0][k]) > coeff_epsilon) else 0
    return (eDual, fDual)




############ Central Element ##############

def leftsum():
    print("@@@@@@@@ Left sum: @@@@@@@@@@", file = fout)
    #sum = ZeroMatrix(2*n, 2*n)
    sum = 0
    for i in range(1,n+1): # L_i
        #secondterm = H[n-2] - H[n-1] # H_{n-1} - H_n
        secondterm = K[n-2] * KNeg[n-1]
        for j in range(i,n):
            #secondterm -= 2 * H[j-1]
            secondterm *= KNeg[j-1]**2
        #summand = q**(2*n - 2*i) * (q**secondterm)
        #summand = q**(-2*n + 2*i) * (q**secondterm)
        summand = q**(-2*n + 2*i) * secondterm
        sum += summand
        #print("+i: ", i, ", summand: ", summand)
        print("+i: ", i)
        print(f"+i: {i}\t summand: {summand}", file = fout)

    for i in range(1,n+1): # -L_i
        #secondterm = - H[n-2] + H[n-1] # H_{n-1} - H_n
        secondterm = KNeg[n-2] * K[n-1]
        for j in range(i,n):
            #secondterm += 2 * H[j-1]
            secondterm *= K[j-1]**2
        #summand = q**(2*n - 2*i) * (q**secondterm)
        summand = q**(2*n - 2*i) * secondterm
        sum += summand
        #print("-i: ", i, ", summand: ", summand)
        print("-i: ", i)
        print(f"-i: {i} \t summand: {summand}", file = fout)
    return sum

def rightsum():
    print("@@@@@@@@ Right sum: @@@@@@@@@@", file = fout)
    #sum = ZeroMatrix(2*n, 2*n)
    sum = 0
  
    for i in range(1,n+1):   # Computing dual elements e* and f* CASE 1
        for j in range(i+1,n+1): # μ = L_i > λ = L_j (i < j)
            pathSet = [x for x in range(i, j)] # e_{μλ} = E_{pathSet} = E_{i, ..., j - 1}
            coeff = q**(1 - 2*n + 2*i) * (q - q**(-1))**(2 * len(pathSet))
            eDual, fDual = dualElements(pathSet, i, j, 1)
            
            #qExponent = H[n-2] - H[n-1] # H_{-L_i - L_j} = H_{-μ - λ}
            kTerm = K[n-2] * KNeg[n-1] # q**(H_{-L_i - L_j}) = q**H_{-μ - λ}
            for k in range(i,n):
                #qExponent -= 2 * H[k-1]
                kTerm *= KNeg[k-1]**2
            for k in range(i,j):
                #qExponent += H[k-1]
                kTerm *= K[k-1]
            print(coeff, eDual, kTerm, fDual)
            #summand = simplify((coeff * eDual * (q**qExponent) * fDual).subs(q**2 - 1, q*r).subs(q - 1/q, r))
            #summand = ((coeff * eDual * (q**qExponent) * fDual))
            summand = ((coeff * eDual * kTerm * fDual))
            #print("i: ", i, ", j: ", j, ", summand: ", summand)
            print("i: ", i, ", j: ", j)
            print(f"i: {i} \t j: {j} \t summand: {summand}", file = fout)
            sum += summand
    

    for i in range(1,n+1): # CASE 2
        for j in range(1,n+1): 
            #if (i == 1 and j == 1):
            #   continue
            if (i != n or j != n): # μ = L_i > λ = -L_j
                pathSet = [] # e_{μλ} = +- E_{pathSet} = +- E_{i, ..., n-2, n, n-1, ..., j} (if i < n) or +- E_{n, n-2, ..., j} (if i = n) or +- E_{i, ..., n-2, n} (if j = n)
                for x in range(i,n-1):
                    pathSet.append(x)
                pathSet.append(n)
                if (i != n and j != n):
                    pathSet.append(n-1)
                for x in reversed(range(j,n-1)):
                    pathSet.append(x) 
                #coeff = q**(2 + 2*n - 2*i) * (q - q**(-1))**(2 * len(pathSet)) if i == j else q**(1 + 2*n - 2*i) * (q - q**(-1))**(2 * len(pathSet))
                coeff = q**(2 - 2*n + 2*i) * (q - q**(-1))**(2 * len(pathSet)) if i == j else q**(1 - 2*n + 2*i) * (q - q**(-1))**(2 * len(pathSet))
                eDual, fDual = dualElements(pathSet, i, j, 2)
                #qExponent = 0 # H_{-L_i + L_j} = H_{-μ - λ}
                kTerm = 1 # q**(H_{-L_i + L_j}) = q**H_{-μ - λ}
                if i < j:
                    for k in range(i, j):
                        kTerm *= KNeg[k-1]
                else:
                    for k in range(j, i):
                        kTerm *= K[k-1]
                #for k in range(min(i,j), max(i,j)):
                    #qExponent += H[k-1]
                    #kTerm *= K[k-1] 
                #if i < j:
                    #qExponent *= -1
                    #kTerm *= -1
                #summand = simplify(((-1) * coeff * eDual * (q**qExponent) * fDual).subs(q**2 - 1, q*r).subs(q - 1/q, r))
                #summand = (((-1) * coeff * eDual * (q**qExponent) * fDual))
                summand = (((-1) * coeff * eDual * kTerm * fDual))
                #print("i: ", i, ", -j: ", j, ", summand: ", summand)
                print("i: ", i, ", -j: ", j)
                print(f"i: {i} \t -j: {j} \t summand: {summand}", file = fout)
                sum += summand

    for j in range(1,n+1): # CASE 3
        for i in range(j+1,n+1): # μ = -L_i > λ = -L_j (i > j)
            pathSet = [x for x in reversed(range(j, i))] # e_{μλ} = E_{pathSet} = E_{i-1, ..., j}
            #coeff = q**(1 + 2*i - 2*n) * (q - q**(-1))**(2 * len(pathSet))
            coeff = q**(1 - 2*i + 2*n) * (q - q**(-1))**(2 * len(pathSet))
            eDual, fDual = dualElements(pathSet, i, j, 3)
            #qExponent = - H[n-2] + H[n-1] # H_{L_i + L_j} = H_{-μ - λ}
            kTerm = KNeg[n-2] * K[n-1] # q**(H_{L_i + L_j}) = q**H_{-μ - λ}
            for k in range(i,n):
                #qExponent += 2 * H[k-1]
                kTerm *= K[k-1]**2
            for k in range(j,i):
                #qExponent += H[k-1]
                kTerm *= K[k-1]
            #summand = simplify((coeff * eDual * (q**qExponent) * fDual).subs(q**2 - 1, q*r).subs(q - 1/q, r))
            #summand = ((coeff * eDual * (q**qExponent) * fDual))
            summand = ((coeff * eDual * kTerm * fDual))
            print("-i: ", i, ", -j: ", j,)
            print(f"-i: {i} \t -j: {j} \t summand: {summand}", file = fout)
            sum += summand

    return sum




################# Testing and Printing ###########################


#print((leftsum() + rightsum()).subs(q-1/q,r))
centralElement = leftsum() + rightsum()

#print(f"\n @@@@@@@@ Total Sum @@@@@@@@{(centralElement)}", file = fout)
#print(result(getPathSet(4,1,1,2)), file = fout)

subList = [(E[i], e(i + 1)) for i in range(n)]
subList += [(F[i], f(i + 1)) for i in range(n)]
subList += [(K[i], k(i + 1)) for i in range(n)]
subList += [(KNeg[i], kNeg(i + 1)) for i in range(n)]


#print(subList)
#print(Ks)

print(centralElement)

'''representation = centralElement.subs(subList).doit()
print(representation)
for i in range(2*n):
    print(representation[i,i])'''

'''C1  = [[1, 2, 4, 3, 2, 1], [1, 2, 4, 3, 1, 2], [1, 2, 4, 2, 3, 1], [1, 2, 4, 1, 2, 3], [1, 2, 3, 2, 4, 1], [1, 2, 3, 1, 2, 4], [1, 2, 2, 4, 3, 1], [1, 2, 1, 2, 4, 3], [1, 4, 2, 3, 2, 1], [1, 4, 2, 3, 1, 2], [1, 4, 2, 1, 2, 3], [1, 4, 3, 2, 1, 2], [1, 3, 2, 4, 1, 2], [1, 3, 2, 1, 2, 4], [1, 1, 2, 4, 3, 2], [2, 1, 4, 3, 2, 1], [2, 1, 4, 2, 3, 1], [2, 1, 3, 2, 4, 1], [2, 1, 2, 4, 3, 1], [4, 2, 1, 3, 2, 1]]



#C1 = getPathSet(4,1,4,1)
#C1 = result(C1)
m = mat(C1)
I = inverse(m)

print(f"m = {m}")
for i in range(len(C1)):
    for j in range(len(C1)):
         I[i,j] = factor(collect(powsimp(powdenest(expand(I[i,j]))), q))
print(f"I = {I}")

#I = m.inv(method  = 'ADJ')

B = I*m
for i in range(len(C1)):
    for j in range(len(C1)):
         B[i,j] = factor(collect(powsimp(powdenest(expand(B[i,j]))), q))
print(B)'''