from sympy import *
import itertools
from collections import deque
f = open('so6symb.txt', 'w')

n = 4
var('q')
var('r')
#q = 10
H = symbols('H(1:5)', commutative = False) # These are 0-indexed (so H[0] is H1)
E = symbols('E(1:5)', commutative = False)
F = symbols('F(1:5)', commutative = False)
eBasesCase2 = [[0 for j in range(n)] for i in range(n)] # Stores eBasis for different values of i and j (eBasis for i, j is the same as fBasis for j, i)
eDualMatricesCase2 = [[0 for j in range(n)] for i in range(n)]
eBasesCase1 = [[0 for j in range(n)] for i in range(n)] 
eDualMatricesCase1 = [[0 for j in range(n)] for i in range(n)]

print(f"q = {q} \t n = {n}", file = f)

def inverse(M):
    detM = M.det()
    size = shape(M)[0]
    coM = zeros(size)
    signRow = 1
    for i in range(size):
        signCol = signRow       #
        for j in range(size):
            minor = M[:,:]
            minor.col_del(j)
            minor.row_del(i)
            coM[i,j] = signCol * minor.det()
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
    N= inverse(M)
    return N.applyfunc(simplify)

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

# Takes in a path set for E and returns (e*, f*)
def dualElements(pathSet, i, j, case):
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

    
    #fBasis = result(list(reversed(pathSet)))
    #eDualMatrix = dual(eBasis)
    #print("eDualCoeffs: ", eDualCoeffs)
    #fDualMatrix = dual(fBasis) 
    eDual = 0 # e*
    fDual = 0 # f*
    for k in range(len(eBasis)):
        eDual += eDualMatrix[0,k] * indicesToF(eBasis[k])
        fDual += fDualMatrix[0,k] * indicesToE(fBasis[k])
    #return (eDual, fDual)
    return (eDual, fDual)

def leftsum():
    print("@@@@@@@@ Left sum: @@@@@@@@@@", file = f)
    sum = 0
    for i in range(1,n+1): # L_i
        secondterm = H[n-2] - H[n-1] # H_{n-1} - H_n
        for j in range(i,n):
            secondterm -= 2 * H[j-1]
        #summand = q**(2*n - 2*i) * (q**secondterm)
        summand = q**(-2*n + 2*i) * (q**secondterm)
        sum += summand
        #print("+i: ", i, ", summand: ", summand)
        print("+i: ", i)
        print(f"+i: {i}\t summand: {summand}", file = f)

    for i in range(1,n+1): # -L_i
        secondterm = - H[n-2] + H[n-1] # H_{n-1} - H_n
        for j in range(i,n):
            secondterm += 2 * H[j-1]
        #summand = q**(2*i - 2*n) * (q**secondterm)
        summand = q**(2*n - 2*i) * (q**secondterm)
        sum += summand
        #print("-i: ", i, ", summand: ", summand)
        print("-i: ", i)
        print(f"-i: {i} \t summand: {summand}", file = f)
    return sum

def rightsum():
    print("@@@@@@@@ Right sum: @@@@@@@@@@", file = f)
    sum = 0
  
    for i in range(1,n+1):   # Computing dual elements e* and f* CASE 1
        for j in range(i+1,n+1): # μ = L_i > λ = L_j (i < j)
            pathSet = [x for x in range(i, j)] # e_{μλ} = E_{pathSet} = E_{i, ..., j - 1}
            #coeff = q**(1 + 2*n - 2*i) * (q - q**(-1))**(2 * len(pathSet))
            coeff = q**(1 - 2*n + 2*i) * (q - q**(-1))**(2 * len(pathSet))
            eDual, fDual = dualElements(pathSet, i, j, 1)
            qExponent = H[n-2] - H[n-1] # H_{-L_i - L_j} = H_{-μ - λ}
            for k in range(i,n):
                qExponent -= 2 * H[k-1]
            for k in range(i,j):
                qExponent += H[k-1]
            summand = simplify((coeff * eDual * (q**qExponent) * fDual).subs(q**2 - 1, q*r).subs(q - 1/q, r))
            #summand = ((coeff * eDual * (q**qExponent) * fDual))
            #print("i: ", i, ", j: ", j, ", summand: ", summand)
            print("i: ", i, ", j: ", j)
            print(f"i: {i} \t j: {j} \t summand: {summand}", file = f)
            sum += summand
    
    for j in range(1,n+1): # CASE 3
        for i in range(j+1,n+1): # μ = -L_i > λ = -L_j (i > j)
            pathSet = [x for x in reversed(range(j, i))] # e_{μλ} = E_{pathSet} = E_{i-1, ..., j}
            #coeff = q**(1 + 2*i - 2*n) * (q - q**(-1))**(2 * len(pathSet))
            coeff = q**(1 - 2*i + 2*n) * (q - q**(-1))**(2 * len(pathSet))
            eDual, fDual = dualElements(pathSet, i, j, 3)
            qExponent = - H[n-2] + H[n-1] # H_{L_i + L_j} = H_{-μ - λ}
            for k in range(i,n):
                qExponent += 2 * H[k-1]
            for k in range(j,i):
                qExponent += H[k-1]
            summand = simplify((coeff * eDual * (q**qExponent) * fDual).subs(q**2 - 1, q*r).subs(q - 1/q, r))
            #summand = ((coeff * eDual * (q**qExponent) * fDual))
            #print("-i: ", i, ", -j: ", j, ", summand: ", summand)
            print("-i: ", i, ", -j: ", j,)
            print(f"-i: {i} \t -j: {j} \t summand: {summand}", file = f)
            sum += summand 


    for i in reversed(range(1,n+1)): # CASE 2
        for j in reversed(range(1,n+1)): 
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
                qExponent = 0 # H_{-L_i + L_j} = H_{-μ - λ}
                for k in range(min(i,j), max(i,j)):
                    qExponent += H[k-1]
                if i < j:
                    qExponent *= -1
                summand = simplify(((-1) * coeff * eDual * (q**qExponent) * fDual).subs(q**2 - 1, q*r).subs(q - 1/q, r))
                #summand = (((-1) * coeff * eDual * (q**qExponent) * fDual))
                #print("i: ", i, ", -j: ", j, ", summand: ", summand)
                print("i: ", i, ", -j: ", j)
                print(f"i: {i} \t -j: {j} \t summand: {summand}", file = f)
                sum += summand

    

    return sum

#print((leftsum() + rightsum()).subs(q-1/q,r))
#print(f"\n @@@@@@@@ Total Sum @@@@@@@@{(leftsum() + rightsum())}", file = f)
print(result([2,4,3]))


