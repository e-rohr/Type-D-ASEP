from sympy import *
import itertools
from collections import deque

n = 3
var('q')
var('r')
#q = 100
H = symbols('H(1:4)', commutative = False) # These are 0-indexed (so H[0] is H1)
E = symbols('E(1:4)', commutative = False)
F = symbols('F(1:4)', commutative = False)

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
    N= M.inv()
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
def dualElements(pathSet):
    eBasis = result(pathSet)
    fBasis = result(list(reversed(pathSet)))
    eDualMatrix = dual(eBasis)
    #print("eDualCoeffs: ", eDualCoeffs)
    fDualMatrix = dual(fBasis) 
    eDual = 0 # e*
    fDual = 0 # f*
    for k in range(len(eBasis)):
        eDual += eDualMatrix[0,k] * indicesToF(eBasis[k])
        fDual += fDualMatrix[0,k] * indicesToE(fBasis[k])
    return (eDual, fDual)


def compareRoots(μ,λ):
    if (μ > 0 and λ > 0):
        return 1 if (μ < λ) else 0
    elif (μ < 0 and λ < 0):
        return 2 if (μ < λ) else 0
    elif (μ == n and λ == -n ):
        return 0
    elif (μ > 0 and λ < 0):
        return 3
    else:
        return 0    



indexSet = [(μ, λ, compareRoots(μ,λ)) for μ in range(-n,n+1) for λ in range(-n,n+1) if compareRoots(μ,λ) in [1,2]]


def leftsum():
    print("Left sum:")
    sum = 0
    for i in range(1,n+1): # L_i
        secondterm = H[n-2] - H[n-1] # H_{n-1} - H_n
        for j in range(i,n):
            secondterm -= 2 * H[j-1]
        #summand = q**(2*n - 2*i) * (q**secondterm)
        summand = q**(-2*n + 2*i) * (q**secondterm)
        sum += summand
        print("+i: ", i, ", summand: ", summand)
        
    for i in range(1,n+1): # -L_i
        secondterm = - H[n-2] + H[n-1] # H_{n-1} - H_n
        for j in range(i,n):
            secondterm += 2 * H[j-1]
        #summand = q**(2*i - 2*n) * (q**secondterm)
        summand = q**(2*n - 2*i) * (q**secondterm)
        sum += summand
        print("-i: ", i, ", summand: ", summand)
    return sum

def rightsum():
    print("Right sum:")
    sum = 0
  
    for i in range(1,n+1):   # Computing dual elements e* and f* CASE 1
        for j in range(i+1,n+1): # μ = L_i > λ = L_j (i < j)
            pathSet = [x for x in range(i, j)] # e_{μλ} = E_{pathSet} = E_{i, ..., j - 1}
            #coeff = q**(1 + 2*n - 2*i) * (q - q**(-1))**(2 * len(pathSet))
            coeff = q**(1 - 2*n + 2*i) * (q - q**(-1))**(2 * len(pathSet))
            eDual, fDual = dualElements(pathSet)
            qExponent = H[n-2] - H[n-1] # H_{-L_i - L_j} = H_{-μ - λ}
            for k in range(i,n):
                qExponent -= 2 * H[k-1]
            for k in range(i,j):
                qExponent += H[k-1]
            summand = simplify((coeff * eDual * (q**qExponent) * fDual).subs(q**2 - 1, q*r).subs(q - 1/q, r))
            print("i: ", i, ", j: ", j, ", summand: ", summand)
            sum += summand
    

    for i in range(1,n+1): # CASE 2
        for j in range(1,n+1): 
            if (i == 1 and j == 1):
               continue
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
                eDual, fDual = dualElements(pathSet)
                qExponent = 0 # H_{-L_i + L_j} = H_{-μ - λ}
                for k in range(min(i,j), max(i,j)):
                    qExponent += H[k-1]
                if i < j:
                    qExponent *= -1
                summand = simplify(((-1) * coeff * eDual * (q**qExponent) * fDual).subs(q**2 - 1, q*r).subs(q - 1/q, r))
                print("i: ", i, ", -j: ", j, ", summand: ", summand)
                sum += summand

    for j in range(1,n+1):
        for i in range(j+1,n+1): # μ = -L_i > λ = -L_j (i > j)
            pathSet = [x for x in reversed(range(j, i))] # e_{μλ} = E_{pathSet} = E_{i-1, ..., j}
            #coeff = q**(1 + 2*i - 2*n) * (q - q**(-1))**(2 * len(pathSet))
            coeff = q**(1 - 2*i + 2*n) * (q - q**(-1))**(2 * len(pathSet))
            eDual, fDual = dualElements(pathSet)
            qExponent = - H[n-2] + H[n-1] # H_{L_i + L_j} = H_{-μ - λ}
            for k in range(i,n):
                qExponent += 2 * H[k-1]
            for k in range(j,i):
                qExponent += H[k-1]
            summand = simplify((coeff * eDual * (q**qExponent) * fDual).subs(q**2 - 1, q*r).subs(q - 1/q, r))
            print("-i: ", i, ", -j: ", j, ", summand: ", summand)
            sum += summand 

    return sum

print((leftsum() + rightsum()).subs(q-1/q,r))
#print(dual(result([1,2])))


'''   To go in rightsum
    for (μ,λ,c) in indexSet:
        match (c):
            case 1: 
                pathSet = [x for x in range(μ, λ)]
                coeff = q**(1 + 2*n - 2*μ + 2*(len(pathSet)))
            case 2: 
                coeff = q**(2 + 2*n - 2*μ) if μ == -1 * λ else q**(1 + 2*n - 2*μ)
                pathSet = []
                for x in range(-1 * λ,n-1):
                    pathSet.append(x)
                for x in reversed(range(μ,n+1)):
                    pathSet.append(x) # e_{μλ} = E_{pathSet} = E_{j, ..., n - 2, n, n - 1, ..., i}
            case 3: 
                coeff = q**(1 - 2*μ - 2*n)
                pathSet = [x for x in range(-1 * λ, -1 * μ)]
'''  
