from sympy import *
from contextlib import redirect_stdout

q = Symbol('q')
a1 = Symbol('a1')
a2 = Symbol('a2')
n = Symbol('n', real=True)
#n = 5
f = open(f'baseCase_{n}.txt', 'w')


def duality(eta, xi):
    term = 1
    
    # site 1, class 1
    if ((eta[0] in [1,3]) and (xi[0] in [1,3])):
        x = 1
        if eta[1] in [1,3]:
            x += 1
        term *= (1 - q**(2*x) / a1)

    # site 1, class 2
    if ((eta[0] in [2, 3]) and (xi[0] in [2,3])):
        x = 1
        if eta[1] in [2,3]:
            x += 1
        term *= (1 - q**(2*x) / a2)

    # site 2, class 1
    if ((eta[1] in [1,3]) and (xi[1] in [1,3])):
        x = 2
        if xi[0] in [1,3]:
            x -= 1
        term *= (1 - q**(2*x) / a1)
        
    # site 2, class 2
    if ((eta[1] in [2, 3]) and (xi[1] in [2,3])):
        x = 2
        if xi[0] in [2,3]:
            x -= 1
        term *= (1 - q**(2*x) / a2)

    return term

def dualityMatrix():
    M = zeros(16)
    for i in range(16):
        for j in range(16):
            M[i, j] = duality((i // 4, i % 4), (j // 4, j % 4))
    return M

def generatorMatrix():
    L = zeros(16)
    
    L[1, 4] = q*(q**(1 - 2*n) + q**(2*n - 1))
    L[1, 1] = - L[1, 4]
    L[2, 8] = q*(q**(1 - 2*n) + q**(2*n - 1))
    L[2, 2] = - L[2, 8]
    L[3, 6] = 2*q**2 + q**(2 - 2*n) - q**(4 - 2*n)
    L[3, 9] = 2*q**2 + q**(2 - 2*n) - q**(4 - 2*n)
    L[3, 12] = q**2*(-q**(1 - n) + q**(n - 1))**2
    L[3, 3] = - L[3, 6] - L[3, 9] - L[3, 12]
    L[4, 1] = (q**(1 - 2*n) + q**(2*n - 1)) * q**(-1)
    L[4, 4] = - L[4, 1]
    L[6, 3] = q**(- 2*n) - q**(2 - 2*n) + 2
    L[6, 9] = (-q**(1 - n) + q**(n - 1))**2
    L[6, 12] = q**(2*n) - q**(2*n - 2) + 2
    L[6, 6] = - L[6, 3] - L[6, 9] - L[6, 12]
    L[7, 13] = q*(q**(1 - 2*n) + q**(2*n - 1))
    L[7, 7] = - L[7, 13]
    L[8, 2] = (q**(1 - 2*n) + q**(2*n - 1)) * q**(-1)
    L[8, 8] = - L[8, 2]
    L[9, 3] = q**(-2*n) - q**(2 - 2*n) + 2
    L[9, 6] = (-q**(1 - n) + q**(n - 1))**2
    L[9, 12] = q**(2*n) - q**(2*n - 2) + 2
    L[9, 9] = - L[9, 3] - L[9, 6] - L[9, 12]
    L[11, 14] = q*(q**(1 - 2*n) + q**(2*n - 1))
    L[11, 11] = - L[11, 14]
    L[12, 3] = (-q**(1 - n) + q**(n - 1))**2*q**(-2)
    L[12, 6] = q**(-2 + 2*n) - q**(-4 + 2*n) + 2*q**(-2)
    L[12, 9] = q**(-2 + 2*n) - q**(-4 + 2*n) + 2*q**(-2)
    L[12, 12] = - L[12, 3] - L[12,6] - L[12,9]
    L[13, 7] = (q**(1 - 2*n) + q**(2*n - 1))*q**(-1)
    L[13, 13] = - L[13, 7]
    L[14, 11] = (q**(1 - 2*n) + q**(2*n - 1))*q**(-1)
    L[14, 14] = - L[14, 11]

    return L
    
D = dualityMatrix()
L = generatorMatrix()

result = L * D - D * L.T
for i in range(16):
    for j in range(16):
        result[i, j] = collect(powsimp(powdenest(expand(result[i,j]))), q)


print(result, file = f)

#for i in range(len(L.row(0))):
#    print(L.row(i))
#print(L)

