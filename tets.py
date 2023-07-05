from sympy import *
from contextlib import redirect_stdout

q = Symbol('q')
a1 = Symbol('a1')
a2 = Symbol('a2')
n = Symbol('n')
f = open('data.txt', 'w')


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
        if xi[1] in [2,3]:
            x -= 1
        term *= (1 - q**(2*x) / a2)

    return term
def dualityMatrix():
    M = zeros(16)
    for i in range(16):
        for j in range(16):
            M[i, j] = duality((i // 4, i % 4), (j // 4, j % 4))
    return M

D = dualityMatrix()

for i in range(16):
    print('|\t', file = f)
    for j in range(16):
        print(f"{D[i,j]}\t", file = f)
    print('\t|\n', file = f)