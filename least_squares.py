from sympy import Symbol
from os import environ
import numpy as np
import numpy.linalg as la
import itertools
from collections import deque
import central_element_functions_numpy as ceNum

n = 3
qNums = [10**k for k in range(1,n + 2)] # for so2n we want powers of q to be {1, q, q^2, q^4, ..., q^2n}
A = np.zeros((n + 1, n + 1))
for i in range(n + 1):
    A[i][0] = 1
    A[i][1] = 10
    for j in range(2, n + 1):
        A[i][j] = qNums[i]**(2*j - 2)

APinv = la.pinv(A)

i = 1
j = 1
case = 2


pathSet = ceNum.getPathSet(n,i,j,case)

ceNum.q = 10
ceNum.n = n
eBasis = ceNum.result(pathSet)
B = np.zeros((len(qNums), len(eBasis)))


for q in range(len(qNums)):
    eDualCoeffs = ceNum.dual(eBasis)[0]
    for dual in range(len(eBasis)):
         ceNum.q = qNums[q]
         B[q][dual] = eDualCoeffs[dual]
print(f"B = \n{B}\n")

X = np.zeros((n + 1, np.shape(B)[1])) # should have a column for each column in B and enough rows to multiply with APinv
print(f"shape of X is {np.shape(X)} and shape of B is {np.shape(B)}\n")
for i in range(np.shape(B)[1]):
    print(B[:,i])
    
    #X[:,i] = APinv*B[:,i] # least squares solution to Ax_i = b_i is given by x_i = (A+)*b_i where A+ = APinv

q = Symbol('q')
Q = [1,q] + [q**k for k in [2,4,6]]


eDual = 0 # e*
for k in range(len(eBasis)):
    coeff = 0
    for (y,i) in zip(Q,range(n+1)):
        coeff += X[i][k]*y

    eDual +=  coeff * ceNum.indicesToF(eBasis[k])

print(eDual)
    



