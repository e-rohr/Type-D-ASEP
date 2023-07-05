from sympy import Symbol
from os import environ
import numpy as np
import numpy.linalg as la
import itertools
from collections import deque
import central_element_functions_numpy as ceNum

n = 3
qNums = [10**k for k in range(1, 2*n + 1)]#[10**6, 10**5, 10**4, 10**3, 10**2, 10] # for so2n we want powers of q to be {1, q, q^2, q^4, ..., q^2n}
A = np.zeros((len(qNums), len(qNums)))
for i in range(len(qNums)):
    for j in range(len(qNums)):
        A[i][j] = qNums[i]**j


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

X = np.zeros((len(qNums), np.shape(B)[1])) # should have a column for each column in B and enough rows to multiply with APinv
for i in range(np.shape(B)[1]):
    X = APinv.dot(B) # least squares solution to Ax_i = b_i is given by x_i = (A+)*b_i where A+ = APinv  


q = Symbol('q')
Q = [q**k for k in range(2*n + 1)]
eDual = 0 # e*
for k in range(len(eBasis)):
    coeff = 0
    for (y,i) in zip(Q,range(len(qNums))):
        coeff += X[i][k]*y

    eDual +=  coeff * ceNum.indicesToF(eBasis[k])
    print(coeff * ceNum.indicesToF(eBasis[k]))
    



