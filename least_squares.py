from sympy import Symbol
from os import environ
import numpy as np
import numpy.linalg as la
import itertools
from collections import deque

n = 3
qNums = [10**k for k in range(1,n + 2)] # for so2n we want powers of q to be {1, q, q^2, q^4, ..., q^2n}
q = Symbol('q')

A = np.zeros((n + 1, n + 1))
for i in range(n + 1):
    A[i][0] = 1
    A[i][1] = 10
    for j in range(2, n + 1):
        A[i][j] = qNums[i]**(2*j - 2)

AInv = la.pinv(A)


