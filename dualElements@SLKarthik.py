from sympy import *
from central_element_functions import *
import itertools
from collections import deque
f = open('eDualSymbolic.txt', '+')



# Takes in e and returns a list of the coefficients in e* (numerical)
# Case I: e_{L_i, L_j} (i < j)
def eDualNumericalI(n, i, j, qNum):
    pathSet = [x for x in range(i, j)] # e_{μλ} = E_{pathSet} = E_{i, ..., j - 1}
    
    eBasis = result(pathSet)
    eDualCoeffs = dual(eBasis).row(0).tolist()[0]
    return eDualCoeffs

# Takes in e and returns a list of the coefficients in e* (numerical)
# Case II: e_{L_i, -L_j}
def eDualNumericalII(n, i, j, qNum):
    pathSet = [] # e_{μλ} = +- E_{pathSet} = +- E_{i, ..., n-2, n, n-1, ..., j} (if i < n) or +- E_{n, n-2, ..., j} (if i = n) or +- E_{i, ..., n-2, n} (if j = n)
    for x in range(i,n-1):
        pathSet.append(x)
    pathSet.append(n)
    if (i != n and j != n):
        pathSet.append(n-1)
    for x in reversed(range(j,n-1)):
        pathSet.append(x) 
    
    eBasis = result(pathSet)
    eDualCoeffs = dual(eBasis).row(0).tolist()[0]
    return eDualCoeffs

# Takes in e and returns a list of the coefficients in e* (numerical)
# Case III: e_{-L_i, -L_j} (j > i)
def eDualNumericalIII(n, i, j, qNum):
    pathSet = [x for x in reversed(range(j, i))] # e_{μλ} = E_{pathSet} = E_{i-1, ..., j}
    eBasis = result(pathSet) # maybe we should store these
    eDualCoeffs = dual(eBasis).row(0).tolist()[0] # maybe we should store these
    return eDualCoeffs
    

# Takes in e and returns a list of the coefficients in e* (symbolic)
# Case I: e_{L_i, L_j} (i < j)
def eDualSymbolicI(n, i, j):
    #if (): # if the line exists in the file
        # return the entry from the file
    #else: 
        # interpolate to find the coefficients
        # write the coefficients to the file
        # return the coefficients