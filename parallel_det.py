import sympy as sp
from multiprocessing import Pool
import itertools
from collections import deque


pool = multiprocess

def exec(matrix, level):
    if (level == 1):
        res = matrix[0,0]
    elif (level == 2):
        res  = matrix[0,0]*matrix[1,1] - matrix[0,1]*matrix[1,0]
    else:
        res = 0
    tasks = []
    for j1 in range(level):
        m = sp.zeros((level - 1, level - 1))
        for i in range(1,level):
            j2 = 0
            for j in range(level):
                if (j == j1):
                    continue
                m[i-1, j2] = matrix[i,j]
                j2 += 1

        



                         




