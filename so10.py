from sympy import *
from sympy.physics.quantum import TensorProduct
import time

start = time.time()
fout = open(f'so10_coproduct.txt', 'w')

######## Variables ########
n = 5
q = Symbol('q')
r = Symbol('r')

####### Matrices ########
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

# 1-indexed, returns Δ(E_index) as a 2n x 2n matrix
def coprodE(index):
    return TensorProduct(e(index), eye(2*n)) + TensorProduct(k(index), e(index))

# 1-indexed, returns Δ(F_index) as a 2n x 2n matrix
def coprodF(index):
    return TensorProduct(eye(2*n), f(index)) + TensorProduct(f(index), kNeg(index))

# 1-indexed, returns Δ(K_index) as a 2n x 2n matrix
def coprodK(index):
    return TensorProduct(k(index), k(index))

# 1-indexed, returns Δ(1/K_index) as a 2n x 2n matrix
def coprodKNeg(index):
    return TensorProduct(kNeg(index), kNeg(index))



########## Use for fundamental representation ###############
'''E1 = e(1)
E2 = e(2)
E3 = e(3)
E4 = e(4)
E5 = e(5)
F1 = f(1)
F2 = f(2)
F3 = f(3)
F4 = f(4)
F5 = f(5)
K1 = k(1)
K2 = k(2)
K3 = k(3)
K4 = k(4)
K5 = k(5)
KNeg1 = kNeg(1)
KNeg2 = kNeg(2)
KNeg3 = kNeg(3)
KNeg4 = kNeg(4)
KNeg5 = kNeg(5)
H1 = h(1)
H2 = h(2)
H3 = h(3)
H4 = h(4)
H5 = h(5)'''

########### Use for symbolic answer ###############
'''E1 = Symbol("E1", commutative = False)
E2 = Symbol("E2", commutative = False)
E3 = Symbol("E3", commutative = False)
E4 = Symbol("E4", commutative = False)
E5 = Symbol("E5", commutative = False)
F1 = Symbol("F1", commutative = False)
F2 = Symbol("F2", commutative = False)
F3 = Symbol("F3", commutative = False)
F4 = Symbol("F4", commutative = False)
F5 = Symbol("F5", commutative = False)
K1 = Symbol("K1", commutative = False)
K2 = Symbol("K2", commutative = False)
K3 = Symbol("K3", commutative = False)
K4 = Symbol("K4", commutative = False)
K5 = Symbol("K5", commutative = False)
KNeg1 = Symbol("KNeg1", commutative = False)
KNeg2 = Symbol("KNeg2", commutative = False)
KNeg3 = Symbol("KNeg3", commutative = False)
KNeg4 = Symbol("KNeg4", commutative = False)
KNeg5 = Symbol("KNeg5", commutative = False)'''

########### Use for coproduct term ##############
E1 = coprodE(1)
E2 = coprodE(2)
E3 = coprodE(3)
E4 = coprodE(4)
E5 = coprodE(5)
F1 = coprodF(1)
F2 = coprodF(2)
F3 = coprodF(3)
F4 = coprodF(4)
F5 = coprodF(5)
K1 = coprodK(1)
K2 = coprodK(2)
K3 = coprodK(3)
K4 = coprodK(4)
K5 = coprodK(5)
KNeg1 = coprodKNeg(1)
KNeg2 = coprodKNeg(2)
KNeg3 = coprodKNeg(3)
KNeg4 = coprodKNeg(4)
KNeg5 = coprodKNeg(5)



########### Left Sum ##########
left1 = q**(-8)*K4*KNeg5*KNeg1**2*KNeg2**2*KNeg3**2*KNeg4**2
left2 = q**(-6)*K4*KNeg5*KNeg2**2*KNeg3**2*KNeg4**2
left3 = q**(-4)*K4*KNeg5*KNeg3**2*KNeg4**2
left4 = q**(-2)*K4*KNeg5*KNeg4**2
left5 = K4*KNeg5
leftNeg1 = q**8*KNeg4*K5*K1**2*K2**2*K3**2*K4**2
leftNeg2 = q**6*KNeg4*K5*K2**2*K3**2*K4**2
leftNeg3 = q**4*KNeg4*K5*K3**2*K4**2
leftNeg4 = q**2*KNeg4*K5*K4**2
leftNeg5 = KNeg4*K5

########### Right Sum #########
# n, i, j, case
right5121 = ((q - 1/q)**2/q**7)*F1*K4*KNeg5*KNeg1**2*KNeg2**2*KNeg3**2*KNeg4**2*K1*E1
right5131 = ((q - 1/q)**4/q**7)*((-q/(q**2 - 1))*F2*F1 + (q**2/(q**2 - 1))*F1*F2)*K4*KNeg5*KNeg1**2*KNeg2**2*KNeg3**2*KNeg4**2*K1*K2*((-q/(q**2 - 1))*E1*E2 + (q**2/(q**2 - 1))*E2*E1)
right5141 = ((q - 1/q)**6/q**7)*((q**2/(q**4 - 2*q**2 + 1))*F3*F2*F1 + (-q**3/(q**4 - 2*q**2 + 1))*F1*F3*F2 + (-q**3/(q**4 - 2*q**2 + 1))*F2*F1*F3 + (q**4/(q**4 - 2*q**2 + 1))*F1*F2*F3)*K4*KNeg5*KNeg1**2*KNeg2**2*KNeg3**2*KNeg4**2*K1*K2*K3*((q**2/(q**4 - 2*q**2 + 1))*E1*E2*E3 + (-q**3/(q**4 - 2*q**2 + 1))*E2*E3*E1 + (-q**3/(q**4 - 2*q**2 + 1))*E3*E1*E2 + (q**4/(q**4 - 2*q**2 + 1))*E3*E2*E1)
right5151 = ((q - 1/q)**8/q**7)*((-q**3/(q**6 - 3*q**4 + 3*q**2 - 1))*F4*F3*F2*F1 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*F1*F4*F3*F2 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*F2*F1*F4*F3 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*F3*F2*F1*F4 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*F1*F2*F4*F3 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*F1*F3*F2*F4 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*F2*F1*F3*F4 + (q**6/(q**6 - 3*q**4 + 3*q**2 - 1))*F1*F2*F3*F4)*K4*KNeg5*KNeg1**2*KNeg2**2*KNeg3**2*KNeg4**2*K1*K2*K3*K4*((-q**3/(q**6 - 3*q**4 + 3*q**2 - 1))*E1*E2*E3*E4 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*E2*E3*E4*E1 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*E3*E4*E1*E2 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*E4*E1*E2*E3 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*E3*E4*E2*E1 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*E4*E2*E3*E1 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*E4*E3*E1*E2 + (q**6/(q**6 - 3*q**4 + 3*q**2 - 1))*E4*E3*E2*E1)
right5231 = ((q - 1/q)**2/q**5)*F2*K4*KNeg5*KNeg2**2*KNeg3**2*KNeg4**2*K2*E2
right5241 = ((q - 1/q)**4/q**5)*((-q/(q**2 - 1))*F3*F2 + (q**2/(q**2 - 1))*F2*F3)*K4*KNeg5*KNeg2**2*KNeg3**2*KNeg4**2*K2*K3*((-q/(q**2 - 1))*E2*E3 + (q**2/(q**2 - 1))*E3*E2)
right5251 = ((q - 1/q)**6/q**5)*((q**2/(q**4 - 2*q**2 + 1))*F4*F3*F2 + (-q**3/(q**4 - 2*q**2 + 1))*F2*F4*F3 + (-q**3/(q**4 - 2*q**2 + 1))*F3*F2*F4 + (q**4/(q**4 - 2*q**2 + 1))*F2*F3*F4)*K4*KNeg5*KNeg2**2*KNeg3**2*KNeg4**2*K2*K3*K4*((q**2/(q**4 - 2*q**2 + 1))*E2*E3*E4 + (-q**3/(q**4 - 2*q**2 + 1))*E3*E4*E2 + (-q**3/(q**4 - 2*q**2 + 1))*E4*E2*E3 + (q**4/(q**4 - 2*q**2 + 1))*E4*E3*E2)
right5341 = ((q - 1/q)**2/q**3)*F3*K4*KNeg5*KNeg3**2*KNeg4**2*K3*E3
right5351 = ((q - 1/q)**4/q**3)*((-q/(q**2 - 1))*F4*F3 + (q**2/(q**2 - 1))*F3*F4)*K4*KNeg5*KNeg3**2*KNeg4**2*K3*K4*((-q/(q**2 - 1))*E3*E4 + (q**2/(q**2 - 1))*E4*E3)
right5451 = ((q - 1/q)**2/q)*F4*K4*KNeg5*KNeg4**2*K4*E4
right5142 = (-(q - 1/q)**10/q**7)*((q**4/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F5*F4*F3*F2*F1 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F1*F5*F4*F3*F2 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F2*F1*F5*F4*F3 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F4*F3*F2*F1*F5 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F5*F3*F2*F1*F4 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F1*F2*F5*F4*F3 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F1*F4*F3*F2*F5 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F1*F5*F3*F2*F4 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F2*F1*F4*F3*F5 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F2*F1*F5*F3*F4 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F3*F2*F1*F5*F4 + (-q**7/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F1*F2*F4*F3*F5 + (-q**7/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F1*F2*F5*F3*F4 + (-q**7/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F1*F3*F2*F5*F4 + (-q**7/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F2*F1*F3*F5*F4 + (q**8/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F1*F2*F3*F5*F4)*KNeg1*KNeg2*KNeg3*((q**4/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E1*E2*E3*E4*E5 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E2*E3*E4*E5*E1 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E3*E4*E5*E1*E2 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E4*E1*E2*E3*E5 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E5*E1*E2*E3*E4 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E3*E4*E5*E2*E1 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E4*E2*E3*E5*E1 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E4*E3*E5*E1*E2 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E4*E5*E1*E2*E3 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E5*E2*E3*E4*E1 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E5*E3*E4*E1*E2 + (-q**7/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E4*E3*E5*E2*E1 + (-q**7/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E4*E5*E2*E3*E1 + (-q**7/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E4*E5*E3*E1*E2 + (-q**7/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E5*E3*E4*E2*E1 + (q**8/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E4*E5*E3*E2*E1)
right5152 = (-(q - 1/q)**8/q**7)*((-q**3/(q**6 - 3*q**4 + 3*q**2 - 1))*F5*F3*F2*F1 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*F1*F5*F3*F2 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*F2*F1*F5*F3 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*F3*F2*F1*F5 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*F1*F2*F5*F3 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*F1*F3*F2*F5 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*F2*F1*F3*F5 + (q**6/(q**6 - 3*q**4 + 3*q**2 - 1))*F1*F2*F3*F5)*KNeg1*KNeg2*KNeg3*KNeg4*((-q**3/(q**6 - 3*q**4 + 3*q**2 - 1))*E1*E2*E3*E5 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*E2*E3*E5*E1 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*E3*E5*E1*E2 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*E5*E1*E2*E3 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*E3*E5*E2*E1 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*E5*E2*E3*E1 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*E5*E3*E1*E2 + (q**6/(q**6 - 3*q**4 + 3*q**2 - 1))*E5*E3*E2*E1)
right5222 = (-(q - 1/q)**12/q**4)*((-q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F2**2*F3*F5*F4*F3 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F2*F3*F5*F2*F3*F4 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F2*F4*F3*F5*F2*F3 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F3*F2*F4*F3*F5*F2 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F5*F3*F2*F4*F3*F2 + ((-q**7 - q**5)/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F2*F3*F4*F3*F5*F2 + ((-q**7 - q**5)/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F2*F5*F3*F4*F3*F2 + ((-q**7 - q**5)/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F3*F2*F5*F4*F3*F2 + ((q**8 + q**6 + q**4)/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F2*F3*F5*F4*F3*F2)*((-q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E2**2*E3*E5*E4*E3 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E2*E3*E5*E2*E3*E4 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E2*E4*E3*E5*E2*E3 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E3*E2*E4*E3*E5*E2 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E5*E3*E2*E4*E3*E2 + ((-q**7 - q**5)/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E2*E3*E4*E3*E5*E2 + ((-q**7 - q**5)/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E2*E5*E3*E4*E3*E2 + ((-q**7 - q**5)/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E3*E2*E5*E4*E3*E2 + ((q**8 + q**6 + q**4)/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E2*E3*E5*E4*E3*E2)
right5232 = (-(q - 1/q)**10/q**5)*((q**4/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F3*F5*F4*F3*F2 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F3*F4*F3*F2*F5 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F5*F3*F4*F3*F2 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F3*F2*F4*F3*F5 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F3*F2*F5*F3*F4 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F4*F3*F2*F5*F3 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F5*F3*F2*F4*F3 + (-q**7/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F2*F3*F4*F3*F5 + (-q**7/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F2*F5*F3*F4*F3 + (q**8/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F2*F3*F5*F4*F3 + ((-q**7 - q**5)/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F3*F2*F5*F4*F3)*KNeg2*((q**4/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E2*E3*E4*E5*E3 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*E3*E5*E3*E4*E2 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*E4*E3*E5*E3*E2 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E3*E2*E3*E4*E5 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E4*E3*E2*E3*E5 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E4*E5*E3*E2*E3 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E5*E3*E2*E3*E4 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E3*E4*E2*E3*E5 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E3*E5*E2*E3*E4 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E4*E3*E5*E2*E3 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E5*E3*E4*E2*E3 + (q**8/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E3*E4*E5*E3*E2 + ((-q**7 - q**5)/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E3*E4*E5*E2*E3)
right5242 = (-(q - 1/q)**8/q**5)*((-q**3/(q**6 - 3*q**4 + 3*q**2 - 1))*F5*F4*F3*F2 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*F2*F5*F4*F3 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*F4*F3*F2*F5 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*F5*F3*F2*F4 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*F2*F4*F3*F5 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*F2*F5*F3*F4 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*F3*F2*F5*F4 + (q**6/(q**6 - 3*q**4 + 3*q**2 - 1))*F2*F3*F5*F4)*KNeg2*KNeg3*((-q**3/(q**6 - 3*q**4 + 3*q**2 - 1))*E2*E3*E4*E5 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*E3*E4*E5*E2 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*E4*E2*E3*E5 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*E5*E2*E3*E4 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*E4*E3*E5*E2 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*E4*E5*E2*E3 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*E5*E3*E4*E2 + (q**6/(q**6 - 3*q**4 + 3*q**2 - 1))*E4*E5*E3*E2)
right5252 = (-(q - 1/q)**6/q**5)*((q**2/(q**4 - 2*q**2 + 1))*F5*F3*F2 + (-q**3/(q**4 - 2*q**2 + 1))*F2*F5*F3 + (-q**3/(q**4 - 2*q**2 + 1))*F3*F2*F5 + (q**4/(q**4 - 2*q**2 + 1))*F2*F3*F5)*KNeg2*KNeg3*KNeg4*((q**2/(q**4 - 2*q**2 + 1))*E2*E3*E5 + (-q**3/(q**4 - 2*q**2 + 1))*E3*E5*E2 + (-q**3/(q**4 - 2*q**2 + 1))*E5*E2*E3 + (q**4/(q**4 - 2*q**2 + 1))*E5*E3*E2)
right5322 = (-(q - 1/q)**10/q**3)*((q**4/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F2*F3*F4*F5*F3 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*F3*F5*F3*F4*F2 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*F4*F3*F5*F3*F2 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F3*F2*F3*F4*F5 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F4*F3*F2*F3*F5 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F4*F5*F3*F2*F3 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F5*F3*F2*F3*F4 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F3*F4*F2*F3*F5 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F3*F5*F2*F3*F4 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F4*F3*F5*F2*F3 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F5*F3*F4*F2*F3 + (q**8/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F3*F4*F5*F3*F2 + ((-q**7 - q**5)/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F3*F4*F5*F2*F3)*K2*((q**4/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E3*E5*E4*E3*E2 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E3*E4*E3*E2*E5 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E5*E3*E4*E3*E2 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E3*E2*E4*E3*E5 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E3*E2*E5*E3*E4 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E4*E3*E2*E5*E3 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E5*E3*E2*E4*E3 + (-q**7/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E2*E3*E4*E3*E5 + (-q**7/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E2*E5*E3*E4*E3 + (q**8/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E2*E3*E5*E4*E3 + ((-q**7 - q**5)/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E3*E2*E5*E4*E3)
right5332 = (-(q - 1/q)**8/q**2)*((-q**3/(q**4 - 2*q**2 + 1))*F3*F4*F3*F5 + (-q**3/(q**4 - 2*q**2 + 1))*F5*F3*F4*F3 + ((q**4 + q**2)/(q**4 - 2*q**2 + 1))*F3*F5*F4*F3)*((-q**3/(q**4 - 2*q**2 + 1))*E3*E4*E3*E5 + (-q**3/(q**4 - 2*q**2 + 1))*E5*E3*E4*E3 + ((q**4 + q**2)/(q**4 - 2*q**2 + 1))*E3*E5*E4*E3)
right5342 = (-(q - 1/q)**6/q**3)*((q**2/(q**4 - 2*q**2 + 1))*F5*F4*F3 + (-q**3/(q**4 - 2*q**2 + 1))*F4*F3*F5 + (-q**3/(q**4 - 2*q**2 + 1))*F5*F3*F4 + (q**4/(q**4 - 2*q**2 + 1))*F3*F5*F4)*KNeg3*((q**2/(q**4 - 2*q**2 + 1))*E3*E4*E5 + (-q**3/(q**4 - 2*q**2 + 1))*E4*E3*E5 + (-q**3/(q**4 - 2*q**2 + 1))*E5*E3*E4 + (q**4/(q**4 - 2*q**2 + 1))*E4*E5*E3)
right5352 = (-(q - 1/q)**4/q**3)*((-q/(q**2 - 1))*F5*F3 + (q**2/(q**2 - 1))*F3*F5)*KNeg3*KNeg4*((-q/(q**2 - 1))*E3*E5 + (q**2/(q**2 - 1))*E5*E3)
right5412 = (-(q - 1/q)**10/q)*((q**4/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F1*F2*F3*F4*F5 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F2*F3*F4*F5*F1 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F3*F4*F5*F1*F2 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F4*F1*F2*F3*F5 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F5*F1*F2*F3*F4 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F3*F4*F5*F2*F1 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F4*F2*F3*F5*F1 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F4*F3*F5*F1*F2 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F4*F5*F1*F2*F3 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F5*F2*F3*F4*F1 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F5*F3*F4*F1*F2 + (-q**7/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F4*F3*F5*F2*F1 + (-q**7/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F4*F5*F2*F3*F1 + (-q**7/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F4*F5*F3*F1*F2 + (-q**7/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F5*F3*F4*F2*F1 + (q**8/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F4*F5*F3*F2*F1)*K1*K2*K3*((q**4/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E5*E4*E3*E2*E1 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E1*E5*E4*E3*E2 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E2*E1*E5*E4*E3 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E4*E3*E2*E1*E5 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E5*E3*E2*E1*E4 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E1*E2*E5*E4*E3 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E1*E4*E3*E2*E5 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E1*E5*E3*E2*E4 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E2*E1*E4*E3*E5 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E2*E1*E5*E3*E4 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E3*E2*E1*E5*E4 + (-q**7/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E1*E2*E4*E3*E5 + (-q**7/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E1*E2*E5*E3*E4 + (-q**7/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E1*E3*E2*E5*E4 + (-q**7/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E2*E1*E3*E5*E4 + (q**8/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E1*E2*E3*E5*E4)
right5422 = (-(q - 1/q)**8/q)*((-q**3/(q**6 - 3*q**4 + 3*q**2 - 1))*F2*F3*F4*F5 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*F3*F4*F5*F2 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*F4*F2*F3*F5 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*F5*F2*F3*F4 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*F4*F3*F5*F2 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*F4*F5*F2*F3 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*F5*F3*F4*F2 + (q**6/(q**6 - 3*q**4 + 3*q**2 - 1))*F4*F5*F3*F2)*K2*K3*((-q**3/(q**6 - 3*q**4 + 3*q**2 - 1))*E5*E4*E3*E2 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*E2*E5*E4*E3 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*E4*E3*E2*E5 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*E5*E3*E2*E4 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*E2*E4*E3*E5 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*E2*E5*E3*E4 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*E3*E2*E5*E4 + (q**6/(q**6 - 3*q**4 + 3*q**2 - 1))*E2*E3*E5*E4)
right5432 = (-(q - 1/q)**6/q)*((q**2/(q**4 - 2*q**2 + 1))*F3*F4*F5 + (-q**3/(q**4 - 2*q**2 + 1))*F4*F3*F5 + (-q**3/(q**4 - 2*q**2 + 1))*F5*F3*F4 + (q**4/(q**4 - 2*q**2 + 1))*F4*F5*F3)*K3*((q**2/(q**4 - 2*q**2 + 1))*E5*E4*E3 + (-q**3/(q**4 - 2*q**2 + 1))*E4*E3*E5 + (-q**3/(q**4 - 2*q**2 + 1))*E5*E3*E4 + (q**4/(q**4 - 2*q**2 + 1))*E3*E5*E4)
right5442 = (-(q - 1/q)**4)*F5*F4*E5*E4
right5452 = (-(q - 1/q)**2/q)*F5*KNeg4*E5
right5512 = (-q*(q - 1/q)**8)*((-q**3/(q**6 - 3*q**4 + 3*q**2 - 1))*F1*F2*F3*F5 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*F2*F3*F5*F1 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*F3*F5*F1*F2 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*F5*F1*F2*F3 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*F3*F5*F2*F1 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*F5*F2*F3*F1 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*F5*F3*F1*F2 + (q**6/(q**6 - 3*q**4 + 3*q**2 - 1))*F5*F3*F2*F1)*K1*K2*K3*K4*((-q**3/(q**6 - 3*q**4 + 3*q**2 - 1))*E5*E3*E2*E1 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*E1*E5*E3*E2 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*E2*E1*E5*E3 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*E3*E2*E1*E5 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*E1*E2*E5*E3 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*E1*E3*E2*E5 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*E2*E1*E3*E5 + (q**6/(q**6 - 3*q**4 + 3*q**2 - 1))*E1*E2*E3*E5)
right5522 = (-q*(q - 1/q)**6)*((q**2/(q**4 - 2*q**2 + 1))*F2*F3*F5 + (-q**3/(q**4 - 2*q**2 + 1))*F3*F5*F2 + (-q**3/(q**4 - 2*q**2 + 1))*F5*F2*F3 + (q**4/(q**4 - 2*q**2 + 1))*F5*F3*F2)*K2*K3*K4*((q**2/(q**4 - 2*q**2 + 1))*E5*E3*E2 + (-q**3/(q**4 - 2*q**2 + 1))*E2*E5*E3 + (-q**3/(q**4 - 2*q**2 + 1))*E3*E2*E5 + (q**4/(q**4 - 2*q**2 + 1))*E2*E3*E5)
right5532 = (-q*(q - 1/q)**4)*((-q/(q**2 - 1))*F3*F5 + (q**2/(q**2 - 1))*F5*F3)*K3*K4*((-q/(q**2 - 1))*E5*E3 + (q**2/(q**2 - 1))*E3*E5)
right5542 = (-q*(q - 1/q)**2)*F5*K4*E5
right5213 = (q**7*(q - 1/q)**2)*F1*KNeg4*K5*K2**2*K3**2*K4**2*K1*E1
right5313 = (q**5*(q - 1/q)**4)*((-q/(q**2 - 1))*F1*F2 + (q**2/(q**2 - 1))*F2*F1)*KNeg4*K5*K3**2*K4**2*K1*K2*((-q/(q**2 - 1))*E2*E1 + (q**2/(q**2 - 1))*E1*E2)
right5413 = (q**3*(q - 1/q)**6)*((q**2/(q**4 - 2*q**2 + 1))*F1*F2*F3 + (-q**3/(q**4 - 2*q**2 + 1))*F2*F3*F1 + (-q**3/(q**4 - 2*q**2 + 1))*F3*F1*F2 + (q**4/(q**4 - 2*q**2 + 1))*F3*F2*F1)*KNeg4*K5*K4**2*K1*K2*K3*((q**2/(q**4 - 2*q**2 + 1))*E3*E2*E1 + (-q**3/(q**4 - 2*q**2 + 1))*E1*E3*E2 + (-q**3/(q**4 - 2*q**2 + 1))*E2*E1*E3 + (q**4/(q**4 - 2*q**2 + 1))*E1*E2*E3)
right5513 = (q*(q - 1/q)**8)*((-q**3/(q**6 - 3*q**4 + 3*q**2 - 1))*F1*F2*F3*F4 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*F2*F3*F4*F1 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*F3*F4*F1*F2 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*F4*F1*F2*F3 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*F3*F4*F2*F1 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*F4*F2*F3*F1 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*F4*F3*F1*F2 + (q**6/(q**6 - 3*q**4 + 3*q**2 - 1))*F4*F3*F2*F1)*KNeg4*K5*K1*K2*K3*K4*((-q**3/(q**6 - 3*q**4 + 3*q**2 - 1))*E4*E3*E2*E1 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*E1*E4*E3*E2 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*E2*E1*E4*E3 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*E3*E2*E1*E4 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*E1*E2*E4*E3 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*E1*E3*E2*E4 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*E2*E1*E3*E4 + (q**6/(q**6 - 3*q**4 + 3*q**2 - 1))*E1*E2*E3*E4)
right5323 = (q**5*(q - 1/q)**2)*F2*KNeg4*K5*K3**2*K4**2*K2*E2
right5423 = (q**3*(q - 1/q)**4)*((-q/(q**2 - 1))*F2*F3 + (q**2/(q**2 - 1))*F3*F2)*KNeg4*K5*K4**2*K2*K3*((-q/(q**2 - 1))*E3*E2 + (q**2/(q**2 - 1))*E2*E3)
right5523 = (q*(q - 1/q)**6)*((q**2/(q**4 - 2*q**2 + 1))*F2*F3*F4 + (-q**3/(q**4 - 2*q**2 + 1))*F3*F4*F2 + (-q**3/(q**4 - 2*q**2 + 1))*F4*F2*F3 + (q**4/(q**4 - 2*q**2 + 1))*F4*F3*F2)*KNeg4*K5*K2*K3*K4*((q**2/(q**4 - 2*q**2 + 1))*E4*E3*E2 + (-q**3/(q**4 - 2*q**2 + 1))*E2*E4*E3 + (-q**3/(q**4 - 2*q**2 + 1))*E3*E2*E4 + (q**4/(q**4 - 2*q**2 + 1))*E2*E3*E4)
right5433 = (q**3*(q - 1/q)**2)*F3*KNeg4*K5*K4**2*K3*E3
right5533 = (q*(q - 1/q)**4)*((-q/(q**2 - 1))*F3*F4 + (q**2/(q**2 - 1))*F4*F3)*KNeg4*K5*K3*K4*((-q/(q**2 - 1))*E4*E3 + (q**2/(q**2 - 1))*E3*E4)
right5543 = (q*(q - 1/q)**2)*F4*KNeg4*K5*K4*E4
########## GUESSES ##########
right5112 = - q**(-6) * r**16 * r**(-6) * (-q**3 *F1*F2*F3*F4*F1*F2*F3*F5 - q * (q**2 + 1)**2 *F1*F2*F3*F4*F3*F5*F2*F1 + q**2 * (q**2 + 1) *F1*F2*F3*F5*F2*F3*F4*F1 + q**2 * (q**2 + 1) *F1*F2*F3*F5*F4*F1*F2*F3 + q**3 *F1*F2*F3*F5*F4*F3*F1*F2 + (q**6 + q**4 + q**2 + 1)*F1*F2*F3*F5*F4*F3*F2*F1 + q**2 * (q**2 + 1) *F1*F2*F4*F3*F5*F2*F3*F1 - q**3 *F1*F2*F5*F3*F4*F1*F2*F3 - q * (q**2 + 1)**2 *F1*F2*F5*F3*F4*F3*F2*F1 - q**2 * (q**2 + 1) *F1*F2**2*F3*F5*F4*F3*F1 + q**3 *F1*F3*F2*F1*F2*F3*F5*F4 + q**2 * (q**2 + 1) *F1*F3*F2*F4*F3*F5*F2*F1 - q**3 *F1*F3*F2*F5*F3*F4*F1*F2 - q * (q**2 + 1)**2 *F1*F3*F2*F5*F4*F3*F2*F1 + q**3 *F1*F3*F5*F4*F3*F2*F1*F2 + q**3 *F1*F4*F3*F2*F1*F2*F3*F5 - q**3 *F1*F4*F3*F2*F5*F3*F1*F2 + q**3 *F1*F5*F3*F2*F1*F2*F3*F4 + q**2 * (q**2 + 1) *F1*F5*F3*F2*F4*F3*F2*F1 + q**3 *F1*F5*F4*F3*F2*F1*F2*F3 - q**2 * (q**2 + 1) *F1**2*F2*F3*F5*F4*F3*F2 + q**3 *F2*F1*F2*F3*F5*F4*F3*F1 + q**2 * (q**2 + 1) *F2*F1*F3*F4*F3*F5*F2*F1 - q**3 *F2*F1*F3*F5*F2*F3*F4*F1 - q*(q**4 + q**2 + 1) *F2*F1*F3*F5*F4*F3*F2*F1 - q**3 *F2*F1*F4*F3*F5*F2*F3*F1 + q**2 * (q**2 + 1) *F2*F1*F5*F3*F4*F3*F2*F1 - q**3 *F3*F2*F1*F4*F3*F5*F2*F1 + q**2 * (q**2 + 1) *F3*F2*F1*F5*F4*F3*F2*F1 - q**3 *F5*F3*F2*F1*F4*F3*F2*F1) * r**(-6) * (-q**3 *E1*E2*E3*E4*E1*E2*E3*E5 - q * (q**2 + 1)**2 *E1*E2*E3*E4*E3*E5*E2*E1 + q**2 * (q**2 + 1) *E1*E2*E3*E5*E2*E3*E4*E1 + q**2 * (q**2 + 1) *E1*E2*E3*E5*E4*E1*E2*E3 + q**3 *E1*E2*E3*E5*E4*E3*E1*E2 + (q**6 + q**4 + q**2 + 1)*E1*E2*E3*E5*E4*E3*E2*E1 + q**2 * (q**2 + 1) *E1*E2*E4*E3*E5*E2*E3*E1 - q**3 *E1*E2*E5*E3*E4*E1*E2*E3 - q * (q**2 + 1)**2 *E1*E2*E5*E3*E4*E3*E2*E1 - q**2 * (q**2 + 1) *E1*E2**2*E3*E5*E4*E3*E1 + q**3 *E1*E3*E2*E1*E2*E3*E5*E4 + q**2 * (q**2 + 1) *E1*E3*E2*E4*E3*E5*E2*E1 - q**3 *E1*E3*E2*E5*E3*E4*E1*E2 - q * (q**2 + 1)**2 *E1*E3*E2*E5*E4*E3*E2*E1 + q**3 *E1*E3*E5*E4*E3*E2*E1*E2 + q**3 *E1*E4*E3*E2*E1*E2*E3*E5 - q**3 *E1*E4*E3*E2*E5*E3*E1*E2 + q**3 *E1*E5*E3*E2*E1*E2*E3*E4 + q**2 * (q**2 + 1) *E1*E5*E3*E2*E4*E3*E2*E1 + q**3 *E1*E5*E4*E3*E2*E1*E2*E3 - q**2 * (q**2 + 1) *E1**2*E2*E3*E5*E4*E3*E2 + q**3 *E2*E1*E2*E3*E5*E4*E3*E1 + q**2 * (q**2 + 1) *E2*E1*E3*E4*E3*E5*E2*E1 - q**3 *E2*E1*E3*E5*E2*E3*E4*E1 - q*(q**4 + q**2 + 1) *E2*E1*E3*E5*E4*E3*E2*E1 - q**3 *E2*E1*E4*E3*E5*E2*E3*E1 + q**2 * (q**2 + 1) *E2*E1*E5*E3*E4*E3*E2*E1 - q**3 *E3*E2*E1*E4*E3*E5*E2*E1 + q**2 * (q**2 + 1) *E3*E2*E1*E5*E4*E3*E2*E1 - q**3 *E5*E3*E2*E1*E4*E3*E2*E1)
#right5112 = - q**(-6) * r**(16) * r**(-6) * (-q**3*F1*F2*F3*F4*F1*F2*F3*F5 - (q**5 + 2*q**3 + q)*F1*F2*F3*F4*F3*F5*F2*F1 + (q**4 + q**2)*F1*F2*F3*F5*F2*F3*F4*F1 + (q**4 + q**2)*F1*F2*F3*F5*F4*F1*F2*F3 + q**3*F1*F2*F3*F5*F4*F3*F1*F2 + (q**6 + q**4 + q**2 + 1)*F1*F2*F3*F5*F4*F3*F2*F1 + (q**4 + q**2)*F1*F2*F4*F3*F5*F2*F3*F1 - q**3*F1*F2*F5*F3*F4*F1*F2*F3 - (q**5 + 2*q**3 + q)*F1*F2*F5*F3*F4*F3*F2*F1 - (q**4 + q**2)*F1*F2**2*F3*F5*F4*F3*F1 + q**3*F1*F3*F2*F1*F2*F3*F5*F4 + (q**4 + q**2)*F1*F3*F2*F4*F3*F5*F2*F1 - q**3*F1*F3*F2*F5*F3*F4*F1*F2 - (q**5 + 2*q**3 + q)*F1*F3*F2*F5*F4*F3*F2*F1 + q**3*F1*F3*F5*F4*F3*F2*F1*F2 + q**3*F1*F4*F3*F2*F1*F2*F3*F5 - q**3*F1*F4*F3*F2*F5*F3*F1*F2 + q**3*F1*F5*F3*F2*F1*F2*F3*F4 + (q**4 + q**2)*F1*F5*F3*F2*F4*F3*F2*F1 + q**3*F1*F5*F4*F3*F2*F1*F2*F3 - (q**4 + q**2)*F1**2*F2*F3*F5*F4*F3*F2 + q**3*F2*F1*F2*F3*F5*F4*F3*F1 + (q**4 + q**2)*F2*F1*F3*F4*F3*F5*F2*F1 - q**3*F2*F1*F3*F5*F2*F3*F4*F1 - (q**5 + q**3 + q)*F2*F1*F3*F5*F4*F3*F2*F1 - q**3*F2*F1*F4*F3*F5*F2*F3*F1 + (q**4 + q**2)*F2*F1*F5*F3*F4*F3*F2*F1 - q**3*F3*F2*F1*F4*F3*F5*F2*F1 + (q**4 + q**2)*F3*F2*F1*F5*F4*F3*F2*F1 - q**3*F5*F3*F2*F1*F4*F3*F2*F1) * r**(-6) * (-q**3*E1*E2*E3*E4*E1*E2*E3*E5 - (q**5 + 2*q**3 + q)*E1*E2*E3*E4*E3*E5*E2*E1 + (q**4 + q**2)*E1*E2*E3*E5*E2*E3*E4*E1 + (q**4 + q**2)*E1*E2*E3*E5*E4*E1*E2*E3 + q**3*E1*E2*E3*E5*E4*E3*E1*E2 + (q**6 + q**4 + q**2 + 1)*E1*E2*E3*E5*E4*E3*E2*E1 + (q**4 + q**2)*E1*E2*E4*E3*E5*E2*E3*E1 - q**3*E1*E2*E5*E3*E4*E1*E2*E3 - (q**5 + 2*q**3 + q)*E1*E2*E5*E3*E4*E3*E2*E1 - (q**4 + q**2)*E1*E2**2*E3*E5*E4*E3*E1 + q**3*E1*E3*E2*E1*E2*E3*E5*E4 + (q**4 + q**2)*E1*E3*E2*E4*E3*E5*E2*E1 - q**3*E1*E3*E2*E5*E3*E4*E1*E2 - (q**5 + 2*q**3 + q)*E1*E3*E2*E5*E4*E3*E2*E1 + q**3*E1*E3*E5*E4*E3*E2*E1*E2 + q**3*E1*E4*E3*E2*E1*E2*E3*E5 - q**3*E1*E4*E3*E2*E5*E3*E1*E2 + q**3*E1*E5*E3*E2*E1*E2*E3*E4 + (q**4 + q**2)*E1*E5*E3*E2*E4*E3*E2*E1 + q**3*E1*E5*E4*E3*E2*E1*E2*E3 - (q**4 + q**2)*E1**2*E2*E3*E5*E4*E3*E2 + q**3*E2*E1*E2*E3*E5*E4*E3*E1 + (q**4 + q**2)*E2*E1*E3*E4*E3*E5*E2*E1 - q**3*E2*E1*E3*E5*E2*E3*E4*E1 - (q**5 + q**3 + q)*E2*E1*E3*E5*E4*E3*E2*E1 - q**3*E2*E1*E4*E3*E5*E2*E3*E1 + (q**4 + q**2)*E2*E1*E5*E3*E4*E3*E2*E1 - q**3*E3*E2*E1*E4*E3*E5*E2*E1 + (q**4 + q**2)*E3*E2*E1*E5*E4*E3*E2*E1 - q**3*E5*E3*E2*E1*E4*E3*E2*E1)
right5122 = - q**(-7) * r**(14) * r**(-6) * (-q**5*F1*F2*F3*F4*F3*F5*F2 + q**4*F1*F2*F3*F5*F2*F3*F4 - q**3*F1*F2*F3*F5*F4*F2*F3 + q**6*F1*F2*F3*F5*F4*F3*F2 + q**4*F1*F2*F4*F3*F5*F2*F3 - q**5*F1*F2*F5*F3*F4*F3*F2 - (q**4 - q**2)*F1*F2**2*F3*F5*F4*F3 + q**4*F1*F3*F2*F4*F3*F5*F2 - q**5*F1*F3*F2*F5*F4*F3*F2 + q**4*F1*F5*F3*F2*F4*F3*F2 + (q**3 - q)*F2*F1*F2*F3*F5*F4*F3 - q**3*F2*F1*F3*F4*F2*F3*F5 + (q**4 + q**2)*F2*F1*F3*F4*F3*F5*F2 - q**3*F2*F1*F3*F5*F2*F3*F4 + (q**4 + q**2)*F2*F1*F3*F5*F4*F2*F3 - (q**5 + q)*F2*F1*F3*F5*F4*F3*F2 - q**3*F2*F1*F4*F3*F5*F2*F3 - q**3*F2*F1*F5*F3*F4*F2*F3 + (q**4 + q**2)*F2*F1*F5*F3*F4*F3*F2 - q*F2*F3*F4*F3*F5*F2*F1 + q**2*F2*F3*F5*F2*F1*F3*F4 - q**3*F2*F3*F5*F4*F2*F1*F3 + 1*F2*F3*F5*F4*F3*F2*F1 + q**2*F2*F4*F3*F5*F2*F1*F3 - q*F2*F5*F3*F4*F3*F2*F1 - q**3*F3*F2*F1*F4*F3*F5*F2 - q**3*F3*F2*F1*F5*F3*F4*F2 + (q**4 + q**2)*F3*F2*F1*F5*F4*F3*F2 + q**2*F3*F2*F4*F3*F5*F2*F1 - q*F3*F2*F5*F4*F3*F2*F1 - q**3*F4*F3*F2*F1*F5*F3*F2 - q**3*F5*F3*F2*F1*F4*F3*F2 + q**2*F5*F3*F2*F4*F3*F2*F1) * KNeg1 * r**(-6) * (1*E1*E2*E3*E5*E4*E3*E2 - q*E2*E1*E2*E3*E5*E4*E3 + q**2*E2*E3*E1*E2*E3*E5*E4 - q**3*E2*E3*E4*E1*E2*E3*E5 + q**4*E2*E3*E4*E3*E5*E1*E2 - (q**5 - q)*E2*E3*E4*E3*E5*E2*E1 - q**3*E2*E3*E5*E1*E2*E3*E4 + (q**4 - q**2)*E2*E3*E5*E2*E3*E4*E1 + (q**4 + q**2)*E2*E3*E5*E4*E1*E2*E3 - (q**5 + q)*E2*E3*E5*E4*E3*E1*E2 + q**6*E2*E3*E5*E4*E3*E2*E1 + q**2*E2*E4*E3*E1*E2*E3*E5 - q**3*E2*E4*E3*E5*E1*E2*E3 + (q**4 - q**2)*E2*E4*E3*E5*E2*E3*E1 + q**2*E2*E5*E3*E1*E2*E3*E4 - q**3*E2*E5*E3*E4*E1*E2*E3 + q**4*E2*E5*E3*E4*E3*E1*E2 - (q**5 - q)*E2*E5*E3*E4*E3*E2*E1 + q**2*E2*E5*E4*E3*E1*E2*E3 - (q**4 - q**2)*E2**2*E3*E5*E4*E3*E1 - q*E3*E2*E1*E2*E3*E5*E4 + q**2*E3*E2*E4*E1*E2*E3*E5 - q**3*E3*E2*E4*E3*E5*E1*E2 + (q**4 - q**2)*E3*E2*E4*E3*E5*E2*E1 + q**2*E3*E2*E5*E1*E2*E3*E4 - q**3*E3*E2*E5*E3*E4*E1*E2 - (q**3 + q)*E3*E2*E5*E4*E1*E2*E3 + (q**4 + q**2)*E3*E2*E5*E4*E3*E1*E2 - (q**5 - q)*E3*E2*E5*E4*E3*E2*E1 + q**2*E3*E4*E3*E2*E5*E1*E2 - q*E3*E5*E4*E3*E2*E1*E2 - q*E4*E3*E2*E1*E2*E3*E5 + q**2*E4*E3*E2*E5*E1*E2*E3 - q**3*E4*E3*E2*E5*E3*E1*E2 - q*E5*E3*E2*E1*E2*E3*E4 + q**2*E5*E3*E2*E4*E1*E2*E3 - q**3*E5*E3*E2*E4*E3*E1*E2 + (q**4 - q**2)*E5*E3*E2*E4*E3*E2*E1 + q**2*E5*E3*E4*E3*E2*E1*E2 - q*E5*E4*E3*E2*E1*E2*E3)
right5212 = - q**(-5) * r**(14) * r**(-6) * (1*F1*F2*F3*F5*F4*F3*F2 - q*F2*F1*F2*F3*F5*F4*F3 + q**2*F2*F3*F1*F2*F3*F5*F4 - q**3*F2*F3*F4*F1*F2*F3*F5 + q**4*F2*F3*F4*F3*F5*F1*F2 - (q**5 - q)*F2*F3*F4*F3*F5*F2*F1 - q**3*F2*F3*F5*F1*F2*F3*F4 + (q**4 - q**2)*F2*F3*F5*F2*F3*F4*F1 + (q**4 + q**2)*F2*F3*F5*F4*F1*F2*F3 - (q**5 + q)*F2*F3*F5*F4*F3*F1*F2 + q**6*F2*F3*F5*F4*F3*F2*F1 + q**2*F2*F4*F3*F1*F2*F3*F5 - q**3*F2*F4*F3*F5*F1*F2*F3 + (q**4 - q**2)*F2*F4*F3*F5*F2*F3*F1 + q**2*F2*F5*F3*F1*F2*F3*F4 - q**3*F2*F5*F3*F4*F1*F2*F3 + q**4*F2*F5*F3*F4*F3*F1*F2 - (q**5 - q)*F2*F5*F3*F4*F3*F2*F1 + q**2*F2*F5*F4*F3*F1*F2*F3 - (q**4 - q**2)*F2**2*F3*F5*F4*F3*F1 - q*F3*F2*F1*F2*F3*F5*F4 + q**2*F3*F2*F4*F1*F2*F3*F5 - q**3*F3*F2*F4*F3*F5*F1*F2 + (q**4 - q**2)*F3*F2*F4*F3*F5*F2*F1 + q**2*F3*F2*F5*F1*F2*F3*F4 - q**3*F3*F2*F5*F3*F4*F1*F2 - (q**3 + q)*F3*F2*F5*F4*F1*F2*F3 + (q**4 + q**2)*F3*F2*F5*F4*F3*F1*F2 - (q**5 - q)*F3*F2*F5*F4*F3*F2*F1 + q**2*F3*F4*F3*F2*F5*F1*F2 - q*F3*F5*F4*F3*F2*F1*F2 - q*F4*F3*F2*F1*F2*F3*F5 + q**2*F4*F3*F2*F5*F1*F2*F3 - q**3*F4*F3*F2*F5*F3*F1*F2 - q*F5*F3*F2*F1*F2*F3*F4 + q**2*F5*F3*F2*F4*F1*F2*F3 - q**3*F5*F3*F2*F4*F3*F1*F2 + (q**4 - q**2)*F5*F3*F2*F4*F3*F2*F1 + q**2*F5*F3*F4*F3*F2*F1*F2 - q*F5*F4*F3*F2*F1*F2*F3) * K1 * r**(-6) * (-q**5*E1*E2*E3*E4*E3*E5*E2 + q**4*E1*E2*E3*E5*E2*E3*E4 - q**3*E1*E2*E3*E5*E4*E2*E3 + q**6*E1*E2*E3*E5*E4*E3*E2 + q**4*E1*E2*E4*E3*E5*E2*E3 - q**5*E1*E2*E5*E3*E4*E3*E2 - (q**4 - q**2)*E1*E2**2*E3*E5*E4*E3 + q**4*E1*E3*E2*E4*E3*E5*E2 - q**5*E1*E3*E2*E5*E4*E3*E2 + q**4*E1*E5*E3*E2*E4*E3*E2 + (q**3 - q)*E2*E1*E2*E3*E5*E4*E3 - q**3*E2*E1*E3*E4*E2*E3*E5 + (q**4 + q**2)*E2*E1*E3*E4*E3*E5*E2 - q**3*E2*E1*E3*E5*E2*E3*E4 + (q**4 + q**2)*E2*E1*E3*E5*E4*E2*E3 - (q**5 + q)*E2*E1*E3*E5*E4*E3*E2 - q**3*E2*E1*E4*E3*E5*E2*E3 - q**3*E2*E1*E5*E3*E4*E2*E3 + (q**4 + q**2)*E2*E1*E5*E3*E4*E3*E2 - q*E2*E3*E4*E3*E5*E2*E1 + q**2*E2*E3*E5*E2*E1*E3*E4 - q**3*E2*E3*E5*E4*E2*E1*E3 + 1*E2*E3*E5*E4*E3*E2*E1 + q**2*E2*E4*E3*E5*E2*E1*E3 - q*E2*E5*E3*E4*E3*E2*E1 - q**3*E3*E2*E1*E4*E3*E5*E2 - q**3*E3*E2*E1*E5*E3*E4*E2 + (q**4 + q**2)*E3*E2*E1*E5*E4*E3*E2 + q**2*E3*E2*E4*E3*E5*E2*E1 - q*E3*E2*E5*E4*E3*E2*E1 - q**3*E4*E3*E2*E1*E5*E3*E2 - q**3*E5*E3*E2*E1*E4*E3*E2 + q**2*E5*E3*E2*E4*E3*E2*E1)
right5132 = - q**(-7) * r**(12) * r**(-5) * (-q**4*F1*F2*F3*F4*F3*F5 + q**5*F1*F2*F3*F5*F4*F3 - q**4*F1*F2*F5*F3*F4*F3 + q**3*F1*F3*F2*F4*F3*F5 + q**3*F1*F3*F2*F5*F3*F4 - (q**4 + q**2)*F1*F3*F2*F5*F4*F3 - q**2*F1*F3*F4*F3*F2*F5 + q*F1*F3*F5*F4*F3*F2 + q**3*F1*F4*F3*F2*F5*F3 + q**3*F1*F5*F3*F2*F4*F3 - q**2*F1*F5*F3*F4*F3*F2 + q**3*F2*F1*F3*F4*F3*F5 - q**4*F2*F1*F3*F5*F4*F3 + q**3*F2*F1*F5*F3*F4*F3 - q**2*F3*F2*F1*F4*F3*F5 - q**2*F3*F2*F1*F5*F3*F4 + (q**3 + q)*F3*F2*F1*F5*F4*F3 + q*F3*F4*F3*F2*F1*F5 - 1*F3*F5*F4*F3*F2*F1 - q**2*F4*F3*F2*F1*F5*F3 - q**2*F5*F3*F2*F1*F4*F3 + q*F5*F3*F4*F3*F2*F1) * KNeg1 * KNeg2 * r**(-5) * (-1*E1*E2*E3*E5*E4*E3 + q*E2*E3*E5*E4*E3*E1 + q*E3*E1*E2*E3*E5*E4 - q**2*E3*E2*E3*E5*E4*E1 - q**2*E3*E4*E1*E2*E3*E5 + q**3*E3*E4*E2*E3*E5*E1 + (q**3 - q)*E3*E4*E3*E5*E1*E2 - (q**4 - q**2)*E3*E4*E3*E5*E2*E1 - q**2*E3*E5*E1*E2*E3*E4 + q**3*E3*E5*E2*E3*E4*E1 + (q**3 + q)*E3*E5*E4*E1*E2*E3 - (q**4 + q**2)*E3*E5*E4*E2*E3*E1 - q**4*E3*E5*E4*E3*E1*E2 + q**5*E3*E5*E4*E3*E2*E1 + q*E4*E3*E1*E2*E3*E5 - q**2*E4*E3*E2*E3*E5*E1 - q**2*E4*E3*E5*E1*E2*E3 + q**3*E4*E3*E5*E2*E3*E1 + q*E5*E3*E1*E2*E3*E4 - q**2*E5*E3*E2*E3*E4*E1 - q**2*E5*E3*E4*E1*E2*E3 + q**3*E5*E3*E4*E2*E3*E1 + (q**3 - q)*E5*E3*E4*E3*E1*E2 - (q**4 - q**2)*E5*E3*E4*E3*E2*E1 + q*E5*E4*E3*E1*E2*E3 - q**2*E5*E4*E3*E2*E3*E1)
right5312 = - q**(-3) * r**(12) * r**(-5) * (-1*F1*F2*F3*F5*F4*F3 + q*F2*F3*F5*F4*F3*F1 + q*F3*F1*F2*F3*F5*F4 - q**2*F3*F2*F3*F5*F4*F1 - q**2*F3*F4*F1*F2*F3*F5 + q**3*F3*F4*F2*F3*F5*F1 + (q**3 - q)*F3*F4*F3*F5*F1*F2 - (q**4 - q**2)*F3*F4*F3*F5*F2*F1 - q**2*F3*F5*F1*F2*F3*F4 + q**3*F3*F5*F2*F3*F4*F1 + (q**3 + q)*F3*F5*F4*F1*F2*F3 - (q**4 + q**2)*F3*F5*F4*F2*F3*F1 - q**4*F3*F5*F4*F3*F1*F2 + q**5*F3*F5*F4*F3*F2*F1 + q*F4*F3*F1*F2*F3*F5 - q**2*F4*F3*F2*F3*F5*F1 - q**2*F4*F3*F5*F1*F2*F3 + q**3*F4*F3*F5*F2*F3*F1 + q*F5*F3*F1*F2*F3*F4 - q**2*F5*F3*F2*F3*F4*F1 - q**2*F5*F3*F4*F1*F2*F3 + q**3*F5*F3*F4*F2*F3*F1 + (q**3 - q)*F5*F3*F4*F3*F1*F2 - (q**4 - q**2)*F5*F3*F4*F3*F2*F1 + q*F5*F4*F3*F1*F2*F3 - q**2*F5*F4*F3*F2*F3*F1) * K1 * K2 * r**(-5) * (-q**4*E1*E2*E3*E4*E3*E5 + q**5*E1*E2*E3*E5*E4*E3 - q**4*E1*E2*E5*E3*E4*E3 + q**3*E1*E3*E2*E4*E3*E5 + q**3*E1*E3*E2*E5*E3*E4 - (q**4 + q**2)*E1*E3*E2*E5*E4*E3 - q**2*E1*E3*E4*E3*E2*E5 + q*E1*E3*E5*E4*E3*E2 + q**3*E1*E4*E3*E2*E5*E3 + q**3*E1*E5*E3*E2*E4*E3 - q**2*E1*E5*E3*E4*E3*E2 + q**3*E2*E1*E3*E4*E3*E5 - q**4*E2*E1*E3*E5*E4*E3 + q**3*E2*E1*E5*E3*E4*E3 - q**2*E3*E2*E1*E4*E3*E5 - q**2*E3*E2*E1*E5*E3*E4 + (q**3 + q)*E3*E2*E1*E5*E4*E3 + q*E3*E4*E3*E2*E1*E5 - 1*E3*E5*E4*E3*E2*E1 - q**2*E4*E3*E2*E1*E5*E3 - q**2*E5*E3*E2*E1*E4*E3 + q*E5*E3*E4*E3*E2*E1)



########## PRINTING ############
#print(right5312.applyfunc(simplify))

leftSumPos = left1 + left2 + left3 + left4 + left5
leftSumNeg = leftNeg1 + leftNeg2 + leftNeg3 + leftNeg4 + leftNeg5
rightSum1 = right5121 + right5131 + right5141 + right5151 + right5231 + right5241 + right5251 + right5341 + right5351 + right5451
rightSum2 = right5112 + right5122 + right5132 + right5142 + right5152 + right5212 + right5222 + right5232 + right5242 + right5252 + right5312 + right5322 + right5332 + right5342 + right5352 + right5412 + right5422 + right5432 + right5442 + right5452 + right5512 + right5522 + right5532 + right5542
rightSum3 = right5213 + right5313 + right5413 + right5513 + right5323 + right5423 + right5523 + right5433 + right5533 + right5543

#representation = leftSumPos + leftSumNeg + rightSum1 + rightSum2 + rightSum3
#representation = rightSum2
#representation = representation.subs(q - 1/q, r).subs(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1, q**4*r**4).subs(q**6 - 3*q**4 + 3*q**2 - 1, q**3*r**3).subs(q**4 - 2*q**2 + 1, q**2*r**2).subs(q**2 - 1, q*r).subs(q**4 - q*r - 1, q**3*r).subs(q**4 + q*r + 1, q**4 + q**2)
#representation = powsimp(powdenest(representation))

#lat = open('so10latex.txt', 'w')
#print(latex(representation), file = lat)

coproductRepresentation = leftSumPos + leftSumNeg + rightSum1 + rightSum2 + rightSum3

indexList = [6, 17, 28, 39, 50, 51, 62, 73, 84, 95] # 10x10 block
for i in range(len(indexList)):
    coproductRepresentation = coproductRepresentation.elementary_row_op(op = "n<->m", row1 = i, row2 = indexList[i] - 1)
    coproductRepresentation = coproductRepresentation.elementary_col_op(op = "n<->m", col1 = i, col2 = indexList[i] - 1) 

fromIndex = [17, 22, 42, 72, 82, 28, 33, 43, 63, 83, 39, 36, 54, 44, 74, 94, 50, 53, 55, 65, 85, 85, 57, 66, 58, 76, 59, 86, 96, 94, 75, 77, 87, 73, 94, 79, 88, 82, 92, 98, 84, 93, 90, 99, 97, 91, 95] # 2x2 blocks
toIndex = [12, 14, 16, 17, 20, 22, 23, 26, 28, 30, 32, 33, 34, 36, 37, 39, 42, 44, 45, 48, 50, 51, 53, 54, 55, 56, 57, 58, 59, 64, 66, 67, 70, 72, 73, 75, 76, 77, 78, 79, 82, 84, 85, 86, 88, 89, 90]
for i in range(len(fromIndex)):
    coproductRepresentation = coproductRepresentation.elementary_row_op(op = "n<->m", row1 = toIndex[i] - 1, row2 = fromIndex[i] - 1)
    coproductRepresentation = coproductRepresentation.elementary_col_op(op = "n<->m", col1 = toIndex[i] - 1, col2 = fromIndex[i] - 1)


#print(coproductRepresentation.shape, file = fout)
#print("\n\n\n", file = fout)
'''simplified = [["0" for i in range((2*n)**2)] for j in range((2*n)**2)]
count = 1
for i in range((2*n)**2):
    for j in range((2*n)**2):
        print(count, i, j)
        if (coproductRepresentation[i,j] != 0):
            simplified[i][j] = "*"
        count += 1
for i in range((2*n)**2):     
    str = ""
    for j in range((2*n)**2):
        str += simplified[i][j]
    print(str, file = fout)'''

#for i in range(90, 100): # print out 1-blocks
#    print(simplify(simplify(coproductRepresentation[i,i]).subs(r, q - 1/q)))

tenByTen = zeros(10) # print out 10-block
count = 0
for i in range(10):
    for j in range(10):
        print(count)
        tenByTen[i,j] = simplify(simplify(coproductRepresentation[i,j]).subs(r, q - 1/q))
        count += 1
#print(tenByTen)
'''twoByTwo = Matrix([[2*q**10 - q**8 + q**6 + q**4 + q**2 + 2 + q**(-2) + q**(-4) + 2*q**(-8) - q**(-10) + q**(-12), q**11 - 2*q**9 + q**7 + q**(-7) - 2*q**(-9) + q**(-11)], [q**11 - 2*q**9 + q**7 + q**(-7) - 2*q**(-9) + q**(-11), 2*q**(-10) - q**(-8) + q**(-6) + q**(-4) + q**(-2) + 2 + q**(2) + q**(4) + 2*q**(8) - q**(10) + q**(12)]])

const = q**12 + q**6 + q**4 + q**2 + 2 + q**(-2) + q**(-4) + q**(-6) + q**(-12)
H = twoByTwo - const * eye(2)
eVect = H.nullspace()[0]
G = zeros(2)
for i in range(2):
    G[i,i] = eVect[i,0]
print(G)'''

# Basis (v1 ⊗ v6, v2 ⊗  v7, v3 ⊗  v8, v4 ⊗  v9, v5 ⊗  v10, v6 ⊗  v1, v7 ⊗  v2, v8 ⊗  v3, v9 ⊗  v4, v10 ⊗  v5)
#tenByTen = Matrix([[(q**20 - q**18 + 2*q**16 + 4*q**10 - q**8 + 2*q**6 + q**4 - q**2 + 3)/q**8, -q**5 + 3*q**3 - 3*q + 1/q - 2/q**5 + 4/q**7 - 2/q**9, -q**6 + 3*q**4 - 3*q**2 + 1 - 2/q**4 + 4/q**6 - 2/q**8, -q**7 + 3*q**5 - 3*q**3 + q - 2/q**3 + 4/q**5 - 2/q**7, (-q**14 + 3*q**12 - 3*q**10 + q**8 - 2*q**4 + 4*q**2 - 2)/q**6, (q**2 - 1)**4*(q**6 + q**4 + q**2 + 1)**2/q**10, -q**11 + 3*q**9 - 3*q**7 + q**5 - 2*q + 4/q - 2/q**3, -q**10 + 3*q**8 - 3*q**6 + q**4 - 2 + 4/q**2 - 2/q**4, -q**9 + 3*q**7 - 3*q**5 + q**3 - 2/q + 4/q**3 - 2/q**5, (-q**14 + 3*q**12 - 3*q**10 + q**8 - 2*q**4 + 4*q**2 - 2)/q**6], [-q**5 + 3*q**3 - 3*q + 1/q - 2/q**5 + 4/q**7 - 2/q**9, (q**22 - q**20 + 2*q**18 - q**16 + 4*q**14 - 2*q**12 + 3*q**10 + q**8 - q**6 + 5*q**4 - 3*q**2 + 2)/q**10, -q**7 + 3*q**5 - 3*q**3 + q - 2/q**3 + 4/q**5 - 2/q**7, (-q**14 + 3*q**12 - 3*q**10 + q**8 - 2*q**4 + 4*q**2 - 2)/q**6, -q**9 + 3*q**7 - 3*q**5 + q**3 - 2/q + 4/q**3 - 2/q**5, -2*q**3 + 4*q - 2/q + q**(-5) - 3/q**7 + 3/q**9 - 1/q**11, q**10 - 2*q**8 + q**6 - 2*q**2 + 4 - 2/q**2 + q**(-6) - 2/q**8 + q**(-10), -q**11 + 3*q**9 - 3*q**7 + q**5 - 2*q + 4/q - 2/q**3, -q**10 + 3*q**8 - 3*q**6 + q**4 - 2 + 4/q**2 - 2/q**4, -q**9 + 3*q**7 - 3*q**5 + q**3 - 2/q + 4/q**3 - 2/q**5], [-q**6 + 3*q**4 - 3*q**2 + 1 - 2/q**4 + 4/q**6 - 2/q**8, -q**7 + 3*q**5 - 3*q**3 + q - 2/q**3 + 4/q**5 - 2/q**7, (q**22 - q**20 + q**18 + 3*q**16 - 2*q**14 + 2*q**12 + 2*q**10 - q**8 + 5*q**6 - q**4 - q**2 + 2)/q**10, -q**9 + 3*q**7 - 3*q**5 + q**3 - 2/q + 4/q**3 - 2/q**5, -q**10 + 3*q**8 - 3*q**6 + q**4 - 2 + 4/q**2 - 2/q**4, -2*q**4 + 4*q**2 - 2 + q**(-4) - 3/q**6 + 3/q**8 - 1/q**10, -2*q**3 + 4*q - 2/q + q**(-5) - 3/q**7 + 3/q**9 - 1/q**11, q**10 - 2*q**8 + q**6 - 2*q**2 + 4 - 2/q**2 + q**(-6) - 2/q**8 + q**(-10), -q**11 + 3*q**9 - 3*q**7 + q**5 - 2*q + 4/q - 2/q**3, -q**10 + 3*q**8 - 3*q**6 + q**4 - 2 + 4/q**2 - 2/q**4], [-q**7 + 3*q**5 - 3*q**3 + q - 2/q**3 + 4/q**5 - 2/q**7, (-q**14 + 3*q**12 - 3*q**10 + q**8 - 2*q**4 + 4*q**2 - 2)/q**6, -q**9 + 3*q**7 - 3*q**5 + q**3 - 2/q + 4/q**3 - 2/q**5, (q**22 - 2*q**20 + 5*q**18 - 3*q**16 + 2*q**14 + q**12 + 5*q**8 - q**6 + q**4 - q**2 + 2)/q**10, -q**11 + 3*q**9 - 3*q**7 + q**5 - 2*q + 4/q - 2/q**3, (-2*q**14 + 4*q**12 - 2*q**10 + q**6 - 3*q**4 + 3*q**2 - 1)/q**9, -2*q**4 + 4*q**2 - 2 + q**(-4) - 3/q**6 + 3/q**8 - 1/q**10, -2*q**3 + 4*q - 2/q + q**(-5) - 3/q**7 + 3/q**9 - 1/q**11, q**10 - 2*q**8 + q**6 - 2*q**2 + 4 - 2/q**2 + q**(-6) - 2/q**8 + q**(-10), -q**11 + 3*q**9 - 3*q**7 + q**5 - 2*q + 4/q - 2/q**3], [(-q**14 + 3*q**12 - 3*q**10 + q**8 - 2*q**4 + 4*q**2 - 2)/q**6, -q**9 + 3*q**7 - 3*q**5 + q**3 - 2/q + 4/q**3 - 2/q**5, -q**10 + 3*q**8 - 3*q**6 + q**4 - 2 + 4/q**2 - 2/q**4, -q**11 + 3*q**9 - 3*q**7 + q**5 - 2*q + 4/q - 2/q**3, (2*q**20 - q**18 + q**16 + q**14 - q**12 + 6*q**10 - q**8 + q**6 + q**4 - q**2 + 2)/q**10, (-2*q**14 + 4*q**12 - 2*q**10 + q**6 - 3*q**4 + 3*q**2 - 1)/q**8, (-2*q**14 + 4*q**12 - 2*q**10 + q**6 - 3*q**4 + 3*q**2 - 1)/q**9, -2*q**4 + 4*q**2 - 2 + q**(-4) - 3/q**6 + 3/q**8 - 1/q**10, -2*q**3 + 4*q - 2/q + q**(-5) - 3/q**7 + 3/q**9 - 1/q**11, q**10 - 2*q**8 + q**6 - 2*q**2 + 4 - 2/q**2 + q**(-6) - 2/q**8 + q**(-10)], [(q**2 - 1)**4*(q**6 + q**4 + q**2 + 1)**2/q**10, -2*q**3 + 4*q - 2/q + q**(-5) - 3/q**7 + 3/q**9 - 1/q**11, -2*q**4 + 4*q**2 - 2 + q**(-4) - 3/q**6 + 3/q**8 - 1/q**10, (-2*q**14 + 4*q**12 - 2*q**10 + q**6 - 3*q**4 + 3*q**2 - 1)/q**9, (-2*q**14 + 4*q**12 - 2*q**10 + q**6 - 3*q**4 + 3*q**2 - 1)/q**8, (3*q**20 - q**18 + q**16 + 2*q**14 - q**12 + 4*q**10 + 2*q**4 - q**2 + 1)/q**12, -2*q**9 + 4*q**7 - 2*q**5 + q - 3/q + 3/q**3 - 1/q**5, -2*q**8 + 4*q**6 - 2*q**4 + 1 - 3/q**2 + 3/q**4 - 1/q**6, (-2*q**14 + 4*q**12 - 2*q**10 + q**6 - 3*q**4 + 3*q**2 - 1)/q**7, (-2*q**14 + 4*q**12 - 2*q**10 + q**6 - 3*q**4 + 3*q**2 - 1)/q**8], [-q**11 + 3*q**9 - 3*q**7 + q**5 - 2*q + 4/q - 2/q**3, q**10 - 2*q**8 + q**6 - 2*q**2 + 4 - 2/q**2 + q**(-6) - 2/q**8 + q**(-10), -2*q**3 + 4*q - 2/q + q**(-5) - 3/q**7 + 3/q**9 - 1/q**11, -2*q**4 + 4*q**2 - 2 + q**(-4) - 3/q**6 + 3/q**8 - 1/q**10, (-2*q**14 + 4*q**12 - 2*q**10 + q**6 - 3*q**4 + 3*q**2 - 1)/q**9, -2*q**9 + 4*q**7 - 2*q**5 + q - 3/q + 3/q**3 - 1/q**5, (2*q**22 - 3*q**20 + 5*q**18 - q**16 + q**14 + 3*q**12 - 2*q**10 + 4*q**8 - q**6 + 2*q**4 - q**2 + 1)/q**12, (-2*q**14 + 4*q**12 - 2*q**10 + q**6 - 3*q**4 + 3*q**2 - 1)/q**7, (-2*q**14 + 4*q**12 - 2*q**10 + q**6 - 3*q**4 + 3*q**2 - 1)/q**8, (-2*q**14 + 4*q**12 - 2*q**10 + q**6 - 3*q**4 + 3*q**2 - 1)/q**9], [-q**10 + 3*q**8 - 3*q**6 + q**4 - 2 + 4/q**2 - 2/q**4, -q**11 + 3*q**9 - 3*q**7 + q**5 - 2*q + 4/q - 2/q**3, q**10 - 2*q**8 + q**6 - 2*q**2 + 4 - 2/q**2 + q**(-6) - 2/q**8 + q**(-10), -2*q**3 + 4*q - 2/q + q**(-5) - 3/q**7 + 3/q**9 - 1/q**11, -2*q**4 + 4*q**2 - 2 + q**(-4) - 3/q**6 + 3/q**8 - 1/q**10, -2*q**8 + 4*q**6 - 2*q**4 + 1 - 3/q**2 + 3/q**4 - 1/q**6, (-2*q**14 + 4*q**12 - 2*q**10 + q**6 - 3*q**4 + 3*q**2 - 1)/q**7, (2*q**22 - q**20 - q**18 + 5*q**16 - q**14 + 2*q**12 + 2*q**10 - 2*q**8 + 3*q**6 + q**4 - q**2 + 1)/q**12, (-2*q**14 + 4*q**12 - 2*q**10 + q**6 - 3*q**4 + 3*q**2 - 1)/q**9, -2*q**4 + 4*q**2 - 2 + q**(-4) - 3/q**6 + 3/q**8 - 1/q**10], [-q**9 + 3*q**7 - 3*q**5 + q**3 - 2/q + 4/q**3 - 2/q**5, -q**10 + 3*q**8 - 3*q**6 + q**4 - 2 + 4/q**2 - 2/q**4, -q**11 + 3*q**9 - 3*q**7 + q**5 - 2*q + 4/q - 2/q**3, q**10 - 2*q**8 + q**6 - 2*q**2 + 4 - 2/q**2 + q**(-6) - 2/q**8 + q**(-10), -2*q**3 + 4*q - 2/q + q**(-5) - 3/q**7 + 3/q**9 - 1/q**11, (-2*q**14 + 4*q**12 - 2*q**10 + q**6 - 3*q**4 + 3*q**2 - 1)/q**7, (-2*q**14 + 4*q**12 - 2*q**10 + q**6 - 3*q**4 + 3*q**2 - 1)/q**8, (-2*q**14 + 4*q**12 - 2*q**10 + q**6 - 3*q**4 + 3*q**2 - 1)/q**9, (2*q**22 - q**20 + q**18 - q**16 + 5*q**14 + q**10 + 2*q**8 - 3*q**6 + 5*q**4 - 2*q**2 + 1)/q**12, -2*q**3 + 4*q - 2/q + q**(-5) - 3/q**7 + 3/q**9 - 1/q**11], [(-q**14 + 3*q**12 - 3*q**10 + q**8 - 2*q**4 + 4*q**2 - 2)/q**6, -q**9 + 3*q**7 - 3*q**5 + q**3 - 2/q + 4/q**3 - 2/q**5, -q**10 + 3*q**8 - 3*q**6 + q**4 - 2 + 4/q**2 - 2/q**4, -q**11 + 3*q**9 - 3*q**7 + q**5 - 2*q + 4/q - 2/q**3, q**10 - 2*q**8 + q**6 - 2*q**2 + 4 - 2/q**2 + q**(-6) - 2/q**8 + q**(-10), (-2*q**14 + 4*q**12 - 2*q**10 + q**6 - 3*q**4 + 3*q**2 - 1)/q**8, (-2*q**14 + 4*q**12 - 2*q**10 + q**6 - 3*q**4 + 3*q**2 - 1)/q**9, -2*q**4 + 4*q**2 - 2 + q**(-4) - 3/q**6 + 3/q**8 - 1/q**10, -2*q**3 + 4*q - 2/q + q**(-5) - 3/q**7 + 3/q**9 - 1/q**11, (2*q**20 - q**18 + q**16 + q**14 - q**12 + 6*q**10 - q**8 + q**6 + q**4 - q**2 + 2)/q**10]])
'''signs = [["?" for i in range(10)] for j in range(10)] # signs of 10-block
for i in range(10):
    for j in range(10):
        for k in range(1, 10):
            if (tenByTen[i,j].subs(q, 10 / k) < 0):
                signs[i][j] = "-"
for i in range(10):     
    str = ""
    for j in range(10):
        str += signs[i][j]
    print(str, file = fout)'''
const = simplify(coproductRepresentation[99, 99].subs(r, q - 1/q))
#const = q**12 + q**6 + q**4 + q**2 + 2 + q**(-2) + q**(-4) + q**(-6) + q**(-12)
H = tenByTen - const * eye(10)
'''HNum = zeros(10)
for i in range(10):
    for j in range(10):
        HNum[i,j] = H[i,j].subs(q, 2)

print(HNum.nullspace())'''


z = Symbol('z') # zero
eVects = [[-q**2, q, z, z, z, -1, q, z, z, z],
          [z, -q**2, q, z, z, z, -1, q, z, z],
          [z, z, -q**2, q, z, z, z, -1, q, z],
          [z, z, z, -q**2, q, z, z, z, -1, q]]
conjs = [zeros(10) for i in range(len(eVects))]
for k in range(len(conjs)):
    for i in range(10):
        for j in range(10):
            conjs[k][i, j] = simplify((H[i, j] * eVects[k][j] / eVects[k][i]).subs(z, 0))


'''for k in range(len(conjs)):
    for i in range(10):     
        str = ""
        for j in range(10):
            eval = conjs[k][i,j].subs(r, q - 1/q).subs(q, 0.5)
            dirichlet = ""
            if (eval == zoo):
                dirichlet = "&"
            elif (eval == 0):
                dirichlet = "0"
            elif (eval < 0):
                dirichlet = "-"
            else:
                dirichlet = "?"
            str += dirichlet
        print(str, file = fout)
    print("\n\n\n\n", file = fout)'''

conjSubmatrices = [zeros(4) for i in range(len(conjs))]
# First conjugate matrix: remove rows/cols 3, 4, 5, 8, 9, 10 (1-indexed)
# Basis (v1 ⊗ v6, v2 ⊗ v7, v6 ⊗ v1, v7 ⊗ v2)
conjSubmatrices[0] = conjs[0][(0, 1, 5, 6), (0, 1, 5, 6)]

# Second conjugate matrix: remove rows/cols 1, 4, 5, 6, 9, 10 (1-indexed)
# Basis (v2 ⊗ v7, v3 ⊗ v8, v7 ⊗ v2, v8 ⊗ v3)
conjSubmatrices[1] = conjs[1][(1, 2, 6, 7), (1, 2, 6, 7)]
 
# Third conjugate matrix: remove rows/cols 1, 2, 5, 6, 7, 10 (1-indexed)
# Basis (v3 ⊗ v8, v4 ⊗ v9, v8 ⊗ v3, v9 ⊗ v4)
conjSubmatrices[2] = conjs[2][(2, 3, 7, 8), (2, 3, 7, 8)]

# Fourth conjugate matrix: remove rows/cols 1, 2, 3, 6, 7, 8 (1-indexed)
# Basis (v4 ⊗ v9, v5 ⊗ v10, v9 ⊗ v4, v10 ⊗ v5)
conjSubmatrices[3] = conjs[3][(3, 4, 8, 9), (3, 4, 8, 9)]

for k in range(4):
    print(f"Conjugated Hamiltonian {k+1}", file = fout)
    for i in range(4):
        print("|", file = fout, end = "\t")
        for j in range(4):
            print(f"{simplify(conjSubmatrices[k][i,j])}", file = fout, end = ",\t")
        print("|", file = fout)
print("\n")



# δ

############ Type D ASEP ###############

# Jump rates:
def rSpeed(q, n):
    return q**(2*n-1)+ q*(-2*n+1)
def rSwap(q, n):
    return (q**(n-1) - q**(-n+1))**2
def rTwist(q, n):
    return 2*q**2 - q**(-2*(n-2)) + q**(-2*(n-1))
def rStick(q, n):
    return q**(2*n) - q**(2*(n-1)) + 2
def rTwistRight(q,n,δ):
    return (1+q**(-2*δ))*rTwist(1/q,n)
def rTwistLeft(q, n, δ):
    return (1+q**(2*δ))*rTwist(q,n)
# Generator matrix with two sites
# Basis: ((3, 0), (2, 1), (0, 3), (1, 2))
# |     *           Twist-Right     Swap        Twist-Right |
# | Stick-Left          *       Stick-Right     Swap        |
# | Swap            Twist-Left      *           Twist-Left  |
# | Stick-Left      Swap        Stick-Right         *       |
def generator(q, n, δ):
    gen = zeros(4)
    gen[0, 1] = rTwistRight(q, n ,δ) * 1 / (1 + q**(2*δ))
    gen[0, 2] = q**(-2) * rSwap(q, n)
    gen[0, 3] = rTwistRight(q, n ,δ) * q**(2*δ) / (1 + q**(2*δ))
    gen[0, 0] = - (gen[0, 1] + gen[0, 2] + gen[0, 3])
    gen[1, 0] = q**(-2*δ) * rStick(q, n) 
    gen[1, 2] = rStick(q**(-1), n)
    gen[1, 3] = rSwap(q, n)
    gen[1, 1] = - (gen[1, 0] + gen[1, 2] + gen[1, 3])
    gen[2, 0] = q**2 * rSwap(q, n)
    gen[2, 1] = rTwistLeft(q, n, δ) * 1 / (1 + q**(2*δ)) # SWITCHED?
    gen[2, 3] = rTwistLeft(q, n, δ) * q**(2*δ) / (1 + q**(2*δ)) # SWITCHED?
    gen[2, 2] = - (gen[2, 0] + gen[2, 1] + gen[2, 3])
    gen[3, 0] = rStick(q, n)
    gen[3, 1] = rSwap(q, n)
    gen[3, 2] = q**(2*δ) * rStick(1/q, n)
    gen[3, 3] = - (gen[3, 0] + gen[3, 1] + gen[3, 2])   
    return gen

gens = [generator(q, 5, i) for i in range(4)]

for k in range(4):
    print(f"Generator {k+1}", file = fout)
    for i in range(4):
        print("|", file = fout, end = "\t")
        for j in range(4):
            print(f"{simplify(gens[k][i,j])}", file = fout, end = ",\t")
        print("|", file = fout)
print("\n")

'''for k in range(len(conjs)):
    for i in range(4):     
        str = ""
        for j in range(4):
            eval = conjSubmatrices[k][i,j].subs(r, q - 1/q).subs(q, 0.5)
            dirichlet = ""
            if (eval == zoo):
                dirichlet = "&"
            elif (eval == 0):
                dirichlet = "0"
            elif (eval < 0):
                dirichlet = "-"
            else:
                dirichlet = "?"
            str += dirichlet
        print(str, file = fout)
    print("\n\n\n\n", file = fout)

for k in range(len(conjs)):    
    print(conjSubmatrices[k], file = fout)
    for i in range(4):
        sum = 0
        for j in range(4):
            sum += conjSubmatrices[k][i,j]
        print(simplify(sum), file = fout)
    print("\n\n\n\n", file = fout)'''

        

'''for i in range(2*n):
    print(simplify(representation[i,i].subs(r, q - 1/q)))'''


'''final = zeros(2*n)
for i in range(2*n):
    for j in range(2*n):
         final[i,j] = collect(powsimp(powdenest(expand(representation[i,j].subs(r, q - 1/q)))), q)
print(final)
for i in range(2*n):
    print(simplify(final[i,i]))'''

    
end = time.time()

print(end-start)  





