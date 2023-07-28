from sympy import *
from sympy.physics.quantum import TensorProduct


######## Variables ########
n = 5
q = Symbol('q')
r = Symbol('r')

######## Matrices #########
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



########## Use for fundamental representation ############
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



########### Case 0 ##########
term1 = q**(-8)*K4*KNeg5*KNeg1**2*KNeg2**2*KNeg3**2*KNeg4**2
term2 = q**(-6)*K4*KNeg5*KNeg2**2*KNeg3**2*KNeg4**2
term3 = q**(-4)*K4*KNeg5*KNeg3**2*KNeg4**2
term4 = q**(-2)*K4*KNeg5*KNeg4**2
term5 = K4*KNeg5
termNeg1 = q**8*KNeg4*K5*K1**2*K2**2*K3**2*K4**2
termNeg2 = q**6*KNeg4*K5*K2**2*K3**2*K4**2
termNeg3 = q**4*KNeg4*K5*K3**2*K4**2
termNeg4 = q**2*KNeg4*K5*K4**2
termNeg5 = KNeg4*K5

########### Case 1 #########
# n, i, j, case
term5121 = ((q - 1/q)**2/q**7)*F1*K4*KNeg5*KNeg1**2*KNeg2**2*KNeg3**2*KNeg4**2*K1*E1
term5131 = ((q - 1/q)**4/q**7)*((-q/(q**2 - 1))*F2*F1 + (q**2/(q**2 - 1))*F1*F2)*K4*KNeg5*KNeg1**2*KNeg2**2*KNeg3**2*KNeg4**2*K1*K2*((-q/(q**2 - 1))*E1*E2 + (q**2/(q**2 - 1))*E2*E1)
term5141 = ((q - 1/q)**6/q**7)*((q**2/(q**4 - 2*q**2 + 1))*F3*F2*F1 + (-q**3/(q**4 - 2*q**2 + 1))*F1*F3*F2 + (-q**3/(q**4 - 2*q**2 + 1))*F2*F1*F3 + (q**4/(q**4 - 2*q**2 + 1))*F1*F2*F3)*K4*KNeg5*KNeg1**2*KNeg2**2*KNeg3**2*KNeg4**2*K1*K2*K3*((q**2/(q**4 - 2*q**2 + 1))*E1*E2*E3 + (-q**3/(q**4 - 2*q**2 + 1))*E2*E3*E1 + (-q**3/(q**4 - 2*q**2 + 1))*E3*E1*E2 + (q**4/(q**4 - 2*q**2 + 1))*E3*E2*E1)
term5151 = ((q - 1/q)**8/q**7)*((-q**3/(q**6 - 3*q**4 + 3*q**2 - 1))*F4*F3*F2*F1 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*F1*F4*F3*F2 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*F2*F1*F4*F3 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*F3*F2*F1*F4 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*F1*F2*F4*F3 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*F1*F3*F2*F4 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*F2*F1*F3*F4 + (q**6/(q**6 - 3*q**4 + 3*q**2 - 1))*F1*F2*F3*F4)*K4*KNeg5*KNeg1**2*KNeg2**2*KNeg3**2*KNeg4**2*K1*K2*K3*K4*((-q**3/(q**6 - 3*q**4 + 3*q**2 - 1))*E1*E2*E3*E4 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*E2*E3*E4*E1 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*E3*E4*E1*E2 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*E4*E1*E2*E3 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*E3*E4*E2*E1 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*E4*E2*E3*E1 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*E4*E3*E1*E2 + (q**6/(q**6 - 3*q**4 + 3*q**2 - 1))*E4*E3*E2*E1)
term5231 = ((q - 1/q)**2/q**5)*F2*K4*KNeg5*KNeg2**2*KNeg3**2*KNeg4**2*K2*E2
term5241 = ((q - 1/q)**4/q**5)*((-q/(q**2 - 1))*F3*F2 + (q**2/(q**2 - 1))*F2*F3)*K4*KNeg5*KNeg2**2*KNeg3**2*KNeg4**2*K2*K3*((-q/(q**2 - 1))*E2*E3 + (q**2/(q**2 - 1))*E3*E2)
term5251 = ((q - 1/q)**6/q**5)*((q**2/(q**4 - 2*q**2 + 1))*F4*F3*F2 + (-q**3/(q**4 - 2*q**2 + 1))*F2*F4*F3 + (-q**3/(q**4 - 2*q**2 + 1))*F3*F2*F4 + (q**4/(q**4 - 2*q**2 + 1))*F2*F3*F4)*K4*KNeg5*KNeg2**2*KNeg3**2*KNeg4**2*K2*K3*K4*((q**2/(q**4 - 2*q**2 + 1))*E2*E3*E4 + (-q**3/(q**4 - 2*q**2 + 1))*E3*E4*E2 + (-q**3/(q**4 - 2*q**2 + 1))*E4*E2*E3 + (q**4/(q**4 - 2*q**2 + 1))*E4*E3*E2)
term5341 = ((q - 1/q)**2/q**3)*F3*K4*KNeg5*KNeg3**2*KNeg4**2*K3*E3
term5351 = ((q - 1/q)**4/q**3)*((-q/(q**2 - 1))*F4*F3 + (q**2/(q**2 - 1))*F3*F4)*K4*KNeg5*KNeg3**2*KNeg4**2*K3*K4*((-q/(q**2 - 1))*E3*E4 + (q**2/(q**2 - 1))*E4*E3)
term5451 = ((q - 1/q)**2/q)*F4*K4*KNeg5*KNeg4**2*K4*E4
########### Case 2 ##########
term5142 = (-(q - 1/q)**10/q**7)*((q**4/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F5*F4*F3*F2*F1 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F1*F5*F4*F3*F2 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F2*F1*F5*F4*F3 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F4*F3*F2*F1*F5 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F5*F3*F2*F1*F4 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F1*F2*F5*F4*F3 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F1*F4*F3*F2*F5 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F1*F5*F3*F2*F4 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F2*F1*F4*F3*F5 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F2*F1*F5*F3*F4 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F3*F2*F1*F5*F4 + (-q**7/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F1*F2*F4*F3*F5 + (-q**7/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F1*F2*F5*F3*F4 + (-q**7/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F1*F3*F2*F5*F4 + (-q**7/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F2*F1*F3*F5*F4 + (q**8/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F1*F2*F3*F5*F4)*KNeg1*KNeg2*KNeg3*((q**4/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E1*E2*E3*E4*E5 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E2*E3*E4*E5*E1 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E3*E4*E5*E1*E2 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E4*E1*E2*E3*E5 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E5*E1*E2*E3*E4 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E3*E4*E5*E2*E1 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E4*E2*E3*E5*E1 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E4*E3*E5*E1*E2 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E4*E5*E1*E2*E3 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E5*E2*E3*E4*E1 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E5*E3*E4*E1*E2 + (-q**7/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E4*E3*E5*E2*E1 + (-q**7/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E4*E5*E2*E3*E1 + (-q**7/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E4*E5*E3*E1*E2 + (-q**7/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E5*E3*E4*E2*E1 + (q**8/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E4*E5*E3*E2*E1)
term5152 = (-(q - 1/q)**8/q**7)*((-q**3/(q**6 - 3*q**4 + 3*q**2 - 1))*F5*F3*F2*F1 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*F1*F5*F3*F2 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*F2*F1*F5*F3 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*F3*F2*F1*F5 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*F1*F2*F5*F3 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*F1*F3*F2*F5 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*F2*F1*F3*F5 + (q**6/(q**6 - 3*q**4 + 3*q**2 - 1))*F1*F2*F3*F5)*KNeg1*KNeg2*KNeg3*KNeg4*((-q**3/(q**6 - 3*q**4 + 3*q**2 - 1))*E1*E2*E3*E5 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*E2*E3*E5*E1 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*E3*E5*E1*E2 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*E5*E1*E2*E3 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*E3*E5*E2*E1 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*E5*E2*E3*E1 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*E5*E3*E1*E2 + (q**6/(q**6 - 3*q**4 + 3*q**2 - 1))*E5*E3*E2*E1)
term5222 = (-(q - 1/q)**12/q**4)*((-q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F2**2*F3*F5*F4*F3 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F2*F3*F5*F2*F3*F4 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F2*F4*F3*F5*F2*F3 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F3*F2*F4*F3*F5*F2 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F5*F3*F2*F4*F3*F2 + ((-q**7 - q**5)/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F2*F3*F4*F3*F5*F2 + ((-q**7 - q**5)/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F2*F5*F3*F4*F3*F2 + ((-q**7 - q**5)/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F3*F2*F5*F4*F3*F2 + ((q**8 + q**6 + q**4)/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F2*F3*F5*F4*F3*F2)*((-q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E2**2*E3*E5*E4*E3 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E2*E3*E5*E2*E3*E4 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E2*E4*E3*E5*E2*E3 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E3*E2*E4*E3*E5*E2 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E5*E3*E2*E4*E3*E2 + ((-q**7 - q**5)/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E2*E3*E4*E3*E5*E2 + ((-q**7 - q**5)/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E2*E5*E3*E4*E3*E2 + ((-q**7 - q**5)/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E3*E2*E5*E4*E3*E2 + ((q**8 + q**6 + q**4)/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E2*E3*E5*E4*E3*E2)
term5232 = (-(q - 1/q)**10/q**5)*((q**4/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F3*F5*F4*F3*F2 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F3*F4*F3*F2*F5 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F5*F3*F4*F3*F2 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F3*F2*F4*F3*F5 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F3*F2*F5*F3*F4 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F4*F3*F2*F5*F3 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F5*F3*F2*F4*F3 + (-q**7/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F2*F3*F4*F3*F5 + (-q**7/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F2*F5*F3*F4*F3 + (q**8/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F2*F3*F5*F4*F3 + ((-q**7 - q**5)/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F3*F2*F5*F4*F3)*KNeg2*((q**4/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E2*E3*E4*E5*E3 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*E3*E5*E3*E4*E2 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*E4*E3*E5*E3*E2 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E3*E2*E3*E4*E5 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E4*E3*E2*E3*E5 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E4*E5*E3*E2*E3 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E5*E3*E2*E3*E4 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E3*E4*E2*E3*E5 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E3*E5*E2*E3*E4 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E4*E3*E5*E2*E3 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E5*E3*E4*E2*E3 + (q**8/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E3*E4*E5*E3*E2 + ((-q**7 - q**5)/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E3*E4*E5*E2*E3)
term5242 = (-(q - 1/q)**8/q**5)*((-q**3/(q**6 - 3*q**4 + 3*q**2 - 1))*F5*F4*F3*F2 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*F2*F5*F4*F3 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*F4*F3*F2*F5 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*F5*F3*F2*F4 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*F2*F4*F3*F5 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*F2*F5*F3*F4 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*F3*F2*F5*F4 + (q**6/(q**6 - 3*q**4 + 3*q**2 - 1))*F2*F3*F5*F4)*KNeg2*KNeg3*((-q**3/(q**6 - 3*q**4 + 3*q**2 - 1))*E2*E3*E4*E5 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*E3*E4*E5*E2 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*E4*E2*E3*E5 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*E5*E2*E3*E4 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*E4*E3*E5*E2 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*E4*E5*E2*E3 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*E5*E3*E4*E2 + (q**6/(q**6 - 3*q**4 + 3*q**2 - 1))*E4*E5*E3*E2)
term5252 = (-(q - 1/q)**6/q**5)*((q**2/(q**4 - 2*q**2 + 1))*F5*F3*F2 + (-q**3/(q**4 - 2*q**2 + 1))*F2*F5*F3 + (-q**3/(q**4 - 2*q**2 + 1))*F3*F2*F5 + (q**4/(q**4 - 2*q**2 + 1))*F2*F3*F5)*KNeg2*KNeg3*KNeg4*((q**2/(q**4 - 2*q**2 + 1))*E2*E3*E5 + (-q**3/(q**4 - 2*q**2 + 1))*E3*E5*E2 + (-q**3/(q**4 - 2*q**2 + 1))*E5*E2*E3 + (q**4/(q**4 - 2*q**2 + 1))*E5*E3*E2)
term5322 = (-(q - 1/q)**10/q**3)*((q**4/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F2*F3*F4*F5*F3 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*F3*F5*F3*F4*F2 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*F4*F3*F5*F3*F2 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F3*F2*F3*F4*F5 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F4*F3*F2*F3*F5 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F4*F5*F3*F2*F3 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F5*F3*F2*F3*F4 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F3*F4*F2*F3*F5 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F3*F5*F2*F3*F4 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F4*F3*F5*F2*F3 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F5*F3*F4*F2*F3 + (q**8/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F3*F4*F5*F3*F2 + ((-q**7 - q**5)/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F3*F4*F5*F2*F3)*K2*((q**4/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E3*E5*E4*E3*E2 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E3*E4*E3*E2*E5 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E5*E3*E4*E3*E2 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E3*E2*E4*E3*E5 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E3*E2*E5*E3*E4 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E4*E3*E2*E5*E3 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E5*E3*E2*E4*E3 + (-q**7/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E2*E3*E4*E3*E5 + (-q**7/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E2*E5*E3*E4*E3 + (q**8/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E2*E3*E5*E4*E3 + ((-q**7 - q**5)/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E3*E2*E5*E4*E3)
term5332 = (-(q - 1/q)**8/q**2)*((-q**3/(q**4 - 2*q**2 + 1))*F3*F4*F3*F5 + (-q**3/(q**4 - 2*q**2 + 1))*F5*F3*F4*F3 + ((q**4 + q**2)/(q**4 - 2*q**2 + 1))*F3*F5*F4*F3)*((-q**3/(q**4 - 2*q**2 + 1))*E3*E4*E3*E5 + (-q**3/(q**4 - 2*q**2 + 1))*E5*E3*E4*E3 + ((q**4 + q**2)/(q**4 - 2*q**2 + 1))*E3*E5*E4*E3)
term5342 = (-(q - 1/q)**6/q**3)*((q**2/(q**4 - 2*q**2 + 1))*F5*F4*F3 + (-q**3/(q**4 - 2*q**2 + 1))*F4*F3*F5 + (-q**3/(q**4 - 2*q**2 + 1))*F5*F3*F4 + (q**4/(q**4 - 2*q**2 + 1))*F3*F5*F4)*KNeg3*((q**2/(q**4 - 2*q**2 + 1))*E3*E4*E5 + (-q**3/(q**4 - 2*q**2 + 1))*E4*E3*E5 + (-q**3/(q**4 - 2*q**2 + 1))*E5*E3*E4 + (q**4/(q**4 - 2*q**2 + 1))*E4*E5*E3)
term5352 = (-(q - 1/q)**4/q**3)*((-q/(q**2 - 1))*F5*F3 + (q**2/(q**2 - 1))*F3*F5)*KNeg3*KNeg4*((-q/(q**2 - 1))*E3*E5 + (q**2/(q**2 - 1))*E5*E3)
term5412 = (-(q - 1/q)**10/q)*((q**4/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F1*F2*F3*F4*F5 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F2*F3*F4*F5*F1 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F3*F4*F5*F1*F2 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F4*F1*F2*F3*F5 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F5*F1*F2*F3*F4 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F3*F4*F5*F2*F1 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F4*F2*F3*F5*F1 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F4*F3*F5*F1*F2 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F4*F5*F1*F2*F3 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F5*F2*F3*F4*F1 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F5*F3*F4*F1*F2 + (-q**7/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F4*F3*F5*F2*F1 + (-q**7/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F4*F5*F2*F3*F1 + (-q**7/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F4*F5*F3*F1*F2 + (-q**7/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F5*F3*F4*F2*F1 + (q**8/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*F4*F5*F3*F2*F1)*K1*K2*K3*((q**4/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E5*E4*E3*E2*E1 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E1*E5*E4*E3*E2 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E2*E1*E5*E4*E3 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E4*E3*E2*E1*E5 + (-q**5/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E5*E3*E2*E1*E4 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E1*E2*E5*E4*E3 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E1*E4*E3*E2*E5 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E1*E5*E3*E2*E4 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E2*E1*E4*E3*E5 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E2*E1*E5*E3*E4 + (q**6/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E3*E2*E1*E5*E4 + (-q**7/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E1*E2*E4*E3*E5 + (-q**7/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E1*E2*E5*E3*E4 + (-q**7/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E1*E3*E2*E5*E4 + (-q**7/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E2*E1*E3*E5*E4 + (q**8/(q**8 - 4*q**6 + 6*q**4 - 4*q**2 + 1))*E1*E2*E3*E5*E4)
term5422 = (-(q - 1/q)**8/q)*((-q**3/(q**6 - 3*q**4 + 3*q**2 - 1))*F2*F3*F4*F5 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*F3*F4*F5*F2 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*F4*F2*F3*F5 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*F5*F2*F3*F4 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*F4*F3*F5*F2 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*F4*F5*F2*F3 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*F5*F3*F4*F2 + (q**6/(q**6 - 3*q**4 + 3*q**2 - 1))*F4*F5*F3*F2)*K2*K3*((-q**3/(q**6 - 3*q**4 + 3*q**2 - 1))*E5*E4*E3*E2 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*E2*E5*E4*E3 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*E4*E3*E2*E5 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*E5*E3*E2*E4 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*E2*E4*E3*E5 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*E2*E5*E3*E4 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*E3*E2*E5*E4 + (q**6/(q**6 - 3*q**4 + 3*q**2 - 1))*E2*E3*E5*E4)
term5432 = (-(q - 1/q)**6/q)*((q**2/(q**4 - 2*q**2 + 1))*F3*F4*F5 + (-q**3/(q**4 - 2*q**2 + 1))*F4*F3*F5 + (-q**3/(q**4 - 2*q**2 + 1))*F5*F3*F4 + (q**4/(q**4 - 2*q**2 + 1))*F4*F5*F3)*K3*((q**2/(q**4 - 2*q**2 + 1))*E5*E4*E3 + (-q**3/(q**4 - 2*q**2 + 1))*E4*E3*E5 + (-q**3/(q**4 - 2*q**2 + 1))*E5*E3*E4 + (q**4/(q**4 - 2*q**2 + 1))*E3*E5*E4)
term5442 = (-(q - 1/q)**4)*F5*F4*E5*E4
term5452 = (-(q - 1/q)**2/q)*F5*KNeg4*E5
term5512 = (-q*(q - 1/q)**8)*((-q**3/(q**6 - 3*q**4 + 3*q**2 - 1))*F1*F2*F3*F5 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*F2*F3*F5*F1 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*F3*F5*F1*F2 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*F5*F1*F2*F3 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*F3*F5*F2*F1 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*F5*F2*F3*F1 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*F5*F3*F1*F2 + (q**6/(q**6 - 3*q**4 + 3*q**2 - 1))*F5*F3*F2*F1)*K1*K2*K3*K4*((-q**3/(q**6 - 3*q**4 + 3*q**2 - 1))*E5*E3*E2*E1 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*E1*E5*E3*E2 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*E2*E1*E5*E3 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*E3*E2*E1*E5 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*E1*E2*E5*E3 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*E1*E3*E2*E5 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*E2*E1*E3*E5 + (q**6/(q**6 - 3*q**4 + 3*q**2 - 1))*E1*E2*E3*E5)
term5522 = (-q*(q - 1/q)**6)*((q**2/(q**4 - 2*q**2 + 1))*F2*F3*F5 + (-q**3/(q**4 - 2*q**2 + 1))*F3*F5*F2 + (-q**3/(q**4 - 2*q**2 + 1))*F5*F2*F3 + (q**4/(q**4 - 2*q**2 + 1))*F5*F3*F2)*K2*K3*K4*((q**2/(q**4 - 2*q**2 + 1))*E5*E3*E2 + (-q**3/(q**4 - 2*q**2 + 1))*E2*E5*E3 + (-q**3/(q**4 - 2*q**2 + 1))*E3*E2*E5 + (q**4/(q**4 - 2*q**2 + 1))*E2*E3*E5)
term5532 = (-q*(q - 1/q)**4)*((-q/(q**2 - 1))*F3*F5 + (q**2/(q**2 - 1))*F5*F3)*K3*K4*((-q/(q**2 - 1))*E5*E3 + (q**2/(q**2 - 1))*E3*E5)
term5542 = (-q*(q - 1/q)**2)*F5*K4*E5
########### Case 3 ##########
term5213 = (q**7*(q - 1/q)**2)*F1*KNeg4*K5*K2**2*K3**2*K4**2*K1*E1
term5313 = (q**5*(q - 1/q)**4)*((-q/(q**2 - 1))*F1*F2 + (q**2/(q**2 - 1))*F2*F1)*KNeg4*K5*K3**2*K4**2*K1*K2*((-q/(q**2 - 1))*E2*E1 + (q**2/(q**2 - 1))*E1*E2)
term5413 = (q**3*(q - 1/q)**6)*((q**2/(q**4 - 2*q**2 + 1))*F1*F2*F3 + (-q**3/(q**4 - 2*q**2 + 1))*F2*F3*F1 + (-q**3/(q**4 - 2*q**2 + 1))*F3*F1*F2 + (q**4/(q**4 - 2*q**2 + 1))*F3*F2*F1)*KNeg4*K5*K4**2*K1*K2*K3*((q**2/(q**4 - 2*q**2 + 1))*E3*E2*E1 + (-q**3/(q**4 - 2*q**2 + 1))*E1*E3*E2 + (-q**3/(q**4 - 2*q**2 + 1))*E2*E1*E3 + (q**4/(q**4 - 2*q**2 + 1))*E1*E2*E3)
term5513 = (q*(q - 1/q)**8)*((-q**3/(q**6 - 3*q**4 + 3*q**2 - 1))*F1*F2*F3*F4 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*F2*F3*F4*F1 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*F3*F4*F1*F2 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*F4*F1*F2*F3 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*F3*F4*F2*F1 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*F4*F2*F3*F1 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*F4*F3*F1*F2 + (q**6/(q**6 - 3*q**4 + 3*q**2 - 1))*F4*F3*F2*F1)*KNeg4*K5*K1*K2*K3*K4*((-q**3/(q**6 - 3*q**4 + 3*q**2 - 1))*E4*E3*E2*E1 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*E1*E4*E3*E2 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*E2*E1*E4*E3 + (q**4/(q**6 - 3*q**4 + 3*q**2 - 1))*E3*E2*E1*E4 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*E1*E2*E4*E3 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*E1*E3*E2*E4 + (-q**5/(q**6 - 3*q**4 + 3*q**2 - 1))*E2*E1*E3*E4 + (q**6/(q**6 - 3*q**4 + 3*q**2 - 1))*E1*E2*E3*E4)
term5323 = (q**5*(q - 1/q)**2)*F2*KNeg4*K5*K3**2*K4**2*K2*E2
term5423 = (q**3*(q - 1/q)**4)*((-q/(q**2 - 1))*F2*F3 + (q**2/(q**2 - 1))*F3*F2)*KNeg4*K5*K4**2*K2*K3*((-q/(q**2 - 1))*E3*E2 + (q**2/(q**2 - 1))*E2*E3)
term5523 = (q*(q - 1/q)**6)*((q**2/(q**4 - 2*q**2 + 1))*F2*F3*F4 + (-q**3/(q**4 - 2*q**2 + 1))*F3*F4*F2 + (-q**3/(q**4 - 2*q**2 + 1))*F4*F2*F3 + (q**4/(q**4 - 2*q**2 + 1))*F4*F3*F2)*KNeg4*K5*K2*K3*K4*((q**2/(q**4 - 2*q**2 + 1))*E4*E3*E2 + (-q**3/(q**4 - 2*q**2 + 1))*E2*E4*E3 + (-q**3/(q**4 - 2*q**2 + 1))*E3*E2*E4 + (q**4/(q**4 - 2*q**2 + 1))*E2*E3*E4)
term5433 = (q**3*(q - 1/q)**2)*F3*KNeg4*K5*K4**2*K3*E3
term5533 = (q*(q - 1/q)**4)*((-q/(q**2 - 1))*F3*F4 + (q**2/(q**2 - 1))*F4*F3)*KNeg4*K5*K3*K4*((-q/(q**2 - 1))*E4*E3 + (q**2/(q**2 - 1))*E3*E4)
term5543 = (q*(q - 1/q)**2)*F4*KNeg4*K5*K4*E4
term5112 = - q**(-6) * r**16 * r**(-6) * (-q**3 *F1*F2*F3*F4*F1*F2*F3*F5 - q * (q**2 + 1)**2 *F1*F2*F3*F4*F3*F5*F2*F1 + q**2 * (q**2 + 1) *F1*F2*F3*F5*F2*F3*F4*F1 + q**2 * (q**2 + 1) *F1*F2*F3*F5*F4*F1*F2*F3 + q**3 *F1*F2*F3*F5*F4*F3*F1*F2 + (q**6 + q**4 + q**2 + 1)*F1*F2*F3*F5*F4*F3*F2*F1 + q**2 * (q**2 + 1) *F1*F2*F4*F3*F5*F2*F3*F1 - q**3 *F1*F2*F5*F3*F4*F1*F2*F3 - q * (q**2 + 1)**2 *F1*F2*F5*F3*F4*F3*F2*F1 - q**2 * (q**2 + 1) *F1*F2**2*F3*F5*F4*F3*F1 + q**3 *F1*F3*F2*F1*F2*F3*F5*F4 + q**2 * (q**2 + 1) *F1*F3*F2*F4*F3*F5*F2*F1 - q**3 *F1*F3*F2*F5*F3*F4*F1*F2 - q * (q**2 + 1)**2 *F1*F3*F2*F5*F4*F3*F2*F1 + q**3 *F1*F3*F5*F4*F3*F2*F1*F2 + q**3 *F1*F4*F3*F2*F1*F2*F3*F5 - q**3 *F1*F4*F3*F2*F5*F3*F1*F2 + q**3 *F1*F5*F3*F2*F1*F2*F3*F4 + q**2 * (q**2 + 1) *F1*F5*F3*F2*F4*F3*F2*F1 + q**3 *F1*F5*F4*F3*F2*F1*F2*F3 - q**2 * (q**2 + 1) *F1**2*F2*F3*F5*F4*F3*F2 + q**3 *F2*F1*F2*F3*F5*F4*F3*F1 + q**2 * (q**2 + 1) *F2*F1*F3*F4*F3*F5*F2*F1 - q**3 *F2*F1*F3*F5*F2*F3*F4*F1 - q*(q**4 + q**2 + 1) *F2*F1*F3*F5*F4*F3*F2*F1 - q**3 *F2*F1*F4*F3*F5*F2*F3*F1 + q**2 * (q**2 + 1) *F2*F1*F5*F3*F4*F3*F2*F1 - q**3 *F3*F2*F1*F4*F3*F5*F2*F1 + q**2 * (q**2 + 1) *F3*F2*F1*F5*F4*F3*F2*F1 - q**3 *F5*F3*F2*F1*F4*F3*F2*F1) * r**(-6) * (-q**3 *E1*E2*E3*E4*E1*E2*E3*E5 - q * (q**2 + 1)**2 *E1*E2*E3*E4*E3*E5*E2*E1 + q**2 * (q**2 + 1) *E1*E2*E3*E5*E2*E3*E4*E1 + q**2 * (q**2 + 1) *E1*E2*E3*E5*E4*E1*E2*E3 + q**3 *E1*E2*E3*E5*E4*E3*E1*E2 + (q**6 + q**4 + q**2 + 1)*E1*E2*E3*E5*E4*E3*E2*E1 + q**2 * (q**2 + 1) *E1*E2*E4*E3*E5*E2*E3*E1 - q**3 *E1*E2*E5*E3*E4*E1*E2*E3 - q * (q**2 + 1)**2 *E1*E2*E5*E3*E4*E3*E2*E1 - q**2 * (q**2 + 1) *E1*E2**2*E3*E5*E4*E3*E1 + q**3 *E1*E3*E2*E1*E2*E3*E5*E4 + q**2 * (q**2 + 1) *E1*E3*E2*E4*E3*E5*E2*E1 - q**3 *E1*E3*E2*E5*E3*E4*E1*E2 - q * (q**2 + 1)**2 *E1*E3*E2*E5*E4*E3*E2*E1 + q**3 *E1*E3*E5*E4*E3*E2*E1*E2 + q**3 *E1*E4*E3*E2*E1*E2*E3*E5 - q**3 *E1*E4*E3*E2*E5*E3*E1*E2 + q**3 *E1*E5*E3*E2*E1*E2*E3*E4 + q**2 * (q**2 + 1) *E1*E5*E3*E2*E4*E3*E2*E1 + q**3 *E1*E5*E4*E3*E2*E1*E2*E3 - q**2 * (q**2 + 1) *E1**2*E2*E3*E5*E4*E3*E2 + q**3 *E2*E1*E2*E3*E5*E4*E3*E1 + q**2 * (q**2 + 1) *E2*E1*E3*E4*E3*E5*E2*E1 - q**3 *E2*E1*E3*E5*E2*E3*E4*E1 - q*(q**4 + q**2 + 1) *E2*E1*E3*E5*E4*E3*E2*E1 - q**3 *E2*E1*E4*E3*E5*E2*E3*E1 + q**2 * (q**2 + 1) *E2*E1*E5*E3*E4*E3*E2*E1 - q**3 *E3*E2*E1*E4*E3*E5*E2*E1 + q**2 * (q**2 + 1) *E3*E2*E1*E5*E4*E3*E2*E1 - q**3 *E5*E3*E2*E1*E4*E3*E2*E1)
term5122 = - q**(-7) * r**(14) * r**(-6) * (-q**5*F1*F2*F3*F4*F3*F5*F2 + q**4*F1*F2*F3*F5*F2*F3*F4 - q**3*F1*F2*F3*F5*F4*F2*F3 + q**6*F1*F2*F3*F5*F4*F3*F2 + q**4*F1*F2*F4*F3*F5*F2*F3 - q**5*F1*F2*F5*F3*F4*F3*F2 - (q**4 - q**2)*F1*F2**2*F3*F5*F4*F3 + q**4*F1*F3*F2*F4*F3*F5*F2 - q**5*F1*F3*F2*F5*F4*F3*F2 + q**4*F1*F5*F3*F2*F4*F3*F2 + (q**3 - q)*F2*F1*F2*F3*F5*F4*F3 - q**3*F2*F1*F3*F4*F2*F3*F5 + (q**4 + q**2)*F2*F1*F3*F4*F3*F5*F2 - q**3*F2*F1*F3*F5*F2*F3*F4 + (q**4 + q**2)*F2*F1*F3*F5*F4*F2*F3 - (q**5 + q)*F2*F1*F3*F5*F4*F3*F2 - q**3*F2*F1*F4*F3*F5*F2*F3 - q**3*F2*F1*F5*F3*F4*F2*F3 + (q**4 + q**2)*F2*F1*F5*F3*F4*F3*F2 - q*F2*F3*F4*F3*F5*F2*F1 + q**2*F2*F3*F5*F2*F1*F3*F4 - q**3*F2*F3*F5*F4*F2*F1*F3 + 1*F2*F3*F5*F4*F3*F2*F1 + q**2*F2*F4*F3*F5*F2*F1*F3 - q*F2*F5*F3*F4*F3*F2*F1 - q**3*F3*F2*F1*F4*F3*F5*F2 - q**3*F3*F2*F1*F5*F3*F4*F2 + (q**4 + q**2)*F3*F2*F1*F5*F4*F3*F2 + q**2*F3*F2*F4*F3*F5*F2*F1 - q*F3*F2*F5*F4*F3*F2*F1 - q**3*F4*F3*F2*F1*F5*F3*F2 - q**3*F5*F3*F2*F1*F4*F3*F2 + q**2*F5*F3*F2*F4*F3*F2*F1) * KNeg1 * r**(-6) * (1*E1*E2*E3*E5*E4*E3*E2 - q*E2*E1*E2*E3*E5*E4*E3 + q**2*E2*E3*E1*E2*E3*E5*E4 - q**3*E2*E3*E4*E1*E2*E3*E5 + q**4*E2*E3*E4*E3*E5*E1*E2 - (q**5 - q)*E2*E3*E4*E3*E5*E2*E1 - q**3*E2*E3*E5*E1*E2*E3*E4 + (q**4 - q**2)*E2*E3*E5*E2*E3*E4*E1 + (q**4 + q**2)*E2*E3*E5*E4*E1*E2*E3 - (q**5 + q)*E2*E3*E5*E4*E3*E1*E2 + q**6*E2*E3*E5*E4*E3*E2*E1 + q**2*E2*E4*E3*E1*E2*E3*E5 - q**3*E2*E4*E3*E5*E1*E2*E3 + (q**4 - q**2)*E2*E4*E3*E5*E2*E3*E1 + q**2*E2*E5*E3*E1*E2*E3*E4 - q**3*E2*E5*E3*E4*E1*E2*E3 + q**4*E2*E5*E3*E4*E3*E1*E2 - (q**5 - q)*E2*E5*E3*E4*E3*E2*E1 + q**2*E2*E5*E4*E3*E1*E2*E3 - (q**4 - q**2)*E2**2*E3*E5*E4*E3*E1 - q*E3*E2*E1*E2*E3*E5*E4 + q**2*E3*E2*E4*E1*E2*E3*E5 - q**3*E3*E2*E4*E3*E5*E1*E2 + (q**4 - q**2)*E3*E2*E4*E3*E5*E2*E1 + q**2*E3*E2*E5*E1*E2*E3*E4 - q**3*E3*E2*E5*E3*E4*E1*E2 - (q**3 + q)*E3*E2*E5*E4*E1*E2*E3 + (q**4 + q**2)*E3*E2*E5*E4*E3*E1*E2 - (q**5 - q)*E3*E2*E5*E4*E3*E2*E1 + q**2*E3*E4*E3*E2*E5*E1*E2 - q*E3*E5*E4*E3*E2*E1*E2 - q*E4*E3*E2*E1*E2*E3*E5 + q**2*E4*E3*E2*E5*E1*E2*E3 - q**3*E4*E3*E2*E5*E3*E1*E2 - q*E5*E3*E2*E1*E2*E3*E4 + q**2*E5*E3*E2*E4*E1*E2*E3 - q**3*E5*E3*E2*E4*E3*E1*E2 + (q**4 - q**2)*E5*E3*E2*E4*E3*E2*E1 + q**2*E5*E3*E4*E3*E2*E1*E2 - q*E5*E4*E3*E2*E1*E2*E3)
term5212 = - q**(-5) * r**(14) * r**(-6) * (1*F1*F2*F3*F5*F4*F3*F2 - q*F2*F1*F2*F3*F5*F4*F3 + q**2*F2*F3*F1*F2*F3*F5*F4 - q**3*F2*F3*F4*F1*F2*F3*F5 + q**4*F2*F3*F4*F3*F5*F1*F2 - (q**5 - q)*F2*F3*F4*F3*F5*F2*F1 - q**3*F2*F3*F5*F1*F2*F3*F4 + (q**4 - q**2)*F2*F3*F5*F2*F3*F4*F1 + (q**4 + q**2)*F2*F3*F5*F4*F1*F2*F3 - (q**5 + q)*F2*F3*F5*F4*F3*F1*F2 + q**6*F2*F3*F5*F4*F3*F2*F1 + q**2*F2*F4*F3*F1*F2*F3*F5 - q**3*F2*F4*F3*F5*F1*F2*F3 + (q**4 - q**2)*F2*F4*F3*F5*F2*F3*F1 + q**2*F2*F5*F3*F1*F2*F3*F4 - q**3*F2*F5*F3*F4*F1*F2*F3 + q**4*F2*F5*F3*F4*F3*F1*F2 - (q**5 - q)*F2*F5*F3*F4*F3*F2*F1 + q**2*F2*F5*F4*F3*F1*F2*F3 - (q**4 - q**2)*F2**2*F3*F5*F4*F3*F1 - q*F3*F2*F1*F2*F3*F5*F4 + q**2*F3*F2*F4*F1*F2*F3*F5 - q**3*F3*F2*F4*F3*F5*F1*F2 + (q**4 - q**2)*F3*F2*F4*F3*F5*F2*F1 + q**2*F3*F2*F5*F1*F2*F3*F4 - q**3*F3*F2*F5*F3*F4*F1*F2 - (q**3 + q)*F3*F2*F5*F4*F1*F2*F3 + (q**4 + q**2)*F3*F2*F5*F4*F3*F1*F2 - (q**5 - q)*F3*F2*F5*F4*F3*F2*F1 + q**2*F3*F4*F3*F2*F5*F1*F2 - q*F3*F5*F4*F3*F2*F1*F2 - q*F4*F3*F2*F1*F2*F3*F5 + q**2*F4*F3*F2*F5*F1*F2*F3 - q**3*F4*F3*F2*F5*F3*F1*F2 - q*F5*F3*F2*F1*F2*F3*F4 + q**2*F5*F3*F2*F4*F1*F2*F3 - q**3*F5*F3*F2*F4*F3*F1*F2 + (q**4 - q**2)*F5*F3*F2*F4*F3*F2*F1 + q**2*F5*F3*F4*F3*F2*F1*F2 - q*F5*F4*F3*F2*F1*F2*F3) * K1 * r**(-6) * (-q**5*E1*E2*E3*E4*E3*E5*E2 + q**4*E1*E2*E3*E5*E2*E3*E4 - q**3*E1*E2*E3*E5*E4*E2*E3 + q**6*E1*E2*E3*E5*E4*E3*E2 + q**4*E1*E2*E4*E3*E5*E2*E3 - q**5*E1*E2*E5*E3*E4*E3*E2 - (q**4 - q**2)*E1*E2**2*E3*E5*E4*E3 + q**4*E1*E3*E2*E4*E3*E5*E2 - q**5*E1*E3*E2*E5*E4*E3*E2 + q**4*E1*E5*E3*E2*E4*E3*E2 + (q**3 - q)*E2*E1*E2*E3*E5*E4*E3 - q**3*E2*E1*E3*E4*E2*E3*E5 + (q**4 + q**2)*E2*E1*E3*E4*E3*E5*E2 - q**3*E2*E1*E3*E5*E2*E3*E4 + (q**4 + q**2)*E2*E1*E3*E5*E4*E2*E3 - (q**5 + q)*E2*E1*E3*E5*E4*E3*E2 - q**3*E2*E1*E4*E3*E5*E2*E3 - q**3*E2*E1*E5*E3*E4*E2*E3 + (q**4 + q**2)*E2*E1*E5*E3*E4*E3*E2 - q*E2*E3*E4*E3*E5*E2*E1 + q**2*E2*E3*E5*E2*E1*E3*E4 - q**3*E2*E3*E5*E4*E2*E1*E3 + 1*E2*E3*E5*E4*E3*E2*E1 + q**2*E2*E4*E3*E5*E2*E1*E3 - q*E2*E5*E3*E4*E3*E2*E1 - q**3*E3*E2*E1*E4*E3*E5*E2 - q**3*E3*E2*E1*E5*E3*E4*E2 + (q**4 + q**2)*E3*E2*E1*E5*E4*E3*E2 + q**2*E3*E2*E4*E3*E5*E2*E1 - q*E3*E2*E5*E4*E3*E2*E1 - q**3*E4*E3*E2*E1*E5*E3*E2 - q**3*E5*E3*E2*E1*E4*E3*E2 + q**2*E5*E3*E2*E4*E3*E2*E1)
term5132 = - q**(-7) * r**(12) * r**(-5) * (-q**4*F1*F2*F3*F4*F3*F5 + q**5*F1*F2*F3*F5*F4*F3 - q**4*F1*F2*F5*F3*F4*F3 + q**3*F1*F3*F2*F4*F3*F5 + q**3*F1*F3*F2*F5*F3*F4 - (q**4 + q**2)*F1*F3*F2*F5*F4*F3 - q**2*F1*F3*F4*F3*F2*F5 + q*F1*F3*F5*F4*F3*F2 + q**3*F1*F4*F3*F2*F5*F3 + q**3*F1*F5*F3*F2*F4*F3 - q**2*F1*F5*F3*F4*F3*F2 + q**3*F2*F1*F3*F4*F3*F5 - q**4*F2*F1*F3*F5*F4*F3 + q**3*F2*F1*F5*F3*F4*F3 - q**2*F3*F2*F1*F4*F3*F5 - q**2*F3*F2*F1*F5*F3*F4 + (q**3 + q)*F3*F2*F1*F5*F4*F3 + q*F3*F4*F3*F2*F1*F5 - 1*F3*F5*F4*F3*F2*F1 - q**2*F4*F3*F2*F1*F5*F3 - q**2*F5*F3*F2*F1*F4*F3 + q*F5*F3*F4*F3*F2*F1) * KNeg1 * KNeg2 * r**(-5) * (-1*E1*E2*E3*E5*E4*E3 + q*E2*E3*E5*E4*E3*E1 + q*E3*E1*E2*E3*E5*E4 - q**2*E3*E2*E3*E5*E4*E1 - q**2*E3*E4*E1*E2*E3*E5 + q**3*E3*E4*E2*E3*E5*E1 + (q**3 - q)*E3*E4*E3*E5*E1*E2 - (q**4 - q**2)*E3*E4*E3*E5*E2*E1 - q**2*E3*E5*E1*E2*E3*E4 + q**3*E3*E5*E2*E3*E4*E1 + (q**3 + q)*E3*E5*E4*E1*E2*E3 - (q**4 + q**2)*E3*E5*E4*E2*E3*E1 - q**4*E3*E5*E4*E3*E1*E2 + q**5*E3*E5*E4*E3*E2*E1 + q*E4*E3*E1*E2*E3*E5 - q**2*E4*E3*E2*E3*E5*E1 - q**2*E4*E3*E5*E1*E2*E3 + q**3*E4*E3*E5*E2*E3*E1 + q*E5*E3*E1*E2*E3*E4 - q**2*E5*E3*E2*E3*E4*E1 - q**2*E5*E3*E4*E1*E2*E3 + q**3*E5*E3*E4*E2*E3*E1 + (q**3 - q)*E5*E3*E4*E3*E1*E2 - (q**4 - q**2)*E5*E3*E4*E3*E2*E1 + q*E5*E4*E3*E1*E2*E3 - q**2*E5*E4*E3*E2*E3*E1)
term5312 = - q**(-3) * r**(12) * r**(-5) * (-1*F1*F2*F3*F5*F4*F3 + q*F2*F3*F5*F4*F3*F1 + q*F3*F1*F2*F3*F5*F4 - q**2*F3*F2*F3*F5*F4*F1 - q**2*F3*F4*F1*F2*F3*F5 + q**3*F3*F4*F2*F3*F5*F1 + (q**3 - q)*F3*F4*F3*F5*F1*F2 - (q**4 - q**2)*F3*F4*F3*F5*F2*F1 - q**2*F3*F5*F1*F2*F3*F4 + q**3*F3*F5*F2*F3*F4*F1 + (q**3 + q)*F3*F5*F4*F1*F2*F3 - (q**4 + q**2)*F3*F5*F4*F2*F3*F1 - q**4*F3*F5*F4*F3*F1*F2 + q**5*F3*F5*F4*F3*F2*F1 + q*F4*F3*F1*F2*F3*F5 - q**2*F4*F3*F2*F3*F5*F1 - q**2*F4*F3*F5*F1*F2*F3 + q**3*F4*F3*F5*F2*F3*F1 + q*F5*F3*F1*F2*F3*F4 - q**2*F5*F3*F2*F3*F4*F1 - q**2*F5*F3*F4*F1*F2*F3 + q**3*F5*F3*F4*F2*F3*F1 + (q**3 - q)*F5*F3*F4*F3*F1*F2 - (q**4 - q**2)*F5*F3*F4*F3*F2*F1 + q*F5*F4*F3*F1*F2*F3 - q**2*F5*F4*F3*F2*F3*F1) * K1 * K2 * r**(-5) * (-q**4*E1*E2*E3*E4*E3*E5 + q**5*E1*E2*E3*E5*E4*E3 - q**4*E1*E2*E5*E3*E4*E3 + q**3*E1*E3*E2*E4*E3*E5 + q**3*E1*E3*E2*E5*E3*E4 - (q**4 + q**2)*E1*E3*E2*E5*E4*E3 - q**2*E1*E3*E4*E3*E2*E5 + q*E1*E3*E5*E4*E3*E2 + q**3*E1*E4*E3*E2*E5*E3 + q**3*E1*E5*E3*E2*E4*E3 - q**2*E1*E5*E3*E4*E3*E2 + q**3*E2*E1*E3*E4*E3*E5 - q**4*E2*E1*E3*E5*E4*E3 + q**3*E2*E1*E5*E3*E4*E3 - q**2*E3*E2*E1*E4*E3*E5 - q**2*E3*E2*E1*E5*E3*E4 + (q**3 + q)*E3*E2*E1*E5*E4*E3 + q*E3*E4*E3*E2*E1*E5 - 1*E3*E5*E4*E3*E2*E1 - q**2*E4*E3*E2*E1*E5*E3 - q**2*E5*E3*E2*E1*E4*E3 + q*E5*E3*E4*E3*E2*E1)



########## CENTRAL ELEMENT ############

Case0Pos = term1 + term2 + term3 + term4 + term5
Case0Neg = termNeg1 + termNeg2 + termNeg3 + termNeg4 + termNeg5
Case1 = term5121 + term5131 + term5141 + term5151 + term5231 + term5241 + term5251 + term5341 + term5351 + term5451
Case2 = term5112 + term5122 + term5132 + term5142 + term5152 + term5212 + term5222 + term5232 + term5242 + term5252 + term5312 + term5322 + term5332 + term5342 + term5352 + term5412 + term5422 + term5432 + term5442 + term5452 + term5512 + term5522 + term5532 + term5542
Case3 = term5213 + term5313 + term5413 + term5513 + term5323 + term5423 + term5523 + term5433 + term5533 + term5543

coproductRepresentation = Case0Pos + Case0Neg + Case1 + Case2 + Case3 # 100 x 100 matrix H


########## REARRANGING BASIS (to get a block diagonal matrix) ##########
indexList = [6, 17, 28, 39, 50, 51, 62, 73, 84, 95] # 10x10 block
for i in range(len(indexList)):
    coproductRepresentation = coproductRepresentation.elementary_row_op(op = "n<->m", row1 = i, row2 = indexList[i] - 1)
    coproductRepresentation = coproductRepresentation.elementary_col_op(op = "n<->m", col1 = i, col2 = indexList[i] - 1) 

fromIndex = [17, 22, 42, 72, 82, 28, 33, 43, 63, 83, 39, 36, 54, 44, 74, 94, 50, 53, 55, 65, 85, 85, 57, 66, 58, 76, 59, 86, 96, 94, 75, 77, 87, 73, 94, 79, 88, 82, 92, 98, 84, 93, 90, 99, 97, 91, 95] # 2x2 blocks
toIndex = [12, 14, 16, 17, 20, 22, 23, 26, 28, 30, 32, 33, 34, 36, 37, 39, 42, 44, 45, 48, 50, 51, 53, 54, 55, 56, 57, 58, 59, 64, 66, 67, 70, 72, 73, 75, 76, 77, 78, 79, 82, 84, 85, 86, 88, 89, 90]
for i in range(len(fromIndex)):
    coproductRepresentation = coproductRepresentation.elementary_row_op(op = "n<->m", row1 = toIndex[i] - 1, row2 = fromIndex[i] - 1)
    coproductRepresentation = coproductRepresentation.elementary_col_op(op = "n<->m", col1 = toIndex[i] - 1, col2 = fromIndex[i] - 1)

# 10 x 10 block, with respect to basis (v1 ⊗ v6, v2 ⊗ v7, v3 ⊗ v8, v4 ⊗ v9, v5 ⊗ v10, v6 ⊗ v1, v7 ⊗ v2, v8 ⊗ v3, v9 ⊗ v4, v10 ⊗ v5)
tenByTen = zeros(10) 
for i in range(10):
    for j in range(10):
        tenByTen[i,j] = simplify(simplify(coproductRepresentation[i,j]).subs(r, q - 1/q))

# 2 x 2 block
twoByTwo = zeros(2) 
for i in range(2):
    for j in range(2):
        twoByTwo[i,j] = simplify(simplify(coproductRepresentation[10 + i, 10 + j]).subs(r, q - 1/q))

# Lambda
Λ = simplify(coproductRepresentation[99, 99].subs(r, q - 1/q))

# Quantum Hamiltonian (of the 10 x 10 part)
Hhat = tenByTen - Λ * eye(10)

########### Conjugating the Hamiltonian to get L_δ ##############
z = Symbol('z') # zero
eVects = [[-q**2, q, z, z, z, -1, q, z, z, z],
          [z, -q**2, q, z, z, z, -1, q, z, z],
          [z, z, -q**2, q, z, z, z, -1, q, z],
          [z, z, z, -q**2, q, z, z, z, -1, q]]
conjs = [zeros(10) for i in range(len(eVects))]
for k in range(len(conjs)):
    for i in range(10):
        for j in range(10):
            conjs[k][i, j] = simplify((Hhat[i, j] * eVects[k][j] / eVects[k][i]).subs(z, 0))

########### Discarding states with negative probability to get \tilde{L}_δ  ###########
conjSubmatrices = [zeros(4) for i in range(len(conjs))]
# First conjugate matrix: remove rows/cols 3, 4, 5, 8, 9, 10 (1-indexed)
# Basis (v1 ⊗ v6, v2 ⊗ v7, v6 ⊗ v1, v7 ⊗ v2)
# Type D ASEP (5, 3)
conjSubmatrices[0] = conjs[0][(0, 1, 5, 6), (0, 1, 5, 6)]

# Second conjugate matrix: remove rows/cols 1, 4, 5, 6, 9, 10 (1-indexed)
# Basis (v2 ⊗ v7, v3 ⊗ v8, v7 ⊗ v2, v8 ⊗ v3)
# Type D ASEP (5, 2)
conjSubmatrices[1] = conjs[1][(1, 2, 6, 7), (1, 2, 6, 7)]
 
# Third conjugate matrix: remove rows/cols 1, 2, 5, 6, 7, 10 (1-indexed)
# Basis (v3 ⊗ v8, v4 ⊗ v9, v8 ⊗ v3, v9 ⊗ v4)
# Type D ASEP (5, 1)
conjSubmatrices[2] = conjs[2][(2, 3, 7, 8), (2, 3, 7, 8)]

# Fourth conjugate matrix: remove rows/cols 1, 2, 3, 6, 7, 8 (1-indexed)
# Basis (v4 ⊗ v9, v5 ⊗ v10, v9 ⊗ v4, v10 ⊗ v5)
# Type D ASEP (5, 0)
conjSubmatrices[3] = conjs[3][(3, 4, 8, 9), (3, 4, 8, 9)]

############ Type D ASEP ###############
# To verify that the conjugated submatrices are indeed the generator matrices multiplied by r^2

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
# Generator matrix with two sites (just the 4x4 block)
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
    gen[2, 1] = rTwistLeft(q, n, δ) * 1 / (1 + q**(2*δ))
    gen[2, 3] = rTwistLeft(q, n, δ) * q**(2*δ) / (1 + q**(2*δ)) 
    gen[2, 2] = - (gen[2, 0] + gen[2, 1] + gen[2, 3])
    gen[3, 0] = rStick(q, n)
    gen[3, 1] = rSwap(q, n)
    gen[3, 2] = q**(2*δ) * rStick(1/q, n)
    gen[3, 3] = - (gen[3, 0] + gen[3, 1] + gen[3, 2])   
    return gen

gens = [generator(q, 5, δ) for δ in range(4)]

