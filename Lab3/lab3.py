import numpy as np
import math

left = 1
right = 2

def function(x):
    return (-1) * (9 + x * x) * (9 + x * x)

def get_discrete_function(function, left, right, n):
    discrete_function = []
    for i in range(1, n):
        discrete_function.append(function(left + i * h))
    return discrete_function

def get_T(n):
    T = np.zeros((n - 1, n - 1))

    T[0][0] = -2
    T[0][1] = 1

    for row in range(1, n - 2):
        T[row][row - 1] = T[row][row + 1] = 1
        T[row][row] = -2

    T[n - 2][n - 3] = 1
    T[n - 2][n - 2] = -2

    return T / h / h

def get_A(n):
    return np.multiply(get_T(n), np.asarray(get_discrete_function(function, left, right, n))).transpose()

def get_B(n):
    return get_T(n) * -1

def get_fk(k, n):
    if k == 0:
        return np.ones(n - 1)
    else:
        fk = []
        for i in range(1, n):
            fk.append(math.pow(left + i * h, k))
        return fk

def get_gk(C, fk):
    return C.dot(fk)

def get_Cm(C):
    fk_arr = []
    for k in range(0, m):
        fk = get_fk(k, n)
        fk_arr.append(fk)

    gk_arr = []
    for k in range(0, m):
        gk = get_gk(C, fk_arr[k])
        gk_arr.append(gk)

    bij = np.zeros((m, m))
    for i in range(m):
        for j in range(m):
            bij[i][j] += np.sum(fk_arr[j] * gk_arr[i])
    bij = np.linalg.inv(bij)

    Cm = np.zeros((n - 1, n - 1))
    for i in range(m):
        gi = gk_arr[i].copy()
        gi.shape = (n - 1, 1)
        for j in range(m):
            Cm += (bij[i][j] * gi) * gk_arr[j]

    return Cm

def Weinstein(B, C, n):
    Cm = get_Cm(C)
    Am = B + Cm

    eigvals = np.linalg.eigvals(Am)
    return min([n for n in eigvals if n > 0])

n = 100
h = (right - left) / (n - 1)
m = 10
m_max = 11

A = get_A(n)
B = get_B(n)
C = A - B

print('m = ' + str(m))
lambd = Weinstein(B, C, n)
print('lambda = ' + str(lambd/100))