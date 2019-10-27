import numpy as np
import math

n = 10
a = 0
b = 2
h = (b - a) / ( n - 1)
l = 2.04756

def f(x):
    return (-1) * (9 + x * x) * (9 + x * x)

def build_2d_mat(size: int) -> np.ndarray:
    mat = np.zeros([size, size], np.float32)
    for i in range(size):
        mat[i, i] = -2
        if i+1 < size:
            mat[i, i+1] = 1
            mat[i+1, i] = 1
    return mat / h / h

def integrate(func: np.ndarray, h: float) -> float:
    areas = []
    for fi in func:
        areas.append(fi * h)
    return np.sum(areas)

space = np.linspace(a, b, n)
f0 = np.vectorize(f)(space)
T = build_2d_mat(f0.size)
fk_1 = f0.copy()
a_list = [integrate(fk_1 * fk_1, (b - a) / n)]
for i in range(20):
    fk = -np.dot(np.linalg.inv(T), fk_1)
    ak = integrate(fk_1 * fk, (b - a) / n)
    a_list.append(ak)
    fk_1 = fk

print("Kolats:")
mk = []
for i in range(len(a_list)-1):
    mk.append(a_list[i]/a_list[i+1])

vk = []
l2 = 4
for i in range(len(mk) - 1):
    vk_ = mk[i + 1] - (mk[i] - mk[i + 1]) / (l2 / mk[i + 1] - 1)
    vk.append(vk_)
for i in range(len(mk) - 1):
    print(str((mk[i]+vk[i])/2*l))

kbl = []
kbu = []
for i in range(len(mk) - 1):
    kbl.append(mk[i + 1] - np.sqrt(np.abs((mk[i] - mk[i + 1]) * mk[i + 1])))
    kbu.append(mk[i + 1] + np.sqrt(np.abs((mk[i] - mk[i + 1]) * mk[i + 1])))

print("\nKrolov-Bogolubov:")
for i in range(len(kbl)-1):
    print(str(kbl[i]*l)+"     " + str((mk[i]+vk[i])/2*l) + "     " + str(kbu[i]*l))