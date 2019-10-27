import numpy
import math

left = 1
right = 2
deviation = 0.01
n = 200
n_max = 3000

def f(x):
    return (-1) * (9 + x * x) * (9 + x * x)

def discrete(f, n):
    function = []
    for index in range(0, n):
        function.append(f(left + index * h))
    return function



def T_matrix(n):
    T = numpy.zeros((n, n))
    T[0][0] = -2
    T[0][1] = 1

    for row in range(1, n - 1):
        T[row][row - 1] = T[row][row + 1] = 1
        T[row][row] = -2

    T[n - 1][n - 2] = 1
    T[n - 1][n - 1] = -2

    return T / h / h

def interpret(l):
    return (l/2-41)/100

def A_matrix(n):
    matrix = numpy.multiply(T_matrix(n), discrete(f, n))

    return matrix.transpose()

lambda_previous = float("inf")
while n <= n_max:
    h = (right - left) / (n - 1)
    print('n: ' + str(n))
    A = A_matrix(n)
    inverse = numpy.linalg.inv(A)
    A_eigen_values = numpy.linalg.eigvalsh(inverse)
    m = max([i for i in A_eigen_values if i > 0])
    lambda1 = 1 / m
    print('lambda = ' + str(interpret(lambda1)))
    if abs(lambda1 - lambda_previous) < deviation:
        break;
    lambda_previous = lambda1
    n += 100
    print('-' * 100)