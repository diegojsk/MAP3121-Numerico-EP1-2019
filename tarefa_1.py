"""
Primeira Tarefa
"""

from main import *
import numpy as np

np.set_printoptions(precision=3, suppress=True)

"""
item a)
"""

n = 64
m = 64
A = np.zeros((n, m))
for i in range(n):
    for j in range(m):
        if i == j:
            A[i][j] = 2
        elif abs(i-j) == 1:
            A[i][j] = 1
        elif abs(i-j) > 1:
            A[i][j] = 0
        else:
            A[i][j] = 0

b = np.ones((n, 1))

print(resolver_sist(A, b))

"""
item b)
"""

n = 20
m = 17
B = np.zeros((n, m))

for i in range(n):
    for j in range(m):
        if abs(i-j) <= 4:
            B[i][j] = 1/((i+1+j+1-1))
        elif abs(i-j) > 4:
            B[i][j] = 0
        else:
            B[i][j] = 0

b = np.zeros((n, 1))
for i in range(n):
    b[i] = i + 1

print(resolver_sist(B, b))
