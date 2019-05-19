"""
Segunda Tarefa
"""

from main import *
import numpy as np

np.set_printoptions(precision=3, suppress=True)

'''
Teste Segunda Tarefa - exemplo do enunciado
'''

A = np.array([[3/10, 3/5, 0],
              [1/2,   0, 1],
              [4/10, 4/5, 0]])

_A = A.copy()

W = np.array([[3/5, 0],
              [0, 1/2],
              [4/5, 0]])

H = np.array([[1/2, 1, 0],
              [1/2, 0, 1]])

np.set_printoptions(precision=3, suppress=True)

P, Q = resolve_mmq(_A, W)

print(P)
print()
print(Q)

'''
item c)
'''

n = 64
p = 64
W = np.zeros((n, p))
for i in range(n):
    for j in range(p):
        if i == j:
            W[i][j] = 2
        elif abs(i-j) == 1:
            W[i][j] = 1
        elif abs(i-j) > 1:
            W[i][j] = 0
        else:
            W[i][j] = 0

m = 3
A = np.zeros((n, m))
for i in range(n):
    for j in range(m):
        if j == 1-1 :
            A[i][j] = 1
        elif j == 2-1:
            A[i][j] = i + 1
        elif j == 3-1:
            A[i][j] = 2*(i+1) - 1

H = resolver_sist(W, A)
print(H)


"""
item d)
"""

n = 20
p = 17
W = np.zeros((n, p))
##print(W)
for i in range(n):
    for j in range(p):
        if abs(i-j) <= 4:
            W[i][j] = 1/((i+1+j+1-1))
        elif abs(i-j) > 4:
            W[i][j] = 0
        else:
            W[i][j] = 0

m = 3
A = np.zeros((n, m))
for i in range(n):
    for j in range(m):
        if j == 1-1:
            A[i][j] = 1
        elif j == 2-1:
            A[i][j] = i + 1
        elif j == 3-1:
            A[i][j] = 2*(i+1) - 1

H = resolver_sist(W, A)

print(H)
