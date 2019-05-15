"""
Segunda Tarefa
"""

from main import *
import numpy as np

'''
Teste Segunda Tarefa - exemplo do enunciado
'''

A = np.array([[3/10, 3/5, 0],
              [1/2,   0, 1],
              [4/10, 4/5, 0]])

_A = A.copy()

W = np.array([[3/5, 0],
              [0, 1],
              [4/5, 0]])

H = np.array([[1/2, 1, 0],
              [1/2, 0, 1]])

np.set_printoptions(precision=3, suppress=True)

# P, Q = resolve_mmq(_A, W)

# print(P)
# print()
# print(Q)

'''
Outros exemplos
'''

W = np.array([[96759699.3224, -2424353.0], [2245.78, 2535.0], [0.34434, 34343.434]])

_A = A.copy()

np.set_printoptions(precision=3, suppress=True)

P, Q = resolve_mmq(_A, W)

print(P)
print()
print(Q)

print(np.matmul(P, Q))


