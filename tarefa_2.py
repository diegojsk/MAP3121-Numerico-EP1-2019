"""
Segunda Tarefa
"""

from main import *
import numpy as np

A = np.array([[3/10, 3/5, 0],
               [1/2,   0, 1],
              [4/10, 4/5, 0]])

W = np.array([[3/5, 0],
                [0, 1],
              [4/5, 0]])

H = np.array([[1/2, 1, 0],
              [1/2, 0, 1]])
                
np.set_printoptions(precision=3, suppress=True)

<<<<<<< HEAD
X,Y = resolve_mmq(A, W)
print(X)
print(Y)
#print(np.matmul(X,Y))


=======
P, Q = resolve_mmq(A, W)

print(P)
print()
print(Q)
>>>>>>> b689866d9de4b835fa3139f44948db9f186a80ad
