"""
Segunda Tarefa
"""

from main import *
import numpy as np

A = np.array([[3/10, 3/5, 0],
               [1/2,   0, 1],
              [4/10, 4/5, 0]])

W = np.array([[3/5,0],
                [0,1],
                [4/5,0]])

H = np.array([[1/2,1,0],
                [1/2,0,1]])
                
np.set_printoptions(precision=3, suppress=True)

print(resolve_mmq(A, W, 1e-5))