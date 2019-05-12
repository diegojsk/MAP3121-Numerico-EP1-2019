from main import resolver_sist
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

_H = resolver_sist(W, A)

print(_H)
# [[0.5 1.  0. ]
#  [0.5 0.  1. ]]
