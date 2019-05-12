from main import fatorar_qr
import numpy as np

A = np.array([ [1, 2, 3, 4],
               [2, 3, 4, 5],
               [3, 4, 5, 6],
               [4, 5, 6, 7] ]).astype(np.double)

B = np.array([ [1, -2, 3, 4],
               [-2, -3, 4, 5],
               [3, 4, -5, -6],
               [1, -9, 7, -10]]).astype(np.double)

np.set_printoptions(precision=3, suppress=True)

print("{0:.1f}".format(np.linalg.det(A)) )
# 0.0
print("{0:.1f}".format(np.linalg.det(B)))
# 30.0

fatorar_qr(A)
fatorar_qr(B)

print("{0:.1f}".format(np.linalg.det(A)))
# 0.0
print("{0:.1f}".format(np.linalg.det(B)))
# 30.0

print(A)
# [[ -5.477  -7.303  -9.129 -10.954]
#  [ -0.      0.816   1.633   2.449]
#  [ -0.      0.      0.      0.   ]
#  [  0.      0.      0.      0.   ]]

print(B)
# [[  3.873   1.807  -3.357  -8.779]
#  [  0.    -10.331   9.189  -5.698]
#  [ -0.     -0.      1.815   8.203]
#  [  0.      0.      0.     -0.413]]