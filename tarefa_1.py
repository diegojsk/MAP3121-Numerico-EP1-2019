"""
Primeira Tarefa
"""

from main import *
import numpy as np

"""
item a) 
"""
    
n = 64
m = 64
A = np.zeros((n,m))
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
                
print(A)
fatorar_qr(A)


for i in range(n):
    for j in range(m):
        if i > j :
            print(A[i][j] , i, j)


b = np.ones((n,1))
print(b)

print(resolver_sist(A,b))


"""
item b)
"""


n = 20
m = 17
B = np.zeros((n,m))
for i in range(n):
    for j in range(m):
        if abs(i-j) <= 4:
            B[i][j] = 1/((i+1+j+1-1))
        elif abs(i-j) > 4:
            B[i][j] = 0
        else:
            B[i][j] = 0

#print(B)
fatorar_qr(B)
for i in range(n):
    for j in range(j):
        if i > j :
            print(B[i][j],i,j)

b = np.zeros((n,1))
for i in range(n):
    b[i] = i + 1
#print(b)

#print(resolver_sist(B,b))

'''
item c)
'''

n = 64 
p = 64
W = np.zeros((n,p))
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
m=3
A = np.zeros((n,m))
for i in range(n):
    for j in range(m):
        if j == 1-1 :
            A[i][j] = 1
        elif j == 2-1:
            A[i][j] = i + 1
        elif j == 3-1:
            A[i][j] = 2*(i+1) - 1

#print(W)
#fatorar_qr(W)
#print(W)
#print(A)
#H = resolver_sist(W,A)
#print(H)



"""
item d)
"""


n = 20
p = 17
W = np.zeros((n,m))
#print(W)
for i in range(n):
    for j in range(m):
        if abs(i-j) <= 4:
            W[i][j] = 1/((i+1+j+1-1))
        elif abs(i-j) > 4:
            W[i][j] = 0
        else:
            W[i][j] = 0

m=3
A = np.zeros((n,m))
for i in range(n):
    for j in range(m):
        if j == 1-1 :
            A[i][j] = 1
        elif j == 2-1:
            A[i][j] = i + 1
        elif j == 3-1:
            A[i][j] = 2*(i+1) - 1

H = resolver_sist(W,A)
#print(H)
