import numpy as np
import math
from cmath import sqrt

def calc_c (a,b):
    """
    Calcula o parâmetro c, definido na página 3 do enunciado 
        :param a: w[i,k] - elemento da matriz W na posição (i, k)
        :param b: w[j,k] - elemento da matriz W na posição (j, k)
    """
    if abs(a) > abs(b):
        T = -b/a
        cos = 1/np.sqrt(1+(T**2))
    else :
        T = -a/b
        cos = calc_s(a,b)*T
    return cos

def calc_s (a,b):
    """
    Calcula o parâmetro s, definido na página 3 do enunciado 
        :param a: w[i,k] - elemento da matriz W na posição (i, k)
        :param b: w[j,k] - elemento da matriz W na posição (j, k)
    """
    if abs(a) > abs(b):
        T = -b/a
        sen = calc_c(a,b)*T
    else :
        T = -a/b
        sen = 1/np.sqrt(1+(T**2))
    return sen

def rot_givens(W,n,m,i,j,c,s):
    """
    Implementa Rotação de Givens para matriz W
        :param W: ndarray
        :param i: linha a ser rotacionada
        :param j: linha a ser zerada
        :param c: 
        :param s: 
    """
    col = 0
    while (W[i][col] == 0 and W[j][col] == 0):
        col += 1

    for r in range(col,m):
        aux = c*W[i][r] - s*W[j][r]
        W[j][r] = s*W[i][r] + c*W[j][r]
        W[i][r] = aux

def zera_elemento(W,i,j,k):
    """
    Realiza uma rotação de Givens de modo a zerar o elemento (j,k)
        :param W: ndarray
        :param i: linha a ser rotacionada
        :param j: linha a ser zerada
        :param k: coluna a ser zerada
    """
    n, m = W.shape
    _s = calc_s(W[i,k], W[j,k])
    _c = calc_c(W[i,k], W[j,k])
    return rot_givens(W, n, m, i, j, _c, _s)

def fatorar_qr (W):
    """
    docstring here
        :param W: ndarray
    """
    n, m = W.shape
    for k in range(m):
        for j in range(n-1,k,-1):
            i = j-1
            if W[j][k] != 0 :
               zera_elemento(W,i,j,k)


def resolver_sist(W,A):
    """
    Resolve sistema simultâneos para W e A
    """
    n,p = W.shape
    n,m = A.shape
    H = np.array([p*[m*[0]]])

    for k in range(1,p):
        for j in range(n-1,k,-1):
            i = j-1
            if W[j][k] != 0 :
                zera_elemento(W,i,j,k)
                zera_elemento(A,i,j,k)

    for k in range(p,1,-1):
        soma = 0
        for i in range(k+1,p):
            soma = soma + W[k][i]*H[i][j]
        for j in range(1,m):
            H[j][k] = (A[k][j] - soma)/W[k][k]

    return H
    
if __name__ == "__main__":

    '''
    Matriz W do enunciado
    '''

    W = np.array([[ 2,  1,  1, -1,  1],
                  [ 0,  3,  0,  1,  2],
                  [ 0,  0,  2,  2, -1],
                  [ 0,  0, -1,  1,  2],
                  [ 0,  0,  0,  3,  1.0]])

    zera_elemento(W, 2, 3, 2)

    print(W*np.sqrt(5))
    print(W)
    
    '''
    Matriz b qualquer
    '''
    
    b = np.array([[1],[1],[1],[1],[1]])
    print(b)
    zera_elemento(b,2,3,0)
    print(b)
    
    
    """
    item a) 
    """
    
    n = 64
    m=64
    A = np.array(n*[m*[0]])
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
    b = np.array(n*[[1]])
    print(b)
    
    fatorar_qr(A)
    fatorar_qr(b)
    print(A)
    print(b)

    """
    item b)
    """
    
    n = 20
    m = 17
    B = np.array(n*[m*[0]])
    print(B)
    for i in range(n):
        for j in range(m):
            if abs(i-j) <= 4:
                B[i][j] = 1/((i+1+j+1-1))
            elif abs(i-j) > 4:
                B[i][j] = 0
            else:
                B[i][j] = 0
    print(B)
    b = np.array(n*[[0]])
    for i in range(n):
        b[i] = i + 1
    print(b)
    
    fatorar_qr(B)
    fatorar_qr(b)
    print(B)
    print(b)


    '''
    item c)
    '''


    n = 64 
    p=64
    W = np.array(n*[p*[0]])
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
    for i in range(n):
        for j in range(m):
            if j == 1-1 :
                A[i][j] = 1
            elif j == 2-1:
                A[i][j] = i + 1
            elif j == 3-1:
                A[i][j] = 2*(i+1) - 1
    
    H = resolver_sist(W,A)
    print(H)

    """
    item d)
    """

    n = 20
    p = 17
    W = np.array(n*[m*[0]])
    print(B)
    for i in range(n):
        for j in range(m):
            if abs(i-j) <= 4:
                W[i][j] = 1/((i+1+j+1-1))
            elif abs(i-j) > 4:
                W[i][j] = 0
            else:
                W[i][j] = 0

    m=3
    for i in range(n):
        for j in range(m):
            if j == 1-1 :
                A[i][j] = 1
            elif j == 2-1:
                A[i][j] = i + 1
            elif j == 3-1:
                A[i][j] = 2*(i+1) - 1
    
    H = resolver_sist(W,A)
    print(H)

    


    