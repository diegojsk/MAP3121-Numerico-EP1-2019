import numpy as np
import math
from cmath import sqrt

#define parâmetros c e s 

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

    for r in range(col, n):
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

def fatorar_qr (W,i,j):
    """
    docstring here
        :param W: ndarray
        :param i: 
        :param j: 
    """
    n, m = W.shape
    for k in range(m):
        for j in range(n-1,-1,k):
            i = j-1
            if W[j][k] != 0 :
                zera_elemento(W,i,j,k)


if __name__ == "__main__":

    W = np.array([[ 2,  1,  1, -1,  1],
                  [ 0,  3,  0,  1,  2],
                  [ 0,  0,  2,  2, -1],
                  [ 0,  0, -1,  1,  2],
                  [ 0,  0,  0,  3,  1.0]])

    zera_elemento(W, 2, 3, 2)

    print(W*np.sqrt(5))
    print(W)
