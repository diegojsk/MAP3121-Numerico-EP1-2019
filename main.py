import numpy as np
import math
from cmath import sqrt

#define parâmetros c e s 

def c (a,b):
    """
    Calcula o parâmetro c, definido na página 3 do enunciado 
        :param a: w[i,j] - elemento da matriz W na posição (i, j)
        :param b: w[j,i] - elemento da matriz W na posição (j, i)
    """
    
    if a>b:
        T = -b/a
        cos = 1/sqrt(1+(T^2))
    else :
        T = -a/b
        cos = s(a,b) *T
    return cos

def s (a,b):
    """
    Calcula o parâmetro s, definido na página 3 do enunciado 
        :param a: w[i,j] - elemento da matriz W na posição (i, j)
        :param b: w[j,i] - elemento da matriz W na posição (j, i)
    """
    if a>b:
        T = -b/a
        sen = c(a,b)*T
    else :
        T = -a/b
        sen = 1/sqrt(1+(T^2))
    return sen

def RotGivens(W,i,j,c,s):
    """
    Implementa Rotação de Givens para matriz W
        :param W: ndarray
        :param i: linha a ser zerada
        :param j: coluna a ser zerada
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

def FatoracaoQR (W,i,j):
    """
    docstring here
        :param W: ndarray
        :param i: 
        :param j: 
    """

    n, m = W.shape
    for k in range(m):
        for j in range(n-1,-1,k):
            if W[j][k] != 0 :
                RotGivens(W,i,j,c(W[i][k],W[j][k]),s(W[i][k],W[j][k]))

