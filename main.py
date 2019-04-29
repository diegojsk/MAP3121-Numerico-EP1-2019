import math
from cmath import sqrt

#define parâmetros c e s 

def c (a,b):
    if a>b:
        T = -b/a
        coss = 1/sqrt(1+(T^2))
    else :
        T = -a/b
        coss = s(a,b) *T
    return coss

def s (a,b):
    if a>b:
        T = -b/a
        sen = c(a,b)*T
    else :
        T = -a/b
        sen = 1/sqrt(1+(T^2))
    return sen

#implementa Rotação de Givens para matriz W

def RotGivens(W,n,m,i,j,c,s):

    for r in range (1,n) :

        aux = c*W[i][r] - s*W[j][r]
        W[j][r] = s*W[i][r] + c*W[j][r]
        W[i][r] = aux

def FatoracaoQR (W,n,m,i,j):

    for k in range(1,m):
        for j in range(n,-1,k+1):
            i = j -1
            if W[j][k] != 0 :
                RotGivens(W,n,m,i,j,c(W[i][k],W[j][k]),s(W[i][k],W[j][k]))

