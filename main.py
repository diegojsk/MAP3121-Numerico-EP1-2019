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
        T = -np.divide(b,a)
        cos = 1/np.sqrt(1+(T**2))
    else :
        T = -np.divide(a,b)
        cos = calc_s(a,b)*T
    return cos

def calc_s (a,b):
    """
    Calcula o parâmetro s, definido na página 3 do enunciado 
        :param a: w[i,k] - elemento da matriz W na posição (i, k)
        :param b: w[j,k] - elemento da matriz W na posição (j, k)
    """
    if abs(a) > abs(b):
        T = -np.divide(b,a)
        sen = calc_c(a,b)*T
    else :
        T = -np.divide(a,b) 
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

def zera_elemento(W,Wc,i,j,k):
    """
    Realiza uma rotação de Givens de modo a zerar o elemento (j,k)
        :param W: ndarray
        :param Wc: array que define seno e cosseno
        :param i: linha a ser rotacionada
        :param j: linha a ser zerada
        :param k: coluna a ser zerada
    """
    n, m = W.shape
    _s = calc_s(Wc[i,k], Wc[j,k])
    _c = calc_c(Wc[i,k], Wc[j,k])
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
               zera_elemento(W,W,i,j,k)


def resolver_sist(W,A):
    """
    Resolve sistema simultâneos para W e A
        :param W: ndarray n;p
        :param A: ndarray n;m
    """
    n1, p = W.shape
    n2, m = A.shape
    n = None
    if n1 != n2:
        raise ValueError("Matrizes de tamanhos incompatíveis!")
    else:
        n = n1

    H = np.zeros((p, m))

    for k in range(p):
        for j in range(n-1,k,-1):
            i = j-1
            if W[j][k] != 0 :
                # n, m = W.shape
                _s = calc_s(W[i,k], W[j,k])
                _c = calc_c(W[i,k], W[j,k])
                rot_givens(W,n,p,i,j,_c,_s)
                rot_givens(A,n,m,i,j,_c,_s)

    for k in range(p-1, 0,-1):
        soma = 0
        for i in range(k,p-1):
            soma = soma + W[k][i]*H[i][j]
        for j in range(m):
            H[k][j] = (A[k][j] - soma)/W[k][k]

    return H

    
    def residuo(A,W,H):
        """
        Calculo residuo para (A-WH)
            :param W: ndarray n;m
            :param A: ndarray n;p
            :param W: ndarray p;m
        """
        na,ma = A.shape
        nw,pw = W.shape
        ph,mh = H.shape

        if na != nw or ma != mh or ph != pw :
            raise ValueError("Matrizes não compatíveis")
        else:
            n = na
            m = mh
            p = ph

        WH = W*H
        erro = 0
        for i in range(n):
            for j in range(m):
                erro = erro + (A[i][j] - WH[i][j])**2

        return erro

    
def normaliza(M):
    """
    Normaliza a matriz M
        :param M: 
    """
    soma_colunas = np.sum(M, axis=0)
    print(soma_colunas)
    n, m = M.shape
    for i in range(n):
        for j in range(m):
            M[i][j] = np.divide(M[i][j], soma_colunas[j])

def calc_transpose(M):
    """
    Calcula a transposta da matriz M
        :param M: 
    """
    n, m = M.shape
    M_t = np.array(m, n)
        for i in range(n):
            for j in range(m):
                M_t[j][i] = M[i][j]
    return M_t

def resolve_mmq(A, H, p, err):

    n, m = A.shape

    _A = A.copy()
    W = np.empty((n, p))

    while residuo(A, W, H) < err:

        H = resolver_sist(W, A)
        H[ H < 0 ] = 0

        A = _A.copy()
        A_t = calc_transpose(A)
        H_t = calc_transpose(H)

        W_t = resolver_sist(H_t, A_t)

        W = calc_transpose(W_t)
        W[ W < 0 ] = 0

    return W

if __name__ == "__main__":

    '''
    Matriz W do enunciado
    '''

    W = np.array([[ 2,  1,  1, -1,  1],
                  [ 0,  3,  0,  1,  2],
                  [ 0,  0,  2,  2, -1],
                  [ 0,  0, -1,  1,  2],
                  [ 0,  0,  0,  3,  1.0]])

    zera_elemento(W,W, 2, 3, 2)
    print(W*np.sqrt(5))
    print(W)
    
    '''
    Matriz b qualquer
    '''
    
    b = np.array([[1],[1],[1],[1],[1]]).astype(np.double)
    print(b)
    zera_elemento(b,W,2,3,0)
    print(b)
    
    """
    Primeira Tarefa
    """
    
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
    b = np.ones((n,1))
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
    B = np.zeros((n,m))
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
    b = np.zeros((n,1))
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
    
    H = resolver_sist(W,A)
    print(H)

    """
    item d)
    """

    n = 20
    p = 17
    W = np.zeros((n,m))
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
    print(H)

    """
    Segunda Tarefa
    """

    A = np.array([[3/10,3/5,0],
                  [1/2,0,1],
                  [4/10,4/5,0]])
    
    W = np.array([[3/5,0],
                  [0,1],
                  [4/5,0]])

    H = W = np.array([[1/2,1,0],
                      [1/2,0,1]])
    
    

    


    