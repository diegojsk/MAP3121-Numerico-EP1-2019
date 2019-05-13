import numpy as np
import math
from cmath import sqrt

ERR = 1e-4
MAX_ITER = 1e2
SUPPORTED_FORMATS = [np.float32, np.float64, np.float_,
                     np.complex64, np.complex128, np.complex_]


def calc_c(a, b):
    """
        Calcula o parâmetro c, definido na página 3 do enunciado
        :param a: w[i,k] - elemento da matriz W na posição (i, k)
        :param b: w[j,k] - elemento da matriz W na posição (j, k)
        :return: Cosseno do ângulo de rotação
    """
    if a.dtype not in SUPPORTED_FORMATS:
        raise TypeError("Matriz de formato não suportado: {}! \
            \nNote que tipos inteiros levam a perdas severas \
            por arredondamento.".format(a.dtype))

    if b.dtype not in SUPPORTED_FORMATS:
        raise TypeError("Matriz de formato não suportado: {}! \
            \nNote que tipos inteiros levam a perdas severas \
            por arredondamento.".format(b.dtype))

    if abs(a) > abs(b):
        T = -np.divide(b, a)
        cos = 1/np.sqrt(1+(T**2))
    else:
        T = -np.divide(a, b)
        cos = calc_s(a, b)*T
    return cos


def calc_s(a, b):
    """
        Calcula o parâmetro s, definido na página 3 do enunciado 
        :param a: double w[i,k] - elemento da matriz W na posição (i, k)
        :param b: double w[j,k] - elemento da matriz W na posição (j, k)

        :return: Seno do ângulo de rotação
    """
    if a.dtype not in SUPPORTED_FORMATS:
        raise TypeError("Matriz de formato não suportado: {}! \
            \nNote que tipos inteiros levam a perdas severas \
            por arredondamento.".format(a.dtype))

    if b.dtype not in SUPPORTED_FORMATS:
        raise TypeError("Matriz de formato não suportado: {}! \
            \nNote que tipos inteiros levam a perdas severas \
            por arredondamento.".format(b.dtype))

    if abs(a) > abs(b):
        T = -np.divide(b, a)
        sen = calc_c(a, b)*T
    else:
        T = -np.divide(a, b)
        sen = 1/np.sqrt(1+(T**2))
    return sen


def rot_givens(W, n, m, i, j, c, s):
    """
        Efetua a rotação de Givens para matriz W
        :param W: ndarray
        :param n: int Número de linhas de W
        :param m: int Número de colunas de W
        :param i: int Linha a ser rotacionada
        :param j: int Linha a ser rotacionada
        :param c: double Cosseno do ângulo de rotação
        :param s: double Seno do ângulo de rotação

        :return: None
    """
    if W.dtype not in SUPPORTED_FORMATS:
        raise TypeError("Matriz de formato não suportado: {}! \
            \nNote que tipos inteiros levam a perdas severas \
            por arredondamento.".format(W.dtype))

    col = 0
    while (W[i][col] == 0 and W[j][col] == 0):
        col += 1

    for r in range(col, m):
        aux = c*W[i][r] - s*W[j][r]
        W[j][r] = s*W[i][r] + c*W[j][r]
        W[i][r] = aux


def zera_elemento(W, i, j, k):
    """
        Realiza uma rotação de Givens de modo a zerar o elemento (j,k)
        :param W: ndarray
        :param i: linha a ser rotacionada
        :param j: linha a ser zerada
        :param k: coluna a ser zerada

        :return: None
    """
    if W.dtype not in SUPPORTED_FORMATS:
        raise TypeError("Matriz de formato não suportado: {}! \
            \nNote que tipos inteiros levam a perdas severas \
            por arredondamento.".format(W.dtype))

    n, m = W.shape
    _s = calc_s(W[i, k], W[j, k])
    _c = calc_c(W[i, k], W[j, k])
    return rot_givens(W, n, m, i, j, _c, _s)


def fatorar_qr(W):
    """
        Encontra a matriz R de modo que Q*R = W, onde Q é o resultado de
    sucessivas matrizes de rotação de Givens e R é uma matriz triangular
    superior.
        Ou seja, a função transforma a matriz W em uma matriz triangular
    superior por meio de sucessivas rotações de Givens
        :param W: ndarray

        :return: None
    """
    if W.dtype not in SUPPORTED_FORMATS:
        raise TypeError("Matriz de formato não suportado: {}! \
            \nNote que tipos inteiros levam a perdas severas \
            por arredondamento.".format(W.dtype))

    n, m = W.shape

    for k in range(m):
        for j in range(n-1, k, -1):
            i = j-1
            if W[j][k] != 0:
                # zera_elemento(W,W,i,j,k)
                _s = calc_s(W[i][k], W[j][k])
                _c = calc_c(W[i][k], W[j][k])
                rot_givens(W, n, m, i, j, _c, _s)


def resolver_sist(W, A):
    """
        Dadas matrizes W e A, encontra a matriz H, tal que W*H = A

        Função Principal da Primeira Tarefa c) d)

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

    H = np.ones((p, m))

    k = 1
    while k <= p:
        j = n
        while j >= k+1:
            i = j-1
            if W[j-1][k-1] != 0:
                # n, m = W.shape
                _s = calc_s(W[i-1][k-1], W[j-1][k-1])
                _c = calc_c(W[i-1][k-1], W[j-1][k-1])
                rot_givens(W, n, p, i-1, j-1, _c, _s)
                rot_givens(A, n, m, i-1, j-1, _c, _s)
            j -= 1
        k += 1

    k = p
    while k >= 1:

        j = 1
        while j <= m:
            soma = 0.0
            i = k + 1
            while i <= p:
                soma += W[k-1][i-1]*H[i-1][j-1]
                i += 1
            H[k-1][j-1] = (A[k-1][j-1] - soma)/W[k-1][k-1]
            j += 1
        k -= 1

    return H


def residuo(A, W, H):
    """
        Calcula o quadrado da norma da matriz E = A-W*H

        :param W: ndarray n;m
        :param A: ndarray n;p
        :param W: ndarray p;m
    """
    na, ma = A.shape
    nw, pw = W.shape
    ph, mh = H.shape

    if na != nw or ma != mh or ph != pw:
        raise ValueError("Matrizes não compatíveis")

    WH = np.matmul(W, H)

    E = A - WH
    erro = np.sum(np.power(E, 2))
    return erro


def normaliza(M):
    """
    Normaliza a matriz M, de modo que a norma de todas as suas colunas
    é igual a zero

        :param M: ndarray
    """
    # soma_colunas = np.power(M.sum(axis=0), 1)
    # n, m = M.shape
    # np.divide(M, soma_colunas, out=M)

    n, m = M.shape
    M_n = np.multiply(M, M)
    soma_colunas = np.sum(M_n, axis=0)
    for i in range(n):
        for j in range(m):
            M[i][j] = np.divide(M[i][j], np.sqrt(soma_colunas[j]))


def calc_transpose(M):
    """
    Calcula a transposta da matriz M
        :param M:
    """
    n, m = M.shape
    M_t = np.empty((m, n))
    for i in range(n):
        for j in range(m):
            M_t[j][i] = M[i][j]
    return M_t


def resolve_mmq(A, W0):

    """
    Resolve o  MMQ
    Função Principal da Segunda Tarefa
        :param A: Matriz a ser fatorada
        :param W0: Matriz W0 para determinação dos fatores
    """

    n1, p = W0.shape
    n2, m = A.shape

    n = None
    if n1 != n2:
        raise ValueError("Matrizes com tamanhos não compatíveis!")
    else:
        n = n1

    H = np.ones((p, m))
    _H = np.ones((p, m))
    W = W0.copy()
    _A = A.copy()
    # W = np.random.rand(n, p)
    # W = np.ones((n, p))

    i = 0
    err = residuo(_A, W0, H)
    prev_err = err + 1

    while (np.power(err - prev_err, 2) > np.power(ERR, 2)) and (i < MAX_ITER):

        i += 1
        print(W)
        normaliza(W)
        print(W)

        H = resolver_sist(W, A)

        H[H < 0] = 0.0
        _H = H.copy()

        A = _A.copy()

        A_t = A.transpose()
        H_t = H.transpose()

        W_t = resolver_sist(H_t, A_t)

        W = W_t.transpose()
        print(W)

        W[W < 0] = 0.0

        A = _A.copy()
        H = _H.copy()

        prev_err = err
        err = residuo(A, W, H)


    return (W, H)


def matriz_arquivo(arquivo):
    '''
    Lê arquivo.txt e transforma em array Matriz
    '''
    arq = open(arquivo, "r+")
    texto = arq.readlines()
    matriz = []

    for linha in texto:
        linha = linha.strip('\n')
        linha = linha.split(' ')
        for i in range(len(linha)):
            linha[i] = float(linha[i])
        matriz.append(linha)
    return matriz

def treinamento(d):
    '''
    Executar treinamento do dígito d gerando a matriz Wd
    '''

    A = np.array(matriz_arquivo('dados_mnist/train_dig'+str(d)+'.txt'))
    n, m = A.shape
    p = 10
    W = np.random.rand(n, p)
    Wd, H = resolve_mmq(A, W)

    np.save(Wd, 'treinamento/W'+str(d)+'.npy')

    return Wd


if __name__ == "__main__":

    '''
    Matriz W do enunciado
    '''
    np.set_printoptions(precision=3, suppress=True)

    W = np.array([[2,  1,  1, -1,  1],
                  [0,  3,  0,  1,  2],
                  [0,  0,  2,  2, -1],
                  [0,  0, -1,  1,  2],
                  [0,  0,  0,  3,  1]]).astype(np.double)

    zera_elemento(W, 2, 3, 2)
    print(W*np.sqrt(5))
    print(W)

    '''
    Verificar se é triangular
    '''

    fatorar_qr(W)
    print(W)
    for i in range(5):
        for j in range(5):
            if i > j:
                print(W[i][j])

    '''
    Matriz b qualquer
    '''

    # b = np.array([[1],[1],[1],[1],[1]]).astype(np.double)
    # print(b)
    # zera_elemento(b,W,2,3,0)
    # print(b)



