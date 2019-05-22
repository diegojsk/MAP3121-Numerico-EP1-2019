'''
Diego Jun Sato Kurashima - 10274231
Felipe Gomes de Melo D'Elia - 10340624
'''
import numpy as np
import math
import time
from cmath import sqrt
import os

ERR = 1e-4
MAX_ITER = 1e2
SUPPORTED_FORMATS = [np.float32, np.float64, np.float_,
                     np.complex64, np.complex128, np.complex_]


def calc_c(a, b):
    """
    Calcula o parâmetro c, correspondente ao cosseno do ângulo da rotação
    de Givens

        :param a: w[i,k] - elemento da matriz W na posição (i, k)
        :param b: w[j,k] - elemento da matriz W na posição (j, k)

        :return: Cosseno do ângulo de rotação
    """

    if abs(a) > abs(b):
        T = -np.divide(b, a)
        return np.divide(1, np.sqrt(1 + np.power(T, 2)))
    else:
        T = -np.divide(a, b)
        return np.divide(T, np.sqrt(1 + np.power(T, 2)))


def calc_s(a, b):
    """
    Calcula o parâmetro s, correspondente ao seno do ângulo da rotação
    de Givens

        :param a: double w[i,k] - elemento da matriz W na posição (i, k)
        :param b: double w[j,k] - elemento da matriz W na posição (j, k)

        :return: Seno do ângulo de rotação
    """

    if abs(a) > abs(b):
        T = -np.divide(b, a)
        return np.divide(T, np.sqrt(1 + np.power(T, 2)))
    else:
        T = -np.divide(a, b)
        return np.divide(1, np.sqrt(1 + np.power(T, 2)))


def rot_givens(W, n, m, i, j, c, s):
    """
    Efetua a rotação de Givens para matriz W
    ! Atenção, as alterações são feitas in place na matriz W

        :param W: ndarray
        :param n: int Número de linhas de W
        :param m: int Número de colunas de W
        :param i: int Linha a ser rotacionada
        :param j: int Linha a ser rotacionada
        :param c: double Cosseno do ângulo de rotação
        :param s: double Seno do ângulo de rotação

        :return: None
    """

    '''
    Implementação mais eficiente da Rotação de Givens
    '''
    aux = c*W[i, :] - s*W[j, :]
    W[j, :] = s*W[i, :] + c*W[j, :]
    W[i, :] = aux

    '''
    Pseudo-código do enunciado
    '''

    # for r in range(col, m):
    # n, m = W.shape
    # col = 0
    # while (col < m) and (W[i][col] == 0 and W[j][col] == 0):
    #     col += 1

    #     aux = c*W[i][r] - s*W[j][r]
    #     W[j][r] = s*W[i][r] + c*W[j][r]
    #     W[i][r] = aux


def zerar_elemento(W, i, j, k):
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
    superior por meio de sucessivas rotações de Givens.

        :param W: ndarray shape(n, m)

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
                # zerar_elemento(W,W,i,j,k)
                _s = calc_s(W[i][k], W[j][k])
                _c = calc_c(W[i][k], W[j][k])
                rot_givens(W, n, m, i, j, _c, _s)


def resolver_sist(W, A):
    """
    Dadas matrizes W e A, encontra a matriz H, tal que W*H = A

    Função Principal da Primeira Tarefa

        :param W: ndarray shape(n, m)
        :param A: ndarray shape(n, p)
    """
    n1, m = W.shape
    n2, p = A.shape
    n = n1

    if n1 != n2:
        raise ValueError("Matrizes de tamanhos incompatíveis!")

    H = np.ones((m, p))
    k = 1
    while k <= m:
        print("Columns: {0:02}/{1:02}".format(k, m), end='\r')
        j = n
        while j >= k+1:
            print("Columns: {2:02}/{3:02} Lines: {0:03}/{1:03}"
                  .format(n - j + 1, n - k + 1, k, m), end='\r')
            i = j-1
            if W[j-1][k-1] != 0:
                _s = calc_s(W[i-1][k-1], W[j-1][k-1])
                _c = calc_c(W[i-1][k-1], W[j-1][k-1])
                rot_givens(W, n, m, i-1, j-1, _c, _s)
                rot_givens(A, n, m, i-1, j-1, _c, _s)
            j -= 1
        k += 1
    print()

    k = m
    while k >= 1:
        soma = np.zeros(p)
        j = k + 1
        while j <= m:
            soma += W[k-1][j-1]*H[j-1, :]
            j += 1
        H[k-1, :] = (A[k-1, :] - soma)/W[k-1][k-1]
        k -= 1

    return H


def residuo(A, W, H):
    """
    Calcula o resíduo quadrático para | A - W*H |

        :param A: ndarray shape(n, p)
        :param W: ndarray shape(p, m)
        :param W: ndarray shape(n, m)

        :return erro: Erro quadrático para a fatoração
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


def normalizar(M):
    """
    Normaliza a matriz M, de modo que a norma de todas as suas colunas
    seja igual a um

        :param M: ndarray

        :return: None
    """

    n, m = M.shape
    M_n = np.multiply(M, M)
    soma_colunas = np.sum(M_n, axis=0)
    np.divide(M, np.sqrt(soma_colunas), out=M)


def fatorar_wh(A, p):

    """
    Encontra uma fatoração para a matriz A na forma A = W*H de modo que
    W e H são matrizes não-negativas

    Função Principal da Segunda Tarefa

        :param A: ndarray Matriz a ser fatorada
        :param p: int Quantidade de colunas de W

        :return W: ndarray Matriz W da fatoração não negativa de A
        :return H: ndarray Matriz H da fatoração não negativa de A
    """

    n, m = A.shape

    H = np.ones((p, m))
    _H = np.ones((p, m))
    _A = A.copy()
    W = np.random.rand(n, p)

    print("Starting MMQ algorithm")

    i = 0
    err = residuo(_A, W, H)
    prev_err = 0.0

    while (np.abs(err - prev_err)/ERR > err) and (i < MAX_ITER):
        begin = time.time()
        i += 1
        normalizar(W)
        print("- Solving for H")
        H = resolver_sist(W, A)

        H[H < 0.0] = 0.0
        _H = H.copy()

        A = _A.copy()

        A_t = A.transpose()
        H_t = H.transpose()

        print("- Solving for W_t")
        W = resolver_sist(H_t, A_t).transpose()

        W[W < 0.0] = 0.0

        A = _A.copy()
        H = _H.copy()

        prev_err = err
        err = residuo(A, W, H)
        end = time.time()
        print("=> Iteration {0} Abs Error {1:.3f} Delta {2:.5f} Time {3:.3f}"
              .format(i, err, (prev_err-err)/err, end - begin))

    return (W, H)


def matriz_arquivo(arquivo, n_train=-1):
    """
    Lê [arquivo] e transforma em array Matriz normalizarda

        :param arquivo: string Nome do arquivo

        :return: ndarray Matriz extraída do arquivo normalizarda
    """

    matriz = []

    with open(arquivo, "r+") as arq:
        for raw_linha in arq:
            raw_linha = raw_linha.strip('\n').split(' ')
            linha = [float(num) for num in raw_linha[:n_train]]
            matriz.append(linha)
    return np.array(matriz)/255.0


def treinar(d, p=10, n_train=100):
    '''
    Executar treinamento do dígito d, gerando a matriz Wd

        :param d: int Dígito d
        :param p: int Quantidade de linhas da matriz W
        :param n_train: int Quantidade de amostras a serem utilizadas
        no treinamento

        :return None:
    '''

    folder = "./{0}-{1}/".format(n_train, p)

    A = matriz_arquivo('dados_mnist/train_dig'+str(d)+'.txt', n_train)
    print("Matrix loaded!")
    n, m = A.shape
    Wd, H = fatorar_wh(A, p)

    if not os.path.isdir(folder):
        os.mkdir(folder)

    np.save('./{0}/W{1}.npy'.format(folder, d), Wd)

    return Wd


def fatorar_digito(d, n_test=1000, n_train=100, p=5):
    '''
    Calcular H a partir de Wd e retorna o erro de A-Wd*H

        :param Wd: ndarray Matriz Wd para certo d
        :param n_test: int Quantidade de amostras utilizadas no teste
        :param n_train: int Quantidade de amostras utilizadas no treinamento
        :param p: int Quantidade de linhas da matriz W

        :return c: double Módulo da matriz A-Wd*H
    '''

    print("[LOG] Testing digit {0}".format(d))
    A = matriz_arquivo("./dados_mnist/test_images.txt", n_test)
    n, n_test_ = A.shape
    if n_test != n_test_:
        raise ValueError("A leitura de A não foi correta ")

    Wd = np.load("./{0}-{1}/W{2}.npy".format(n_train, p, d))

    H = resolver_sist(Wd.copy(), A.copy())

    p, k = H.shape
    if k != n_test:
        raise ValueError("Matriz H não é compatível")

    res = A - np.matmul(Wd, H)
    c = np.zeros(n_test)
    c[:] = np.sum(np.power(res, 2), axis=0)

    return c


def classificar(n_test=1000, n_train=100, p=5):
    '''
    Classificar as imagens do dataset de teste

        :param n_test: int Quantidade de amostras utilizadas no teste
        :param n_train: int Quantidade de amostras utilizadas no treinamento
        :param p: int Quantidade de linhas da matriz W

        :return d: ndarray (n_test,) Digito calculado para cada imagem
    '''

    folder = "./estimador/"

    if not os.path.isdir(folder):
        os.mkdir(folder)

    digito = np.zeros(n_test)
    menor_erro = np.zeros(n_test)

    erros = np.empty((n_test, 10))

    raw_err = [fatorar_digito(d, n_test=n_test, n_train=n_train, p=p)
               for d in range(10)]

    err = np.array(raw_err)

    erros = err.transpose().copy()

    digito[:] = np.argmin(erros, axis=1)
    np.save(folder + "{0}-{1}.npy".format(n_train, p), digito)
    # analisar(digito)
    return digito


def analisar(estimativa, n_test=1000):
    '''
    Analisar o percentual de acertos

        :param estimativa: ndarray Lista de classificações do conjunto de teste
        :param n_test: int Quantidade de amostras utilizadas no teste

        :return acertos: list Lista de acertos para cada digito
        :return permil: list Fração de acertos relativa a cada mil amostras
    '''

    '''
    Leitura dos índices corretos a partir de test_index.txt
    '''

    index = np.zeros(n_test).astype(np.int)
    with open("./dados_mnist/test_index.txt", "r+") as arq:
        for i, raw_linha in enumerate(arq):
            if i < n_test:
                index[i] = int(raw_linha.strip('\n'))
            else:
                break

    total = np.array([np.sum(index == i) for i in range(10)])

    acertos = np.zeros(10)
    for resp, gaba in zip(estimativa, index):
        if int(gaba) == int(resp):
            acertos[int(gaba)] += 1

    permil = (acertos*1000/total)

    return acertos, permil
