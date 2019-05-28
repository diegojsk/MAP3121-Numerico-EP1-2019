[![forthebadge](https://forthebadge.com/images/badges/built-with-science.svg)](https://forthebadge.com) [![forthebadge](https://forthebadge.com/images/badges/made-with-python.svg)](https://forthebadge.com)


# Numerico-EP1-2019

Exercício Programa 1 para a disciplina MAP 3121 - Métodos Numéricos e Aplicações
Poli - USP - 2019

[Link](https://www.ime.usp.br/~map3121/2019/map3121/programas/EP1-MachineLearning_v2.pdf) do enunciado

[Link](https://docs.google.com/document/d/1__LbmVL0IIN13Hf8ZoNnkv6He_3kZvP12aE3JNr4U-k/edit) do relatório EP1

## Estrutura dos arquivos

 - __main.py__ - Biblioteca com as funções necessárias para a resolução do EP
 - __tarefa_1.py__ - Execução dos casos de teste da tarefa 1
 - __tarefa_2.py__ - Execução dos casos de teste da tarefa 2
 - __train.py__ - Script de treinamento para a obtenção das matrizes Wd.  
As matrizes obtidas são salvas na forma de arquivos .npy dentro de pastas na raiz do projeto.
As pastas são nomeadas na forma ./[n_train]-[p]/
   - Uso ```python3 train.py digit n_train p```
   - **digit** - Dígito para começar o treinamento
   - **n_train** - Quantidade de imagens a serem uitlizadas no treinamento
   - **p** - Quantidade de linhas da matriz Wd
   - i.e. ```python3 train.py 0 4000 15```
 - __test.py__ - Classifica *n_test* imagens do conjunto de teste e afere o desempenho do modelo  
O script acessa as matrizes salvas nas pastas na raiz do projeto, de acordo com os 
parâmetros *n_train* e *p* especificados na chamada do script
   - Uso ```python3 test.py n_train p```
   - **n_train** - Quantidade de imagens uitlizadas no treinamento
   - **p** - Quantidade de linhas da matriz Wd
   - i.e. ```python3 test.py 4000 15```
 - __utils.py__ - Script para visualização das imagens na forma de vetor coluna


## Estrutura da main.py

 - def calc_c(a, b)
   - Calcula o parâmetro c, correspondente ao cosseno do ângulo da rotação
    de Givens

        - :param a: w[i,k] - elemento da matriz W na posição (i, k)
        - :param b: w[j,k] - elemento da matriz W na posição (j, k)

        - :return: Cosseno do ângulo de rotação

 - def calc_s(a, b)
   - Calcula o parâmetro s, correspondente ao seno do ângulo da rotação
    de Givens

        - :param a: double w[i,k] - elemento da matriz W na posição (i, k)
        - :param b: double w[j,k] - elemento da matriz W na posição (j, k)

        - :return: Seno do ângulo de rotação

 - def rot_givens(W, n, m, i, j, c, s)
   - Efetua a rotação de Givens para matriz W
    ! Atenção, as alterações são feitas in place na matriz W

        - :param W: ndarray
        - :param n: int Número de linhas de W
        - :param m: int Número de colunas de W
        - :param i: int Linha a ser rotacionada
        - :param j: int Linha a ser rotacionada
        - :param c: double Cosseno do ângulo de rotação
        - :param s: double Seno do ângulo de rotação

        - :return: None

 - def zerar_elemento(W, i, j, k)    
   - Realiza uma rotação de Givens de modo a zerar o elemento (j,k)

        - :param W: ndarray
        - :param i: linha a ser rotacionada
        - :param j: linha a ser zerada
        - :param k: coluna a ser zerada

        - :return: None

 - def fatorar_qr(W)
   - Encontra a matriz R de modo que Q*R = W, onde Q é o resultado de
    sucessivas matrizes de rotação de Givens e R é uma matriz triangular
    superior.  
    Ou seja, a função transforma a matriz W em uma matriz triangular
    superior por meio de sucessivas rotações de Givens.

        - :param W: ndarray shape(n, m)

        - :return: None

 - def resolver_sist(W, A)

   - Dadas matrizes W e A, encontra a matriz H, tal que W*H = A  
    Função Principal da Primeira Tarefa

        - :param W: ndarray shape(n, m)
        - :param A: ndarray shape(n, p)

 - def residuo(A, W, H)

   - Calcula o resíduo quadrático para | A - W*H |

        - :param A: ndarray shape(n, p)
        - :param W: ndarray shape(p, m)
        - :param W: ndarray shape(n, m)

        - :return erro: Erro quadrático para a fatoração

 - def normalizar(M)

   - Normaliza a matriz M, de modo que a norma de todas as suas colunas
    seja igual a um

        - :param M: ndarray

        - :return: None

 - def fatorar_wh(A, p)

   - Encontra uma fatoração para a matriz A na forma A = W*H de modo que
    W e H são matrizes não-negativas  
    Função Principal da Segunda Tarefa

        - :param A: ndarray Matriz a ser fatorada
        - :param p: int Quantidade de colunas de W

        - :return W: ndarray Matriz W da fatoração não negativa de A
        - :return H: ndarray Matriz H da fatoração não negativa de A

 - def matriz_arquivo(arquivo, n_train=-1)

   - Lê [arquivo] e transforma em array Matriz normalizarda

        - :param arquivo: string Nome do arquivo

        - :return: ndarray Matriz extraída do arquivo normalizarda

 - def treinar(d, p=10, n_train=100)

   - Executar treinamento do dígito d, gerando a matriz Wd

        - :param d: int Dígito d
        - :param p: int Quantidade de linhas da matriz W
        - :param n_train: int Quantidade de amostras a serem utilizadas
        no treinamento

        - :return None:

 - def fatorar_digito(d, n_test=1000, n_train=100, p=5)

   - Calcular H a partir de Wd e retorna o erro de A-Wd*H

        - :param Wd: ndarray Matriz Wd para certo d
        - :param n_test: int Quantidade de amostras utilizadas no teste
        - :param n_train: int Quantidade de amostras utilizadas no treinamento
        - :param p: int Quantidade de linhas da matriz W

        - :return c: double Módulo da matriz A-Wd*H

 - def classificar(n_test=1000, n_train=100, p=5)

   - Classificar as imagens do dataset de teste

        - :param n_test: int Quantidade de amostras utilizadas no teste
        - :param n_train: int Quantidade de amostras utilizadas no treinamento
        - :param p: int Quantidade de linhas da matriz W

        - :return d: ndarray (n_test,) Digito calculado para cada imagem

 - def analisar(estimativa, n_test=1000)

   - Analisar o percentual de acertos

        - :param estimativa: ndarray Lista de classificações do conjunto de teste
        - :param n_test: int Quantidade de amostras utilizadas no teste

        - :return acertos: list Lista de acertos para cada digito
        - :return permil: list Fração de acertos relativa a cada mil amostras

