"Executar treinamentos dos digitos"

from main import *
import numpy as np
import time

begin = time.time()
for d in range(5, 10):
    partial = time.time()
    print("[LOG] Loading file ({}/10)".format(d+1))
    treinamento(d, p=10, ndig_treino=1000)
    end = time.time()
    print("[LOG] Finished file {0} in {1} s".format(d+1, end - partial))
    print("[LOG] Total elapsed time {0} in {1} s".format(d+1, end - begin))