"Executar treinamentos dos digitos"

from main import *
import numpy as np
import time

begin = time.time()
for d in range(10):
    partial = time.time()
    print("Loading file ({}/10)".format(d+1))
    treinamento(d)
    end = time.time()
    print("Finished file {0} in {1} s".format(d+1, end - partial))
    print("Total elapsed time {0} in {1} s".format(d+1, end - partial))