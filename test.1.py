from main import classificar, analisar
import numpy as np
import argparse
import time

train = [1000, 4000]
ps = [5, 10, 15]


if __name__ == "__main__":

    np.set_printoptions(precision=3, suppress=True)

    output = open('resultados.csv', 'w')

    for n_train in train:
        for p in ps:
            digits = classificar(n_train=n_train, p=p)

            T, A = analisar(digits)
            print("{}-{}".format(n_train, p), end=', ')
            output.write("{}-{}, ".format(n_train, p))
            for d, i in enumerate(A):
                print("{1:.2f}%".format(d, i/10), end=', ')
                output.write("{1:.2f}%, ".format(d, i/10))
            print("{}%\n".format(np.sum(T)/10))
            output.write("{}%\n".format(np.sum(T)/10))

    output.close()
