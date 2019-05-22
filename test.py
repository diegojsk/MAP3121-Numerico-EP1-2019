from main import classificar, analisar
import numpy as np
import argparse
import time

train = [100, 1000, 4000]
ps = [5, 10, 15]


if __name__ == "__main__":

    np.set_printoptions(precision=3, suppress=True)

    parser = argparse.ArgumentParser(description="Test Wd matrix")

    parser.add_argument('n_train', type=int,
                        help='Amount of images used from training dataset')
    parser.add_argument('p', type=int, help='Amount of lines of W matrix')

    args = parser.parse_args()

    digits = classificar(n_train=args.n_train, p=args.p)

    T, A = analisar(digits)
    print("N_train = {} P = {}".format(args.n_train, args.p))
    for d, i in enumerate(A):
        print("{0}: {1:.2f}%".format(d, i/10))
    print("Total: {}%\n".format(np.sum(T)/10))
