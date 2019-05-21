import argparse
from main import *
import numpy as np
import time

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Train Wd matrix")

    parser.add_argument('digit', type=int, help='Digit to be trained')
    parser.add_argument('n_dig', type=int, help='Digit to be trained')
    parser.add_argument('p', type=int, help='Digit to be trained')
    parser.add_argument('-s', '--sequence', type=int, default=1,
                        help='How many digits to train')

    args = parser.parse_args()

    begin = time.time()
    for d in range(args.digit, args.digit + args.sequence):
        partial = time.time()
        print("[LOG] Loading file ({}/10)".format(d))
        treinar(d, p=args.p, ndig_treino=args.n_dig)
        end = time.time()
        print("[LOG] Finished file {0} in {1} s".format(d+1, end - partial))
        print("[LOG] Total elapsed time {0} in {1} s".format(d+1, end - begin))
