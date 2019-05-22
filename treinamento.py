import argparse
from main import *
import numpy as np
import time

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Train Wd matrix")

    parser.add_argument('digit', type=int, help='Digit to start training')
    parser.add_argument('n_train', type=int,
                        help='Amount of images to use from training dataset')
    parser.add_argument('p', type=int, help='Amount of lines of W matrix')

    args = parser.parse_args()

    begin = time.time()
    for d in range(args.digit, 10):
        partial = time.time()
        print("[LOG] Loading file ({}/10)".format(d+1s))
        treinar(d, p=args.p, n_train=args.n_train)
        end = time.time()
        print("[LOG] Finished file {0} in {1} s".format(d+1, end - partial))
        print("[LOG] Total elapsed time {0} in {1} s".format(d+1, end - begin))
