from main import classificar
import argparse
import time

train = [100, 1000, 4000]
ps = [5, 10, 15]


#         classificar(n_train=n_train, p=p)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Test Wd matrix")

    # parser.add_argument('n_train', type=int, help='Digit to be trained')
    parser.add_argument('p', type=int, help='Digit to be trained')

    args = parser.parse_args()

    # N_TRAIN = args.n_train
    P = args.p

    begin = time.time()
    # for p in ps:
    for n_train in train:
        N_TRAIN = n_train
        # P = p
        partial = time.time()
        print("[LOG] Wd {0}-{1}".format(N_TRAIN, P))
        classificar(n_train=N_TRAIN, p=P)
        end = time.time()
        print("[LOG] Finished in {0} s".format(end - partial))
        print("[LOG] Total elapsed time {0} s".format(end - begin))
