import matplotlib.pyplot as plt
import numpy as np


def load_data(dig=0, model_path="./dados_mnist/train_dig{0}.txt", col_i=0,
              cols=1):

    path = model_path.format(dig)
    output = np.empty((784, cols))
    with open(path, 'r') as data:
        for i, line in enumerate(data):
            interest = line.strip('\n').split(' ')[col_i: col_i+cols]
            output[i, :] = [int(num) for num in interest]
    return output.reshape((28, 28, cols))

LINS = 2
COLS = 4

# A = load_data(7, model_path='./output/matmul-4-4000-10.npy', cols=LINS*COLS)

<<<<<<< HEAD
A = np.load('./output/matmul-4-4000-10.npy')[:,:LINS*COLS]
=======
A = np.load('./output/matmul-1-4000-10.npy')[:, :LINS*COLS]
>>>>>>> 338cc2a48db2618581804354540d12c958d35a9e
A = A.reshape((28, 28, LINS*COLS))

fig = plt.figure()
axs = fig.subplots(LINS, COLS)

n_test = LINS*COLS
index = np.zeros(n_test).astype(np.int)
with open("./dados_mnist/test_index.txt", "r+") as arq:
    for i, raw_linha in enumerate(arq):
        if i < n_test:
            index[i] = int(raw_linha.strip('\n'))
        else:
            break

for i, ax in enumerate(axs.flatten()):
    ax.imshow(-A[:, :, i], cmap='Greys')
    ax.set_title("Digit: {}".format(index[i]))

plt.tight_layout()
plt.show()
