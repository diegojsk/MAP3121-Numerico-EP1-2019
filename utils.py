import matplotlib.pyplot as plt
import numpy as np


def load_data(dig=0, model_path="./dados_mnist/train_dig{0}.txt", col_i=0, cols=1):

    path = model_path.format(dig)
    output = np.empty((784, cols))
    with open(path, 'r') as data:
        for i, line in enumerate(data):
            interest = line.strip('\n').split(' ')[col_i: col_i+cols]
            output[i, :] = [int(num) for num in interest]
    return output.reshape((28, 28, cols))

LINS = 2
COLS = 4

A = load_data(7, cols=LINS*COLS)

fig = plt.figure()
axs = fig.subplots(LINS, COLS)

for i, ax in enumerate(axs.flatten()):
    ax.imshow(-A[:, :, i], cmap='Greys')

plt.show()
