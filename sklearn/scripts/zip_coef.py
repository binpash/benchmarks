import sys
import numpy as np
import pickle
import os

fold_coefs_ = []

tmp = os.environ.get('TMP')

with open(sys.argv[1], 'r+b') as file1:
    model = pickle.load(file1)

    for i in range(1, 8):
        filepath = os.path.join(tmp, f'result_{i}.obj')
        with open(filepath, 'r+b') as file2:
            fold_coefs_.append(pickle.load(file2))

    fold_coefs_, _, n_iter_ = zip(*fold_coefs_)

    model.n_iter = np.asarray(n_iter_, dtype=np.int32)[:, 0]
    pickle.dump(model, file1)

filepath = os.path.join(tmp, 'fold_coef.obj')
with open(filepath, 'w+b') as file:
    pickle.dump(fold_coefs_, file)
