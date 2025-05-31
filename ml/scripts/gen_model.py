#!/usr/bin/env python3

from sklearn.linear_model import LogisticRegression
import pickle
import sys
import os

reg = LogisticRegression(max_iter=int(sys.argv[1]), 
                         solver='newton-cholesky', 
                         multi_class='ovr')

tmp = os.environ.get('TMP')
filepath = os.path.join(tmp, 'model.obj')
with open(filepath, 'w+b') as file:
    pickle.dump(reg, file)