#!/usr/bin/env python3

from sklearn.linear_model import _logistic
import numpy as np
import sys
import os
import pickle
X = np.load(sys.argv[1])
max_squared_sum = _logistic.row_norms(X, squared=True).max()

out_dir = os.environ.get("OUT", ".")
with open(os.path.join(out_dir, "max_squared_sum.obj"), "wb") as f:
    pickle.dump(max_squared_sum, f)