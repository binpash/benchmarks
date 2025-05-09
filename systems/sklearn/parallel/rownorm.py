from sklearn.linear_model import _logistic
import sys
import pickle
import os

with open(sys.argv[1], 'rb') as file:
    X = pickle.load(file)
    
max_squared_sum = _logistic.row_norms(X, squared=True).max()

tmp = os.environ.get('TMP')
filepath = os.path.join(tmp, 'max_squared_sum.obj')
with open(filepath, 'w+b') as file:
    pickle.dump(max_squared_sum, file)
