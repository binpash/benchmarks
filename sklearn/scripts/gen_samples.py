#!/usr/bin/env python3

import argparse
import os
import pickle
from sklearn.model_selection import train_test_split
from sklearn import datasets
from imblearn.over_sampling import SMOTE
import numpy as np


parser = argparse.ArgumentParser()
parser.add_argument('--small', action='store_true', help='Use a small subset of the data')
parser.add_argument('--min', action='store_true', help='Use a very small subset (~1MB) of the data')
args = parser.parse_args()

raw_data = datasets.fetch_covtype(data_home="inputs", download_if_missing=False)
X, y = raw_data.data, raw_data.target

if args.min:
    target_rows = int(1e6 / (X.shape[1] * 8))
    X, _, y, _ = train_test_split(X, y, train_size=target_rows, stratify=y, random_state=0)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=0)

elif args.small:
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=0)
else:
    n_target = int(8e8 / X.shape[1] / 8)
    if n_target > len(X):
        smote = SMOTE(random_state=0)
        X, y = smote.fit_resample(X, y)
        while len(X) < n_target:
            X_extra, y_extra = smote.fit_resample(X, y)
            X = np.vstack((X, X_extra))
            y = np.hstack((y, y_extra))
        X, y = X[:n_target], y[:n_target]
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=0)

filenames = ['X_train', 'X_test', 'y_train', 'y_test']
tmp = os.environ.get('TMP')
data = [X_train, X_test, y_train, y_test]

filepath = os.path.join(tmp, 'model.obj')
for datum, name in zip(data, filenames):
    filepath = os.path.join(tmp, f'{name}.obj')
    with open(filepath, 'w+b') as file:
        pickle.dump(datum, file)
