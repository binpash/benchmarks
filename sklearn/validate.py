#!/usr/bin/env python3
# ‑‑ minimal changes: only flag parsing + matching data‑sizing logic

import argparse
import os
import pickle
import numpy as np
from sklearn import datasets
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from imblearn.over_sampling import SMOTE

parser = argparse.ArgumentParser()
parser.add_argument('--small', action='store_true', help='Use a small subset of the data')
parser.add_argument('--min',   action='store_true', help='Use a very small subset (~1 MB)')
args = parser.parse_args()

max_iter = 100

raw = datasets.fetch_covtype(data_home="inputs", download_if_missing=False)
X, y = raw.data, raw.target

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

control_model = LogisticRegression(max_iter=max_iter,
                                   solver='newton-cholesky',
                                   multi_class='ovr')
control_model.fit(X_train, y_train)
control_score = control_model.score(X_test, y_test)

# fall back to ./tmp if $TMP is unset, mirroring the training script
tmp_dir = os.environ.get('TMP', 'tmp')
model_path = os.path.join(tmp_dir, 'trained_model.obj')

with open(model_path, 'rb') as fh:
    experiment_model = pickle.load(fh)

experiment_score = experiment_model.score(X_test, y_test)

# strict equality as before
assert experiment_score == control_score
assert np.array_equal(control_model.coef_, experiment_model.coef_)
