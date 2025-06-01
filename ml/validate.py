#!/usr/bin/env python3
# ‑‑ minimal changes: only flag parsing + matching data‑sizing logic

import argparse
import os
import pickle
import numpy as np
from sklearn import datasets
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.datasets import fetch_covtype
from imblearn.over_sampling import SMOTE


def balanced_sample(X, y, target_rows, *, random_state=0):
    smote = SMOTE(random_state=random_state)
    X_bal, y_bal = smote.fit_resample(X, y)

    if target_rows <= len(X_bal):
        sel = np.random.RandomState(random_state).choice(
            len(X_bal), target_rows, replace=False
        )
        return X_bal[sel], y_bal[sel]

    def extra_samples(n):
        rng = np.random.RandomState(random_state + 1)
        idx = rng.choice(len(X_bal), n, replace=True)
        return X_bal[idx], y_bal[idx]

    need = target_rows - len(X_bal)
    X_extra, y_extra = extra_samples(need)
    X_big = np.vstack((X_bal, X_extra))
    y_big = np.hstack((y_bal, y_extra))
    return X_big, y_big


parser = argparse.ArgumentParser()
g = parser.add_mutually_exclusive_group()
g.add_argument("--small", action="store_true", help="20 % test split of original")
g.add_argument("--min", action="store_true", help="≈ 2 k rows for sanity checks")
parser.add_argument(
    "--rows",
    type=int,
    default=None,
    help="Exact row count to generate (implies balanced oversampling)",
)
args = parser.parse_args()

X, y = fetch_covtype(data_home="inputs", download_if_missing=True, return_X_y=True)

if args.min:
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=0
    )

elif args.small:
    n_target = int(5e9 / (X.shape[1] * 8))
    X_bal, y_bal = balanced_sample(X, y, n_target)
    X_train, X_test, y_train, y_test = train_test_split(
        X_bal, y_bal, test_size=0.2, random_state=0
    )

else:
    n_target = int(10e9 / (X.shape[1] * 8))
    X_bal, y_bal = balanced_sample(X, y, n_target)
    X_train, X_test, y_train, y_test = train_test_split(
        X_bal, y_bal, test_size=0.2, random_state=0
    )

control_model = LogisticRegression(
    max_iter=1000, solver="newton-cholesky", multi_class="ovr"
)
control_model.fit(X_train, y_train)
control_score = control_model.score(X_test, y_test)

tmp_dir = os.environ.get("OUT")
model_path = os.path.join(tmp_dir, "trained_model.obj")
with open(model_path, "rb") as fh:
    experiment_model = pickle.load(fh)

experiment_score = experiment_model.score(X_test, y_test)

assert experiment_score == control_score
assert np.array_equal(control_model.coef_, experiment_model.coef_)
