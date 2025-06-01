#!/usr/bin/env python3

import argparse
import os
import tempfile
from collections import Counter

import pickle
import numpy as np
from imblearn.over_sampling import SMOTE
from sklearn.datasets import fetch_covtype
from sklearn.model_selection import train_test_split


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


def main():
    parser = argparse.ArgumentParser()
    g = parser.add_mutually_exclusive_group()
    g.add_argument("--small", action="store_true", help="20 % test split of original")
    g.add_argument("--min", action="store_true", help="≈ 2 k rows for sanity checks")
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



    out_dir = os.getenv("TMP", tempfile.gettempdir())
    os.makedirs(out_dir, exist_ok=True)

    filenames = ['X_train', 'X_test', 'y_train', 'y_test']
    tmp = os.environ.get('TMP')
    data = [X_train, X_test, y_train, y_test]

    for datum, name in zip(data, filenames):
        filepath = os.path.join(tmp, f'{name}.obj')
        with open(filepath, 'w+b') as file:
            pickle.dump(datum, file)


if __name__ == "__main__":
    main()
