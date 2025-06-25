#!/usr/bin/env python3

import sys
import numpy as np
import pickle
from sklearn.linear_model import LogisticRegression

def main():
    model_path = sys.argv[1]
    X_path = sys.argv[2]
    y_path = sys.argv[3]

    X = np.load(X_path)
    y = np.load(y_path)

    with open(model_path, "rb") as model_file:
        model: LogisticRegression = pickle.load(model_file)

    try:
        model._validate_data(
            X,
            y,
            accept_sparse="csr",
            dtype=np.float64 if model.solver == "lbfgs" else [np.float64, np.float32],
            order="C",
            accept_large_sparse=model.solver not in ["liblinear", "sag", "saga"],
        )
        exit(0)

    except ValueError:
        exit(1)


if __name__ == "__main__":
    main()
