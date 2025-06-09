#!/usr/bin/env python3

import os
import sys
import pickle
import numpy as np
from sklearn.linear_model import _logistic


def main():
    model_path = sys.argv[1]
    X_path = sys.argv[2]
    y_path = sys.argv[3]
    C_path = sys.argv[4]
    warm_coef_path = sys.argv[5]
    max_sq_sum_path = sys.argv[6]
    multi_class = sys.argv[7]
    penalty = sys.argv[8]
    class_ = int(sys.argv[9])

    # Load non-array objects
    with open(model_path, "rb") as f:
        model = pickle.load(f)
    with open(C_path, "rb") as f:
        C_ = pickle.load(f)
    with open(warm_coef_path, "rb") as f:
        warm_start_coef_ = pickle.load(f)
    with open(max_sq_sum_path, "rb") as f:
        max_squared_sum = pickle.load(f)

    X = np.load(X_path)
    y = np.load(y_path)

    result = _logistic._logistic_regression_path(
        X,
        y,
        pos_class=class_,
        Cs=[C_],
        l1_ratio=model.l1_ratio,
        fit_intercept=model.fit_intercept,
        tol=model.tol,
        verbose=model.verbose,
        solver=model.solver,
        multi_class=multi_class,
        max_iter=model.max_iter,
        class_weight=model.class_weight,
        check_input=False,
        random_state=model.random_state,
        coef=None,
        penalty=penalty,
        max_squared_sum=max_squared_sum,
        sample_weight=None,
    )

    out_dir = os.environ.get("OUT", ".")
    out_path = os.path.join(out_dir, f"result_{class_}.obj")
    with open(out_path, "wb") as f:
        pickle.dump(result, f)


if __name__ == "__main__":
    main()
