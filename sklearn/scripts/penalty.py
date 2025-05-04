#!/usr/bin/env python3

import warnings
import sys
import pickle
import numpy as np
import os

with open(sys.argv[1], 'rb') as file:
    model = pickle.load(file)

if model.penalty != "elasticnet" and model.l1_ratio is not None:
    warnings.warn(
        "l1_ratio parameter is only used when penalty is "
        "'elasticnet'. Got "
        "(penalty={})".format(model.penalty)
    )

if model.penalty is None or model.penalty == "none":
    if model.C != 1.0:  # default values
        warnings.warn(
            "Setting penalty=None will ignore the C and l1_ratio parameters"
        )
        # Note that check for l1_ratio is done right above
    C_ = np.inf
    penalty = "l2"
else:
    C_ = model.C
    penalty = model.penalty

tmp = os.environ.get('TMP')
filepath = os.path.join(tmp, 'C_.obj')
with open(filepath, 'w+b') as file:
    pickle.dump(C_, file)
print(penalty)
