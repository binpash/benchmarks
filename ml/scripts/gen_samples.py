#!/usr/bin/env python3
"""
Generate (optionally up-sampled) train / test splits of the Covertype data set
without exceeding RAM, even when the target size is in the multi-GB range.
"""

import argparse
import os
import tempfile
import numpy as np
from imblearn.over_sampling import SMOTE
from numpy.lib.format import open_memmap
from sklearn.datasets import fetch_covtype
from sklearn.model_selection import train_test_split


def balanced_sample_chunked(
    X,
    y,
    target_rows,
    *,
    chunk_size: int = 200_000,
    random_state: int = 0,
    work_dir: str | None = None,
):
    if target_rows <= len(X):
        rng = np.random.default_rng(random_state)
        sel = rng.choice(len(X), target_rows, replace=False)
        return X[sel], y[sel]

    if work_dir is None:
        work_dir_obj = tempfile.TemporaryDirectory()
        work_dir = work_dir_obj.name
    else:
        work_dir_obj = None

    n_features = X.shape[1]
    X_mm = np.memmap(
        os.path.join(work_dir, "X_balanced.dat"),
        dtype=X.dtype,
        mode="w+",
        shape=(target_rows, n_features),
    )
    y_mm = np.memmap(
        os.path.join(work_dir, "y_balanced.dat"),
        dtype=y.dtype,
        mode="w+",
        shape=(target_rows,),
    )

    rng = np.random.default_rng(random_state)
    smote = SMOTE(random_state=random_state)

    filled = 0
    progress_every = max(target_rows // 20, chunk_size)

    while filled < target_rows:
        need = target_rows - filled
        this_chunk = min(chunk_size, need)

        idx = rng.choice(len(X), this_chunk, replace=True)
        X_chunk, y_chunk = X[idx], y[idx]

        X_s, y_s = smote.fit_resample(X_chunk, y_chunk)

        take = min(len(X_s), target_rows - filled)
        assert take > 0, f"SMOTE returned too few samples ({len(X_s)}); filled={filled}, target_rows={target_rows}"
        sel = rng.choice(len(X_s), take, replace=False)

        X_mm[filled : filled + take] = X_s[sel]
        y_mm[filled : filled + take] = y_s[sel]
        filled += take

        if filled % progress_every == 0 or filled == target_rows:
            print(f"[balanced_sample_chunked] {filled:,}/{target_rows:,} rows ready")

    del X_mm, y_mm
    X_bal = np.memmap(
        os.path.join(work_dir, "X_balanced.dat"),
        dtype=X.dtype,
        mode="r",
        shape=(target_rows, n_features),
    )
    y_bal = np.memmap(
        os.path.join(work_dir, "y_balanced.dat"),
        dtype=y.dtype,
        mode="r",
        shape=(target_rows,),
    )

    if work_dir_obj is not None:
        X_bal._tmp_dir_keeper = work_dir_obj
        y_bal._tmp_dir_keeper = work_dir_obj

    return X_bal, y_bal



def atomic_save(arr: np.ndarray, path: str, *, block: int = 100_000) -> None:
    tmp = f"{path}.tmp"
    fp = open_memmap(tmp, mode="w+", dtype=arr.dtype, shape=arr.shape)
    step = slice(None) if arr.ndim == 1 else (slice(None), ...)
    for start in range(0, len(arr), block):
        end = start + block
        fp[start:end][step] = arr[start:end]
    fp.flush()
    assert fp.shape == arr.shape, f"Mismatch saving array to {path}"
    os.replace(tmp, path)


parser = argparse.ArgumentParser()
g = parser.add_mutually_exclusive_group()
g.add_argument("--small", action="store_true", help="≈5 GB up-sampled dataset")
g.add_argument("--min", action="store_true", help="No oversampling – use original rows")
parser.add_argument(
    "--rows",
    type=int,
    default=None,
    help="Exact row count to generate (overrides preset sizes)",
)
args = parser.parse_args()

X, y = fetch_covtype(data_home="inputs", download_if_missing=True, return_X_y=True)

if args.min:
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=0
    )
else:
    scale_bytes = 5_000_000_000 if args.small else 10_000_000_000
    n_target = int(scale_bytes / (X.shape[1] * 8))
    if args.rows is not None:
        n_target = args.rows

    X_bal, y_bal = balanced_sample_chunked(
        X,
        y,
        n_target,
        random_state=0,
        work_dir=os.environ.get("TMP"),
    )
    X_train, X_test, y_train, y_test = train_test_split(
        X_bal, y_bal, test_size=0.2, random_state=0
    )

out_dir = os.getenv("TMP", tempfile.gettempdir())
os.makedirs(out_dir, exist_ok=True)

for name, arr in {
    "X_train": X_train,
    "X_test": X_test,
    "y_train": y_train,
    "y_test": y_test,
}.items():
    atomic_save(arr, os.path.join(out_dir, f"{name}.npy"))

del X_train, X_test, y_train, y_test
