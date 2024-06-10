# sklearn benchmark
This benchmark runs a series of scripts that trains a model from sklearn (Scikit-Learn). I got the series of scripts via decomposing the sklearn source code by hand. [Original](https://github.com/scikit-learn/scikit-learn/blob/289326704e13f7a5bf4c6c594c038051e968e1fd/sklearn/linear_model/_logistic.py)

## Purpose
I think this benchmark shows two things for a system like hS - viability in AI workflows and correctness. The first is quite self explanatory. If hS can run this benchmark, then it has proven that hS can handle the task of gluing together a nontrivial ML training workflow.
The second is correctness. There is a very clear ground truth (the model trained by pure sklearn) to compare hS's output to. Assuming the random seeds are set to the same value across the board, hS's model should produce the *exactly* the same weights as the ground truth.

## Usage
Running fit.sh will generate temporary files in a ./tmp folder

Before running, the user need to install the packages (possibly in a
virtual environment) by `pip install -r requirements.txt`, and make
sure the result direcotry exists (`mkdir -p result`). Then run
`run.sh` with appropriate environment where python is aliased to the
correct python3 installation (a.k.a. in a virtual environment).

To parallelize, we want one-vs-rest classification, where we generate multiple models.
Additionally, the forest cover dataset has much more samples than it has features.
This makes the Newton-Cholesky solver ideal for this task.

