# sklearn benchmark
This benchmark runs a series of scripts that trains a model from sklearn (Scikit-Learn). I got the series of scripts via decomposing the sklearn source code by hand.

## Purpose
I think this benchmark shows two things for a system like hS - viability in AI workflows and correctness. The first is quite self explanatory. If hS can run this benchmark, then it has proven that hS can handle the task of gluing together a nontrivial ML training workflow.
The second is correctness. There is a very clear ground truth (the model trained by pure sklearn) to compare hS's output to. Assuming the random seeds are set to the same value across the board, hS's model should produce the *exactly* the same weights as the ground truth.

## Usage
Running fit.sh will generate temporary files in a ./tmp folder

To parallelize, we want one-vs-rest classification, where we generate multiple models.
Additionally, the forest cover dataset has much more samples than it has features.
This makes the Newton-Cholesky solver ideal for this task.

