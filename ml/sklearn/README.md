# sklearn

This benchmark runs a series of scripts that trains a model from sklearn (Scikit-Learn).

## sklearn

This benchmark replicates the training process of a Scikit-Learn `LogisticRegression` model using manually decomposed steps derived from Scikit-Learn’s internal source code.

### Inputs

The input dataset is the `Covtype` dataset from the `sklearn.datasets` module, which is a synthetic dataset with 58.1012 samples and 54 features.
Input size is inflated using the `SMOTE` algorithm to create a balanced dataset with 800M samples.

### Running

The benchmark simulates the entire `fit` pipeline using individual Python scripts that replicate Scikit-Learn internals:
- Generates training/test data.
- Constructs model, computes row norms, warm-start coefficients.
- Runs coefficient calculation in parallel over class labels.
- Assembles and stores the final model in `result/trained_model.obj`.

### References

- https://archive.ics.uci.edu/dataset/31/covertype
