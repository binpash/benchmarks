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

## dpt

This benchmark segments and classifies hieroglyphs in a set of input images using pre-trained models.

### Inputs

- `inputs/images/`: A directory containing image files to be processed.

### Running

Each image is processed in two steps:

1. Segmentation:
   -  Regions of interest are extracted from the image using a pre-trained segmentation model.
  
2. Classification:
   -  Regions of interest are extracted from the image using a pre-trained segmentation model.

Results are written to a single output file and include:
- Image filename
- Bounding box coordinates for each detected region
- Hieroglyph classification labels
- Confidence scores for each prediction
All outputs are deterministic across runs.

### Validation
Correctness is validated by comparing each generated output file against a reference stored in the hashes/ directory.

### References
- “Digital Pyramids Text Project,” 2025. [http://dpt.cs.brown.edu/](http://dpt.cs.brown.edu/)
- [https://arxiv.org/abs/2304.02643](https://arxiv.org/abs/2304.02643)
