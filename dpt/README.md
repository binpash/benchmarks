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
- “Digital Pyramids Text Project,” 2025. [https://example.org/dpt](http://dpt.cs.brown.edu/)
- [https://arxiv.org/abs/2304.02643](https://arxiv.org/abs/2304.02643)
