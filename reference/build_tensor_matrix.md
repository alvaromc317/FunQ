# Build Tensor Design Matrix

Inner function to construct the tensor product design matrix required
for functional data modeling. It efficiently scales the B-spline basis
evaluations by the subject-specific functional principal component
scores across all observed time points.

## Usage

``` r
build_tensor_matrix(scores, Y.mask, spline.basis, intercept = TRUE)
```

## Arguments

- scores:

  A matrix of functional Principal Component Analysis (fPCA) scores for
  each subject.

- Y.mask:

  A logical mask matrix of dimensions (N x T) indicating which time
  points are observed (`TRUE`) or missing (`FALSE`) for each subject.

- spline.basis:

  The evaluation matrix of the B-spline basis functions at the time
  points.

- intercept:

  Boolean. If TRUE (default) the design starts with an intercept block
  built from `spline.basis[, -1]`; if FALSE only the score-scaled spline
  blocks are produced.

## Value

A matrix containing the stacked tensor products of the (optional)
intercept bases and the scaled spline bases for all valid observations.
