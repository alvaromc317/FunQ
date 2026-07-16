# Rotation of loadings and scores

Performs the SVD-based rotation of loadings and scores (with sign
alignment) in order to ensure the solution is unique.

## Usage

``` r
rotate_factors(scores, loadings, intercept = NULL)
```

## Arguments

- scores:

  Matrix of scores.

- loadings:

  Matrix of loadings.

- intercept:

  Optional population intercept. When provided, the scores are
  mean-centered and the effect of their mean is moved into the
  intercept; when NULL (e.g. the mfqpca within level) the scores are
  rotated as they are.

## Value

The rotated matrices of loadings and scores, the (possibly updated)
intercept and the rotation matrix.
