# Compute scores between

Inner function to compute the between level scores of the mfqpca
methodology.

## Usage

``` r
mfqpca_compute_scores_between(
  Y,
  Y.mask,
  group,
  intercept,
  loadings,
  quantile.value
)
```

## Arguments

- Y:

  The \\(N \times T)\\ matrix of observed time instants.

- Y.mask:

  Mask matrix of the same dimensions as Y indicating which observations
  in Y are known.

- group:

  An N dimensional array indicating the hierarchical structure of the
  data. Elements in the array with the same value indicate they are
  repeated measures of the same individual.

- intercept:

  A T dimensional vector of intercept values.

- loadings:

  Matrix of loading coefficients.

- quantile.value:

  The quantile considered.

## Value

The matrix of between level scores.
