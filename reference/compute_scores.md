# Compute scores

Inner function to compute the per-subject scores via independent
quantile regressions on the loading functions.

## Usage

``` r
compute_scores(Y, Y.mask, loadings, quantile.value, offset = NULL)
```

## Arguments

- Y:

  The \\(N \times T)\\ matrix of observed time instants.

- Y.mask:

  Mask matrix of the same dimensions as Y indicating which observations
  in Y are known.

- loadings:

  Matrix of loading coefficients.

- quantile.value:

  The quantile considered.

- offset:

  Optional offset subtracted from Y before fitting: NULL (no offset), a
  length-T vector (e.g. a population intercept), or an \\(N \times T)\\
  matrix (e.g. fitted values from a regression part).

## Value

The matrix of scores.
