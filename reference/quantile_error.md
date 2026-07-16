# Quantile error computation

Quantile error metric computation

## Usage

``` r
quantile_error(Y, Y.pred, quantile.value)
```

## Arguments

- Y:

  \\(N \times T)\\ matrix of observed time instants.

- Y.pred:

  \\(N \times T)\\ matrix of predicted time instants.

- quantile.value:

  The quantile considered.

## Value

The quantile error between the two matrices
