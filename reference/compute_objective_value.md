# Compute objective value function.

Inner function to compute the objective value of the alternating
quantile algorithms at each iteration.

## Usage

``` r
compute_objective_value(Y, quantile.value, scores, intercept, loadings)
```

## Arguments

- Y:

  \\(N \times T)\\ matrix of observed time instants.

- quantile.value:

  The quantile considered.

- scores:

  The matrix of estimated scores.

- intercept:

  population intercept

- loadings:

  The matrix of estimated loadings.

## Value

The objective value function.
