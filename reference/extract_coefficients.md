# Extract the coefficient vector from a quantile regression model object

Returns the flat coefficient vector regardless of whether the model was
fitted with
[`conquer::conquer`](https://rdrr.io/pkg/conquer/man/conquer.html)
(`$coeff`) or [`quantreg::rq`](https://rdrr.io/pkg/quantreg/man/rq.html)
(`$coefficients`).

## Usage

``` r
extract_coefficients(model)
```

## Arguments

- model:

  A fitted model object.

## Value

The coefficient vector.
