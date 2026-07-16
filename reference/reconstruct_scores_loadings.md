# Truncated reconstruction from scores and loadings

Helper function to compute fitted values from the first `n.comp`
components.

## Usage

``` r
reconstruct_scores_loadings(scores, loadings, n.comp)
```

## Arguments

- scores:

  matrix of scores

- loadings:

  matrix of loadings

- n.comp:

  number of components used. If 0 then a matrix of 0s is returned.
