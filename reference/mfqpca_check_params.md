# Check mfqpca input parameters

Check input parameters of the mfqpca methodology (also used by its CV
driver).

## Usage

``` r
mfqpca_check_params(
  npc.between,
  npc.within,
  quantile.value,
  periodic,
  splines.df,
  splines.method,
  tol,
  max.iters,
  verbose,
  seed
)
```

## Arguments

- npc.between:

  The number of estimated between components.

- npc.within:

  The number of estimated within components.

- quantile.value:

  The quantile considered.

- periodic:

  Boolean indicating if the data is expected to be periodic (start
  coincides with end) or not.

- splines.df:

  Degrees of freedom for the splines.

- splines.method:

  Method used in the resolution of the splines quantile regression
  model. It currently accepts the methods `c('conquer', 'quantreg')`.

- tol:

  Tolerance on the convergence of the algorithm.

- max.iters:

  Maximum number of iterations.

- verbose:

  Boolean indicating the verbosity.

- seed:

  Seed for the random generator number.

## Value

No return
