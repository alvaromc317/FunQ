# Check fqpca input parameters

Check input parameters of the fqpca methodology (also used by its CV
drivers).

## Usage

``` r
check_fqpca_params(
  npc,
  quantile.value,
  periodic,
  splines.df,
  splines.method,
  penalized,
  lambda.ridge,
  tol,
  max.iters,
  verbose,
  seed
)
```

## Arguments

- npc:

  The number of estimated components.

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

- penalized:

  Boolean indicating if the smoothness should be controlled using a
  second derivative penalty. This functionality is experimental and is
  much slower than the control of the smoothness using the degrees of
  freedom.

- lambda.ridge:

  Hyper parameter controlling the penalization on the second derivative
  of the splines. It has effect only with `penalized=TRUE`.

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
