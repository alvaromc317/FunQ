# Cross validation of degrees of freedom

Performs cross validation on degrees of freedom parameter of fqpca

## Usage

``` r
fqpca_cv_df(
  data,
  colname = NULL,
  npc = 2,
  pve = NULL,
  quantile.value = 0.5,
  lambda.ridge = 0,
  n.folds = 3,
  return.models = TRUE,
  criteria = "points",
  periodic = TRUE,
  splines.df.grid = c(5, 10, 15, 20),
  tol = 0.001,
  max.iters = 20,
  splines.method = "conquer",
  penalized = FALSE,
  verbose.fqpca = FALSE,
  verbose.cv = TRUE,
  seed = NULL
)
```

## Arguments

- data:

  An \\(N \times T)\\ matrix, a tf object from the tidyfun package or a
  data.frame containing the functional data as a tf column.

- colname:

  The name of the column containing the functional data. Use only if
  data is a dataframe and colname is a column in the dataframe.

- npc:

  The number of estimated components.

- pve:

  Float between 0 and 1. Percentage of variability explained by
  components. This affects the number of components used in the curve
  reconstruction and error estimation. Set to NULL to avoid this
  behavior.

- quantile.value:

  The quantile considered.

- lambda.ridge:

  Hyper parameter controlling the penalization on the second derivative
  of the splines. It has effect only with `penalized=TRUE` and
  `method='conquer'`.

- n.folds:

  Number of folds to be used on cross validation.

- return.models:

  Should the list of all the models built be returned?

- criteria:

  Criteria used to divide the data. Valid values are `'rows'`, which
  considers the division based on full rows, or `'points'`, which
  considers the division based on points within the matrix.

- periodic:

  Boolean indicating if the data is expected to be periodic (start
  coincides with end) or not.

- splines.df.grid:

  Grid of possible values for the degrees of freedom.

- tol:

  Tolerance on the convergence of the algorithm.

- max.iters:

  Maximum number of iterations.

- splines.method:

  Method used in the resolution of the splines quantile regression
  model. It currently accepts the methods `c('conquer', 'quantreg')`.

- penalized:

  Boolean indicating if the smoothness should be controlled using a
  second derivative penalty. This functionality is experimental.

- verbose.fqpca:

  Boolean indicating verbosity of the fqpca function.

- verbose.cv:

  Boolean indicating verbosity of the cross-validation process.

- seed:

  Seed for the random generator number.

## Value

A list containing the matrix of scores, the matrix of loadings, and a
secondary list with extra information.

## Examples

``` r
n.obs = 150
n.time = 144

# Generate scores
c1.vals = rnorm(n.obs)
c2.vals = rnorm(n.obs)

# Generate pc's
pc1 = sin(seq(0, 2*pi, length.out = 144))
pc2 = cos(seq(0, 2*pi, length.out = 144))

# Generate data
Y <- c1.vals * matrix(pc1, nrow = n.obs, ncol=n.time, byrow = TRUE) +
c2.vals * matrix(pc2, nrow = n.obs, ncol=n.time, byrow = TRUE)

# Add noise
Y <- Y + matrix(rnorm(n.obs * n.time, 0, 0.4), nrow = n.obs)

# Add missing observations
Y[sample(n.obs*n.time, as.integer(0.2*n.obs*n.time))] <- NA

cv_result <- fqpca_cv_df(data=Y, splines.df.grid = c(5, 10, 15), n.folds = 2)
#> Degrees of freedom: 5 ---------------------
#> 2026-07-16 19:48:47. Fold: 1
#> 2026-07-16 19:48:47. Fold: 2
#> Degrees of freedom: 5. Execution completed in: 0.397 seconds.
#> Degrees of freedom: 10 ---------------------
#> 2026-07-16 19:48:47. Fold: 1
#> 2026-07-16 19:48:48. Fold: 2
#> Degrees of freedom: 10. Execution completed in: 0.383 seconds.
#> Degrees of freedom: 15 ---------------------
#> 2026-07-16 19:48:48. Fold: 1
#> 2026-07-16 19:48:48. Fold: 2
#> Degrees of freedom: 15. Execution completed in: 0.805 seconds.
```
