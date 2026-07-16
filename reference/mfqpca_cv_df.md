# Cross validation of degrees of freedom

Performs cross validation on degrees of freedom parameter of mfqpca

## Usage

``` r
mfqpca_cv_df(
  data,
  group,
  colname = NULL,
  npc.between = 2,
  npc.within = 2,
  pve.between = NULL,
  pve.within = NULL,
  quantile.value = 0.5,
  n.folds = 3,
  return.models = TRUE,
  criteria = "points",
  periodic = TRUE,
  splines.df.grid = c(5, 10, 15, 20),
  tol = 0.001,
  max.iters = 20,
  splines.method = "conquer",
  verbose.mfqpca = FALSE,
  verbose.cv = TRUE,
  seed = NULL
)
```

## Arguments

- data:

  An \\(N \times T)\\ matrix, a tf object from the tidyfun package or a
  data.frame containing the functional data as a tf column.

- group:

  Either a string or an array. If it is a string, it must point to the
  grouping variable in data only if data is a dataframe. If an array, it
  must be the N dimensional array indicating the hierarchical structure
  of the data. Elements in the array with the same value indicate they
  are repeated measures of the same individual.

- colname:

  The name of the column containing the functional data. Use only if
  data is a dataframe and colname is a column in the dataframe.

- npc.between:

  The number of estimated between level components.

- npc.within:

  The number of estimated within level components.

- pve.between:

  Float between 0 and 1. Percentage of variability explained by between
  level components. This affects the number of components used in the
  curve reconstruction and error estimation. Set to NULL to avoid this
  behavior.

- pve.within:

  Float between 0 and 1. Percentage of variability explained by within
  level components. This affects the number of components used in the
  curve reconstruction and error estimation. Set to NULL to avoid this
  behavior.

- quantile.value:

  The quantile considered.

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

- verbose.mfqpca:

  Boolean indicating verbosity of the fqpca function.

- verbose.cv:

  Boolean indicating verbosity of the cross-validation process.

- seed:

  Seed for the random generator number.

## Value

A list containing the matrices of scores, the matrices of loadings, and
a secondary list with extra information.

## Examples

``` r
n.individuals <- 20
n.repeated <- 10
n.time = 144
N <- n.repeated * n.individuals

group <- rep(1:n.individuals, each=n.repeated)

# Define score values using a normal distribution
c1.vals <- rnorm(n.individuals)
c1.vals <- c1.vals[match(group, unique(group))]
c2.vals <- rnorm(N)

# Define principal components
pcb <- sin(seq(0, 2*pi, length.out = n.time))
pcw <- cos(seq(0, 2*pi, length.out = n.time))

# Generate a data matrix and add missing observations
Y <- c1.vals * matrix(pcb, nrow = N, ncol=n.time, byrow = TRUE) +
c2.vals * matrix(pcw, nrow = N, ncol=n.time, byrow = TRUE)
Y <- Y + matrix(rnorm(N*n.time, 0, 0.4), nrow = N)
Y[sample(N*n.time, as.integer(0.2*N))] <- NA

cv_results <- mfqpca_cv_df(data = Y, group = group,  splines.df.grid = c(5, 10, 15), n.folds = 2)
#> Degrees of freedom: 5 ---------------------
#> 2026-07-16 19:48:55. Fold: 1
#> 2026-07-16 19:48:56. Fold: 2
#> Degrees of freedom: 5. Execution completed in: 0.79 seconds.
#> Degrees of freedom: 10 ---------------------
#> 2026-07-16 19:48:56. Fold: 1
#> 2026-07-16 19:48:57. Fold: 2
#> Degrees of freedom: 10. Execution completed in: 1.191 seconds.
#> Degrees of freedom: 15 ---------------------
#> 2026-07-16 19:48:57. Fold: 1
#> 2026-07-16 19:48:58. Fold: 2
#> Degrees of freedom: 15. Execution completed in: 1.412 seconds.
```
