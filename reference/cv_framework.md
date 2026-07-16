# Generic cross-validation framework

Inner engine shared by all the cross-validation drivers: splits the data
into folds, loops over a hyper-parameter grid fitting a model per (grid
value, fold) pair via `fit_fun`, reconstructs the test predictions via
`reconstruct_fun`, and collects the quantile error matrix.

## Usage

``` r
cv_framework(
  Y,
  grid,
  n.folds,
  criteria,
  quantile.value,
  fit_fun,
  reconstruct_fun,
  return.models = TRUE,
  verbose.cv = TRUE,
  seed = NULL,
  grid.label = "Degrees of freedom: ",
  model.prefix = "df_idx"
)
```

## Arguments

- Y:

  The \\(N \times T)\\ matrix of functional data.

- grid:

  Numeric vector of hyper-parameter values to cross-validate.

- n.folds:

  Number of folds to be used on cross validation.

- criteria:

  Criteria used to divide the data. Valid values are `'rows'` or
  `'points'`.

- quantile.value:

  The quantile considered.

- fit_fun:

  Function `(Y.train, grid.value, train.idx)` returning a fitted model.
  `train.idx` contains the training row indices when `criteria='rows'`
  and is NULL otherwise.

- reconstruct_fun:

  Function `(model, Y.test, train.idx)` returning the predicted \\(N
  \times T)\\ matrix used to score the fold.

- return.models:

  Should the list of all the models built be returned?

- verbose.cv:

  Boolean indicating verbosity of the cross-validation process.

- seed:

  Seed for the random generator number (used to build the folds).

- grid.label:

  Label prefixing the grid value in verbose messages.

- model.prefix:

  Prefix used to name the stored models (`'<prefix>=<i>.fold=<j>'`).

## Value

A list containing the error matrix, the execution time, the grid, the
criteria and the list of models.
