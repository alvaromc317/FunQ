# Split a given matrix Y into k-folds

Splits a given matrix Y into k-folds. based on two possible criteria:
either based on the total number of rows or the total number of data
points.

## Usage

``` r
create_folds(Y, criteria = "points", folds = 3, seed = NULL)
```

## Arguments

- Y:

  \\(N \times P)\\ matrix of predictive variables.

- criteria:

  Criteria used to divide the data. Valid values are `'rows'`, which
  considers the division based on full rows, or `'points'`, which
  considers the division based on points within the matrix.

- folds:

  Integer number indicating number of folds

- seed:

  Seed for the random generator number.

## Value

A list containing two inside lists, one for training and one for
testing. The length of the inside lists is equal to the number of folds

## Examples

``` r
# Generate a small matrix
Y <- matrix(rnorm(50), nrow = 10)

# Divide based on full rows
kfolds_rows <- create_folds(Y, criteria = 'rows', folds = 3, seed = 1)
Y.train.list <- kfolds_rows$Y.train.list
Y.test.list <- kfolds_rows$Y.test.list

kfolds_points <- create_folds(Y, criteria='points', folds=3, seed=1)
Y.train.list <- kfolds_points$Y.train.list
Y.test.list <- kfolds_points$Y.test.list
```
