# Split a given matrix Y into k-folds

Splits full rows of a given matrix Y into k-folds

## Usage

``` r
kfold_cv_rows(Y, folds = 3, seed = NULL)
```

## Arguments

- Y:

  \\(N \times P)\\ matrix of predictive variables.

- folds:

  Integer number indicating number of folds.

- seed:

  Seed for the random generator number.

## Value

A list containing three inside lists, one for training, one for testing
and one with training indices.. The length of the inside lists is equal
to the number of folds
