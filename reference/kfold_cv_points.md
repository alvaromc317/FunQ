# Split a given matrix Y into k-folds

Splits a given matrix Y into k-folds. This function keeps the dimensions
of the original Y in both train and test and substitutes the values in
both split by NAs

## Usage

``` r
kfold_cv_points(Y, folds = 3, seed = NULL)
```

## Arguments

- Y:

  \\(N \times P)\\ matrix of predictive variables.

- folds:

  Integer number indicating number of folds.

- seed:

  Seed for the random generator number.

## Value

A list containing two inside lists, one for training and one for
testing. The length of the inside lists is equal to the number of folds
