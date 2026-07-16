# Split points of a given matrix Y into train / test

Split points of a given matrix Y into train / test based on parameters
train.pct and train.size. If both are informed train.pct takes
preference. This function keeps the dimensions of the original Y in both
train and test and substitutes the values in both splits by NAs

## Usage

``` r
train_test_split_points(Y, train.pct = NULL, train.size = NULL, seed = NULL)
```

## Arguments

- Y:

  \\(N \times P)\\ matrix of predictive variables.

- train.pct:

  Float number indicating the % of rows used for training.

- train.size:

  Integer number indicating number of rows used for training.
  `train.size` is superseded by `train.pct`.

- seed:

  Seed for the random generator number.

## Value

A list containing a matrix Y.train and a matrix Y.test
