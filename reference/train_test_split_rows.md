# Split rows of a given matrix Y into train / test

Split rows of a given matrix Y into train / test based on parameters
train.pct and train.size. If both are informed train.pct takes
preference

## Usage

``` r
train_test_split_rows(Y, train.pct = NULL, train.size = NULL, seed = NULL)
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

A list containing a matrix Y.train, a matrix Y.test and a list of
training indices.
