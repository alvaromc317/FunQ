# Split a given matrix Y into train / test

Splits a given matrix into a train / test split based on two possible
criteria: either based on the total number of rows or the total number
of data points.

## Usage

``` r
train_test_split(
  Y,
  criteria = "points",
  train.pct = NULL,
  train.size = NULL,
  seed = NULL
)
```

## Arguments

- Y:

  \\(N \times P)\\ matrix of predictive variables.

- criteria:

  Criteria used to divide the data. Valid values are `'rows'`, which
  considers the division based on full rows, or `'points'`, which
  considers the division based on points within the matrix.

- train.pct:

  Float number indicating the % of rows used for training. This takes
  precedence over `train.size`.

- train.size:

  Integer number indicating number of rows used for training.
  `train.size` is superseded by `train.pct`.

- seed:

  Seed for the random generator number.

## Value

A list containing a matrix Y.train and a matrix Y.test

## Examples

``` r
# Generate a small matrix
Y <- matrix(rnorm(50), nrow = 10)

# Divide based on full rows
train_test_rows <- train_test_split(Y, criteria = 'rows', train.pct = 0.5)
Y.train <- train_test_rows$Y.train
Y.test <- train_test_rows$Y.test

train_test_points <- train_test_split(Y, criteria = 'points', train.pct = 0.5)
Y.train <- train_test_points$Y.train
Y.test <- train_test_points$Y.test
```
