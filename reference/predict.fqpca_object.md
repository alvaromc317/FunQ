# Predict fqpca scores

S3 method for class 'fqpca_object' Given a new matrix Y, predicts the
value of the scores associated to the given matrix.

## Usage

``` r
# S3 method for class 'fqpca_object'
predict(object, newdata, ...)
```

## Arguments

- object:

  An object output of the fqpca function.

- newdata:

  The N by T matrix of observed time instants to be tested, or the
  dataframe storing the tf functional vector using the same colname as
  the one used in the fqpca function.

- ...:

  further arguments passed to or from other methods.

## Value

The predicted matrix of scores.

## Examples

``` r

n.obs = 150
n.time = 144

# Generate scores
c1.vals = rnorm(n.obs)
c2.vals = rnorm(n.obs)

# Generate pc's
pc1 = sin(seq(0, 2*pi, length.out = n.time))
pc2 = cos(seq(0, 2*pi, length.out = n.time))

# Generate data
Y <- c1.vals * matrix(pc1, nrow = n.obs, ncol=n.time, byrow = TRUE) +
c2.vals * matrix(pc2, nrow = n.obs, ncol=n.time, byrow = TRUE)

# Add noise
Y <- Y + matrix(rnorm(n.obs * n.time, 0, 0.4), nrow = n.obs)

# Add missing observations
Y[sample(n.obs*n.time, as.integer(0.2*n.obs*n.time))] <- NA

results <- fqpca(data = Y, npc = 2, quantile.value = 0.5, seed=1)

predictions <- predict(object = results, newdata = Y[101:150,])
```
