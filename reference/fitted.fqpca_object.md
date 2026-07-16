# Fit Yhat

S3 method for class 'fqpca_object'. Given an fqpca_object model,
estimates Yhat for different pve values.

## Usage

``` r
# S3 method for class 'fqpca_object'
fitted(object, pve = 0.95, ...)
```

## Arguments

- object:

  An object output of the fqpca function.

- pve:

  Percentage of explained variability (between 0 and 1) used to select
  the number of components in Yhat estimation. Set to NULL to use all
  components.

- ...:

  further arguments passed to or from other methods.

## Value

The matrix of fitted values.

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

Yhat <- fitted(object = results, pve=0.99)
```
