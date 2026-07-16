# FQPCA (Functional Quantile Principal Component Analysis)

Solves the functional quantile principal component analysis methodology

## Usage

``` r
fqpca(
  data,
  colname = NULL,
  npc = 2,
  quantile.value = 0.5,
  periodic = TRUE,
  splines.df = 10,
  splines.method = "conquer",
  penalized = FALSE,
  lambda.ridge = 0,
  tol = 0.001,
  max.iters = 20,
  verbose = FALSE,
  seed = NULL
)
```

## Arguments

- data:

  An \\(N \times T)\\ matrix, a tf object from the tidyfun package or a
  data.frame containing the functional data as a tf column.

- colname:

  The name of the column containing the functional data. Use only if
  data is a dataframe and colname is a column in the dataframe.

- npc:

  The number of estimated components.

- quantile.value:

  The quantile considered.

- periodic:

  Boolean indicating if the data is expected to be periodic (start
  coincides with end) or not.

- splines.df:

  Degrees of freedom for the splines.

- splines.method:

  Method used in the resolution of the splines quantile regression
  model. It currently accepts the methods `c('conquer', 'quantreg')`.

- penalized:

  Boolean indicating if the smoothness should be controlled using a
  second derivative penalty. This functionality is experimental.

- lambda.ridge:

  Hyper parameter controlling the penalization on the second derivative
  of the splines. It has effect only with `penalized=TRUE` and
  `method='conquer'`.

- tol:

  Tolerance on the convergence of the algorithm.

- max.iters:

  Maximum number of iterations.

- verbose:

  Boolean indicating the verbosity.

- seed:

  Seed for the random generator number.

## Value

fqpca_object

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

intercept <- results$intercept
loadings <- results$loadings
scores <- results$scores
```
