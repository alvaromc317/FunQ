# MFQPCA (Multilevel Functional Quantile Principal Component Analysis)

Solves the multilevel functional quantile principal component analysis
methodology

## Usage

``` r
mfqpca(
  data,
  group,
  colname = NULL,
  npc.between = 2,
  npc.within = 2,
  quantile.value = 0.5,
  periodic = TRUE,
  splines.df = 10,
  splines.method = "conquer",
  tol = 0.01,
  max.iters = 10,
  verbose = FALSE,
  seed = NULL
)
```

## Arguments

- data:

  An \\(N \times T)\\ matrix, a tf object from the tidyfun package or a
  data.frame containing the functional data as a tf column.

- group:

  Either a string or an array. If it is a string, it must point to the
  grouping variable in data only if data is a dataframe. If an array, it
  must be the N dimensional array indicating the hierarchical structure
  of the data. Elements in the array with the same value indicate they
  are repeated measures of the same individual.

- colname:

  The name of the column containing the functional data. Use only if
  data is a dataframe and colname is a column in the dataframe.

- npc.between:

  The number of estimated between level components.

- npc.within:

  The number of estimated within level components.

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

- tol:

  Tolerance on the convergence of the algorithm.

- max.iters:

  Maximum number of iterations.

- verbose:

  Boolean indicating the verbosity.

- seed:

  Seed for the random generator number.

## Value

mfqpca_object

## Examples

``` r

n.individuals <- 20
n.repeated <- 10
n.time = 144
N <- n.repeated * n.individuals

group <- rep(1:n.individuals, each=n.repeated)

# Define score values using a normal distribution
c1.vals <- rnorm(n.individuals)
c1.vals <- c1.vals[match(group, unique(group))]
c2.vals <- rnorm(N)

# Define principal components
pcb <- sin(seq(0, 2*pi, length.out = n.time))
pcw <- cos(seq(0, 2*pi, length.out = n.time))

# Generate a data matrix and add missing observations
Y <- c1.vals * matrix(pcb, nrow = N, ncol=n.time, byrow = TRUE) +
c2.vals * matrix(pcw, nrow = N, ncol=n.time, byrow = TRUE)

# Add noise
Y <- Y + matrix(rnorm(N*n.time, 0, 0.4), nrow = N)
Y[sample(N*n.time, as.integer(0.2*N))] <- NA

results <- mfqpca(data = Y, group = group, npc.between = 1, npc.within=1, quantile.value = 0.5)
```
