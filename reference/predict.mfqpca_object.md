# Predict mfqpca scores

S3 method for class 'mfqpca_object' Given a new matrix Y, predicts the
value of the scores associated to the given matrix.

## Usage

``` r
# S3 method for class 'mfqpca_object'
predict(object, newdata, newdata.group, ...)
```

## Arguments

- object:

  An object output of the fqpca function.

- newdata:

  The N by T matrix of observed time instants to be tested

- newdata.group:

  An N dimensional array indicating the hierarchical structure of the
  data. Elements in the array with the same value indicate they are
  repeated measures of the same individual.

- ...:

  further arguments passed to or from other methods.

## Value

List of matrices of predicted scores.

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
Y <- Y + matrix(rnorm(N*n.time, 0, 0.4), nrow = N)
Y[sample(N*n.time, as.integer(0.2*N))] <- NA

results <- mfqpca(data = Y, group=group, npc.between = 1, npc.within=1, quantile.value = 0.5)
predictions <- predict(object = results, newdata = Y[101:150,], newdata.group = group[101:150])
```
