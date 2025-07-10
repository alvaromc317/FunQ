
test_data_fqpca = function(n_obs=20, n_time=12, proportion_na=0.2, seed=5)
{
  if(!is.null(seed)){set.seed(seed)}
  Y = matrix(rep(sin(seq(0, 2*pi, length.out=n_time)), n_obs), byrow=TRUE, nrow=n_obs)
  Y = Y + matrix(rnorm(n_obs*n_time, mean=0, sd=0.4), nrow=n_obs)
  Y[sample(n_obs*n_time, as.integer(proportion_na*n_obs*n_time))] = NA
  return(Y)
}

test_data_cvxr = function(n_obs=20, n_features=2, seed=5)
{
  if(!is.null(seed)){set.seed(seed)}
  x = matrix(rnorm(n_obs*n_features), nrow=n_obs)
  y = x %*% matrix(1:n_features, nrow=n_features) +  rnorm(n_obs)
  return(list('x'=x, 'y'=y))
}

test_data_mfqpca = function(proportion_na=0.2, seed=5)
{
  if(!is.null(seed)){set.seed(seed)}
  n.individuals <- 10
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
  Y <- Y + matrix(rnorm(N*n.time, 0, proportion_na), nrow = N)
  Y[sample(N*n.time, as.integer(0.2*N))] <- NA
  return(list(Y, group))
}
