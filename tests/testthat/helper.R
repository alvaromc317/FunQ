
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
