# QUANTILE FACTOR ANALYSIS ALGORITHM ------------------------------------------

#' QFA
#'
#' Solves the Quantile factor analysis methodology.
#'
#' @param x The N by T matrix of observed time instants
#' @param n.components The number of estimated components. Default is 2.
#' @param quantile.value The quantile considered. Default is the median 0.5.
#' @param tol Tolerance on the convergence of the algorithm. Smaller values can spped up computation but may affect the quality of the estimations. Default is 1e-3.
#' @param n.iters Maximum number of iterations. Default is 30.
#' @param verbose Boolean indicating verbosity of the function. Default is FALSE.
#' @param seed Seed for the random generator number. Parameter included for reproducibility purposes. Default is NULL (meaning no seed is assigned).
#'
#' @return A list containing the matrix of scores, the matrix of loadings, and a secondary list with extra information.
#' @export
#'
#' @examples
#' # Generate fake dataset with 200 observations and 144 time points
#'
#' x = matrix(rep(sin(seq(0, 2*pi, length.out=144)), 200), byrow=TRUE, nrow=200)
#' x = x + matrix(rnorm(200*144, 0, 0.4), nrow=200)
#'
#' results = qfa(x=x, n.components=1, quantile.value=0.5)
#'
#' loadings = results$loadings
#' scores = results$scores
qfa = function(x, n.components, quantile.value=0.5, tol=1e-3, n.iters=100, verbose=TRUE, seed=NULL)
{
  # param X: is the N by T matrix of observed variables
  # param n.components: is the number of estimated factors. Default is 2
  # param quantile.value: quantile level. Defaul is 0.5
  # param tol: convergence threshold.
  # param n.iters: max number of iterations
  # param verbose: boolean, should execution time per iteration be displayed?
  # param seed: seed for the random generator number

  global_start_time = base::Sys.time()
  if(!is.null(seed)){base::set.seed(seed)}
  x = base::as.matrix(x)
  # Compute the IQR algorithm defined in paper from Chen2020
  # Step 1: Initialize values
  n_obs = base::dim(x)[1]
  n_time = base::dim(x)[2]
  F0 = base::matrix(c(stats::rnorm(n_time*n.components)), nrow = n_time)
  Lambda0 = c()
  for(j in 1:n_obs){ Lambda0 = base::rbind(Lambda0, quantreg::rq(x[j,] ~ -1 + F0, tau=quantile.value)$coefficients) }
  # Given Lambda0 and F0 compute of convergence criteria
  of_value0 = base::mean((quantile.value - (x - Lambda0 %*% t(F0) < 0)) * (x - Lambda0 %*% t(F0)))
  convergence_criteria = c()
  of_value = c(of_value0)
  for(i in 1:n.iters)
  {
    loop_start_time = base::Sys.time()
    F1 = c()
    for(k in 1:n_time){ F1 = base::rbind(F1, quantreg::rq(x[,k] ~ -1 + Lambda0, tau=quantile.value)$coefficients)}
    # Given F1 compute Lambda1
    Lambda1 = c()
    for(j in 1:n_obs){ Lambda1 = base::rbind(Lambda1, quantreg::rq(x[j,] ~ -1 + F1, tau=quantile.value)$coefficients) }
    # Given Lambda1 and F1 compute of convergence criteria
    of_value1 = base::mean((quantile.value - (x - Lambda1 %*% t(F1) < 0)) * (x - Lambda1 %*% t(F1)))
    of_value = c(of_value, of_value1)
    convergence_criteria[i] = base::abs(of_value1 - of_value0)
    loop_end_time = base::Sys.time()
    loop_execution_time = loop_end_time - loop_start_time
    if(verbose)
    {
      message('Iteration ', i, ' completed in ', base::round(loop_execution_time, 3), ' seconds',
              '\n', 'Convergence criteria value: ', base::round(convergence_criteria[i], 4),
              '\n', 'OF0 value: ', base::round(of_value0, 4), ' OF1 value: ', round(of_value1, 4))
      message('___________________________________________________')
    }
    if(convergence_criteria[i] < tol)
    {
      if(verbose){message('Convergence criteria is satisfied. Value: ', base::round(convergence_criteria[i], 4))}
      break
    }
    Lambda0 = Lambda1
    of_value0 = of_value1
    # Measure computation time
  }
  if(i == n.iters){message('Algorithm did not converge in ', n.iters, ' iterations')}
  colnames(F1) = NULL
  Fhat = F1
  Lhat = Lambda1
  # Step 4: Normalize variables (code extracted from IQR.m file (original))
  sigmaF = (1/n_time) * t(Fhat) %*% Fhat
  sigmaA = (1/n_obs) * t(Lhat) %*% Lhat
  dum1 = expm::sqrtm(sigmaF) %*% sigmaA %*% expm::sqrtm(sigmaF)
  svd_decomp = base::svd(dum1)
  R = expm::sqrtm(base::solve(sigmaF)) %*% svd_decomp$u
  Fhat = Fhat %*% R;
  Lhat = Lhat %*% t(base::solve(R))
  global_end_time = base::Sys.time()
  global_execution_time = global_end_time - global_start_time
  results = list(
    loadings = Fhat,
    scores = Lhat,
    extra_information = list(
      unnormalized_loadings = F1,
      normalization_matrix = R,
      quantile.value = quantile.value,
      objective_function_value = of_value1,
      list_objective_function_values = of_value,
      n_iter=i,
      execution_time=global_execution_time
    )
  )
  return(results)
}

# PREDICTIONS -----------------------------------------------------------------

#' Predict QFA
#' @param qfa_object An object output of the qfa function.
#' @param x The N by T matrix of observed time instants to be tested.
#'
#' @return The normalized matrix of scores.
#' @export
#'
#' @examples
#' # Generate fake dataset with 200 observations and 144 time points
#'
#' x = matrix(rep(sin(seq(0, 2*pi, length.out=144)), 200), byrow=TRUE, nrow=200)
#' x = x + matrix(rnorm(200*144, 0, 0.4), nrow=200)
#'
#' results = qfa(x=x[1:150,], n.components=1, quantile.value=0.5)
#'
#' predictions = predict_scores_qfa(results, x=x[151:200,])
predict_scores_qfa = function(qfa_object, x)
{
  # Unload required information from fqfa object
  normalization_matrix=qfa_object$extra_information$normalization_matrix
  quantile.value = qfa_object$extra_information$quantile.value
  unnormalized_loadings = qfa_object$extra_information$unnormalized_loadings
  x = as.matrix(x)
  n_obs = nrow(x)
  n.components = ncol(unnormalized_loadings)
  Lambda1 =  matrix(0, n_obs, n.components)
  for(i in seq(n_obs))
  {
    xi = x[i, ]
    Lambda1[i,] =  quantreg::rq(xi ~ -1 + unnormalized_loadings, tau=quantile.value)$coefficients
  }
  Lambda1 = Lambda1 %*% t(base::solve(normalization_matrix))
  return(Lambda1)
}
