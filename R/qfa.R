# QUANTILE FACTOR ANALYSIS ALGORITHM ------------------------------------------

new_qfa <- function(loadings, scores, unnormalized_loadings, normalization_matrix,
                    objective_function_value, list_objective_function_values,
                    quantile_value, n_components=n_components, n_iters, execution_time)
{
  structure(list(
    loadings = loadings,
    scores = scores,
    unnormalized_loadings = unnormalized_loadings,
    normalization_matrix = normalization_matrix,
    objective_function_value = objective_function_value,
    list_objective_function_values = list_objective_function_values,
    quantile_value = quantile_value,
    n_components = n_components,
    n_iters = n_iters,
    execution_time = execution_time
  ),
  class = "qfa")
}


#' QFA
#'
#' Solves the Quantile factor analysis methodology.
#'
#' @param x The N by T matrix of observed time instants
#' @param n_components The number of estimated components. Default is 2.
#' @param quantile_value The quantile considered. Default is the median 0.5.
#' @param tol Tolerance on the convergence of the algorithm. Smaller values can spped up computation but may affect the quality of the estimations. Default is 1e-3.
#' @param n_iters Maximum number of iterations. Default is 30.
#' @param verbose Boolean indicating verbosity of the function. Default is FALSE.
#' @param seed Seed for the random generator number. Parameter included for reproducibility purposes. Default is NULL (meaning no seed is assigned).
#' @param method Parameter from the quantreg package. Method used to solve the quantile regression models.
#' @return A list containing the matrix of scores, the matrix of loadings, and a secondary list with extra information.
#' @export
#'
#' @examples
#' # Generate fake dataset with 200 observations and 144 time points
#'
#' x = matrix(rep(sin(seq(0, 2*pi, length.out=144)), 200), byrow=TRUE, nrow=200)
#' x = x + matrix(rnorm(200*144, 0, 0.4), nrow=200)
#'
#' results = qfa(x=x, n_components=1, quantile_value=0.5)
#'
#' loadings = results$loadings
#' scores = results$scores
qfa = function(x, n_components, quantile_value=0.5, tol=1e-3, n_iters=100, verbose=FALSE, seed=NULL, method='br')
{
  if(!is.data.frame(x) & !is.matrix(x)){stop('x is not of a valid type. Object provided: ', typeof(x))}
  if(! n_components == floor(n_components)){stop('n_components must be an integer number. Value provided: ', n_components)}
  if(quantile_value<0 | quantile_value>1){stop('quantile_value must be a value between 0 and 1. Value provided: ', quantile_value)}
  if(! n_iters == floor(n_iters)){stop('n_iters must be an integer number. Value provided: ', n_iters)}

  global_start_time = base::Sys.time()
  if(!is.null(seed)){base::set.seed(seed)}
  x = unname(base::as.matrix(x))
  if(sum(is.na(x)) > 0){stop('x contains missing data. QFA algorithm is unable to deal with these. Use FQPCA algorithm instead')}

  # Step 1: Initialize values
  n_obs = base::dim(x)[1]
  n_time = base::dim(x)[2]
  loadings0 = base::matrix(c(stats::rnorm(n_time*n_components)), nrow = n_time)
  scores0 = base::matrix(0, nrow=n_obs, ncol=n_components)
  for(j in 1:n_obs){ scores0[j,] = quantreg::rq(x[j,] ~ -1 + loadings0, tau=quantile_value, method=method)$coefficients }
  # Given scores0 and loadings0 compute convergence criteria value
  of_value0 = base::mean((quantile_value - (x - scores0 %*% t(loadings0) < 0)) * (x - scores0 %*% t(loadings0)))
  convergence_criteria = c()
  of_value = c(of_value0)

  for(i in 1:n_iters)
  {
    # Given scores from previous step, compute loadings
    loop_start_time = base::Sys.time()
    loadings = base::matrix(0, nrow=n_time, ncol=n_components)
    for(k in 1:n_time){ loadings[k,] = quantreg::rq(x[,k] ~ -1 + scores0, tau=quantile_value, method=method)$coefficients}
    # Given loadings compute scores
    scores = base::matrix(0, nrow=n_obs, ncol=n_components)
    for(j in 1:n_obs){ scores[j,] = quantreg::rq(x[j,] ~ -1 + loadings, tau=quantile_value, method=method)$coefficients }
    # Given scores and loadings compute convergence criteria
    of_value1 = base::mean((quantile_value - (x - scores %*% t(loadings) < 0)) * (x - scores %*% t(loadings)))
    of_value = c(of_value, of_value1)
    convergence_criteria[i] = base::abs(of_value1 - of_value0)
    loop_end_time = base::Sys.time()
    loop_execution_time = difftime(loop_end_time, loop_start_time, units='secs')
    if(verbose)
    {
      message('Iteration ', i, ' completed in ', base::round(loop_execution_time, 3), ' seconds',
              '\n', 'Convergence criteria value: ', base::round(convergence_criteria[i], 4),
              '\n', 'Objective function (i-1) value: ', base::round(of_value0, 4), ' Objective function i value: ', round(of_value1, 4))
      message('___________________________________________________')
    }
    if(convergence_criteria[i] < tol)
    {
      if(verbose){message('Convergence criteria is satisfied. Value: ', base::round(convergence_criteria[i], 4))}
      break
    }
    scores0 = scores
    of_value0 = of_value1
    # Measure computation time
  }
  if(i == n_iters){warning('Algorithm did not converge in ', n_iters, ' iterations')}
  # Step 4: Normalize variables (process extracted from paper Chen et al. 2020)
  Fhat = loadings
  Lhat = scores
  sigmaF = (1/n_time) * t(Fhat) %*% Fhat
  sigmaA = (1/n_obs) * t(Lhat) %*% Lhat
  dum1 = expm::sqrtm(sigmaF) %*% sigmaA %*% expm::sqrtm(sigmaF)
  svd_decomp = base::svd(dum1)
  R = expm::sqrtm(base::solve(sigmaF)) %*% svd_decomp$u
  Fhat = Fhat %*% R
  Lhat = Lhat %*% t(base::solve(R))
  global_end_time = base::Sys.time()
  global_execution_time = difftime(global_end_time, global_start_time, units='secs')
  results = new_qfa(
    loadings = Fhat,
    scores = Lhat,
    unnormalized_loadings = loadings,
    normalization_matrix = R,
    objective_function_value = of_value1,
    list_objective_function_values = of_value,
    quantile_value = quantile_value,
    n_components = n_components,
    n_iters = i,
    execution_time = global_execution_time
  )
  return(results)
}

# PREDICTIONS -----------------------------------------------------------------

#' Predict QFA
#'
#' ## S3 method for class 'qfa'
#' Given a new matrix x, predicts the value of the scores associated to the given matrix
#'
#' @param object An object output of the qfa function.
#' @param newdata The N by T matrix of observed time instants to be tested.
#' @param ... further arguments passed to or from other methods.
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
#' results = qfa(x=x[1:150,], n_components=1, quantile_value=0.5)
#'
#' predictions = predict(results, newdata=x[151:200,])
predict.qfa = function(object, newdata, ...)
{
  if (!inherits(object, "qfa"))
  {
    stop('The object must be of class qfa')
  }
  # Unload required information from fqfa object
  if(!is.data.frame(newdata) & !is.matrix(newdata)){stop('x is not of a valid type. Object provided: ', typeof(newdata))}
  newdata = unname(as.matrix(newdata))
  if(ncol(newdata) == 1){ newdata = t(newdata)}

  normalization_matrix=object$normalization_matrix
  quantile_value = object$quantile_value
  unnormalized_loadings = object$unnormalized_loadings
  n_components = ncol(unnormalized_loadings)
  n_obs = nrow(newdata)
  scores =  matrix(0, n_obs, n_components)

  for(i in seq(n_obs))
  {
    scores[i,] =  quantreg::rq(newdata[i, ] ~ -1 + unnormalized_loadings, tau=quantile_value)$coefficients
  }
  scores = scores %*% t(base::solve(normalization_matrix))
  return(scores)
}

#' Plot QFA loadings
#'
#' ## S3 method for class 'qfa'
#' Given a fqpca object, plot the loadings
#'
#' @param x An object output of the fqpca function.
#' @param ... further arguments passed to or from other methods.
#'
#' @return The plot of loadings
#' @export
#'
#' @examples
#' # Generate fake dataset with 200 observations and 144 time points
#'
#' x = matrix(rep(sin(seq(0, 2*pi, length.out=144)), 200), byrow=TRUE, nrow=200)
#' x = x + matrix(rnorm(200*144, 0, 0.4), nrow=200)
#'
#' results = qfa(x=x[1:150,], n_components=1, quantile_value=0.5)
#'
#' plot(results)
plot.qfa = function(x, ...)
{
  loadings = x$loadings
  graphics::matplot(loadings, type = "l", lty = 1, lwd = 2)
  graphics::legend("topright", legend = 1:ncol(loadings), col=1:ncol(loadings), lty = 1, lwd = 2)
  graphics::title('Loadings')
}
