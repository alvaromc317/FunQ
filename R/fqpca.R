# FUNCTIONAL QUANTILE PCA ALGORITHM -------------------------------------------

# INNER LOOP FUNCTIONS --------------------------------------------------------

#' Compute objective value function
#'
#' Inner function to compute the objective value
#'
#' @param x The N by T matrix of observed time instants.
#' @param quantile_value The quantile considered.
#' @param scores The matrix of estimated scores
#' @param loadings The matrix of estimated loadings
#'
#' @return The objective value function
compute_objective_value = function(x, quantile_value, scores, loadings)
{
  x_reconstruction = scores %*% t(loadings)
  objective_value = base::sum(quantile_function(quantile_value=quantile_value, x=(x - x_reconstruction)), na.rm=T)
  objective_value = objective_value / sum(!is.na(x))
  return(objective_value)
}

#' Compute splines coefficients
#'
#' Inner function to compute the splines coefficients of the fqpca methodology.
#'
#' @param x The N by T matrix of observed time instants.
#' @param x_mask Mask matrix of the same dimensions as x indicating wich observations in x are known.
#' @param scores Initial matrix of scores.
#' @param spline_basis The spline basis matrix.
#' @param quantile_value The quantile considered. Default is the median 0.5.
#' @param method Method used in the resolution of the quantile regression model.
#'                It currently accepts the methods c('br', 'fn', 'pfn', 'sfn')
#'                along with any available solver in CVXR package like the free
#'                'SCS' or the comercial 'MOSEK'.
#' @param alpha_ridge The hyper-parameter controlling the penalization on the splines. This parameter has no effect if parameter penalized is set to FALSE. Default is 1e-12.
#' @param R_block Block diagonal matrix made of quadratic matrices used to apply a ridge based penalty.  This only has effect if parameter 'penalized' is TRUE.
#' @param valid_in_quantreg Array of valid quantreg methods
#'
#' @return The matrix of spline coefficients splines coefficients
compute_spline_coefficients = function(x, x_mask, scores, spline_basis, quantile_value, method, alpha_ridge, R_block, valid_in_quantreg)
{
  # At each iteration, call this function to obtain the estimation of F and Lambda
  n_obs = base::nrow(x)
  n_components = base::ncol(scores)
  n_basis = base::ncol(spline_basis)
  # vectorize X [x11, ..., x1T, x21, ..., x2T,...XNT]
  x_vector = c()
  for(i in base::seq(n_obs))
  {
    xi = x[i, x_mask[i,]]
    x_vector = c(x_vector, xi)
  }
  # Compute tensor product of Lambda and spline basis
  tensor_matrix = base::matrix(0, nrow=base::length(x_vector), ncol=n_components * n_basis)
  row_idx = 1
  for(i in base::seq(n_obs))
  {
    lambda_i = base::matrix(scores[i,], nrow=1)
    tmp_splines = spline_basis[x_mask[i,], ]
    n_time_i = base::nrow(tmp_splines)
    tensor_matrix[row_idx:(row_idx+n_time_i-1),] = base::kronecker(lambda_i, tmp_splines)
    row_idx = row_idx+n_time_i
  }

  if(method %in% valid_in_quantreg)
  {
    B_vector = quantreg::rq(x_vector ~ -1+tensor_matrix, tau=quantile_value, method=method)$coefficients
  } else{
    B_vector = quantile_regression_ridge(x=tensor_matrix, y=x_vector,  R=R_block, quantile_value=quantile_value, lambda=alpha_ridge, solver=method)
  }
  spline_coefficients = base::matrix(B_vector, byrow=FALSE, ncol=n_components)
  return(spline_coefficients)
}

#' Compute loadings (aka principal components)
#'
#' Inner function to compute the loadings of the fqpca methodology.
#'
#' @param spline_basis The spline basis matrix.
#' @param spline_coefficients the matrix of spline coefficients
#'
#' @return The matrix of loadings
compute_loadings = function(spline_basis, spline_coefficients)
{
  loadings = spline_basis %*% spline_coefficients
  return(loadings)
}

#' Compute scores
#'
#' Inner function to compute the scores of the fqpca methodology.
#'
#' @param x The N by T matrix of observed time instants.
#' @param x_mask Mask matrix of the same dimensions as x indicating wich observations in x are known.
#' @param loadings Matrix of loading coefficients
#' @param quantile_value The quantile considered.
#'
#' @return The matrix of scores
compute_scores = function(x, x_mask, loadings, quantile_value)
{
  n_obs = base::nrow(x)
  n_components = base::ncol(loadings)
  # Initialize the matrix of scores
  scores =  base::matrix(0, n_obs, n_components-1)
  for(i in base::seq(n_obs))
  {
    # Obtain the loadings associated to observation i with a specific time grid
    loadings_obs = loadings[x_mask[i, ], ]
    # Obtain the observation removing missing data
    xi = x[i, x_mask[i, ]]
    # Substract intercept
    xi = xi - loadings_obs[,1]

    scores[i,] =  quantreg::rq(xi ~ -1 + loadings_obs[,2:n_components], tau=quantile_value)$coefficients
  }
  # Add intercept column
  scores = base::cbind(1, scores)
  return(scores)
}

# OUTER LOOP FUNCTIONS --------------------------------------------------------

#' FQPCA normalization
#'
#' Performs the normalization of matrices of loadings and scores in order to ensure the solution is unique.
#'
#' @param loadings Matrix of loadings
#' @param scores Matrix of scores
#'
#' @return The normalized matrices of loadings and scores and the normalization matrix.
algorithm_normalization = function(loadings, scores)
{
  # Normalization process to ensure an identifiable solution
  n_time = nrow(loadings)
  n_obs = nrow(scores)
  n_components = ncol(loadings)-1

  # Remove intercept from normalization
  loadings_no_intercept = matrix(loadings[, 2:ncol(loadings)], ncol=n_components)
  scores_no_intercept = matrix(scores[, 2:ncol(scores)], ncol=n_components)

  # Ensure scores are mean-centered and move the effect of the mean scores to the intercept
  means_scores = matrix(colMeans(scores_no_intercept), ncol=1)
  scores_no_intercept = scale(scores_no_intercept, center=TRUE, scale=FALSE)
  loadings[, 1] = loadings[, 1] + loadings_no_intercept %*% means_scores

  # Make rotations to ensure identifiability
  sigmaF = (1/n_time) * t(loadings_no_intercept) %*% loadings_no_intercept
  sigmaA = (1/n_obs) * t(scores_no_intercept) %*% scores_no_intercept
  dum1 = expm::sqrtm(sigmaF) %*% sigmaA %*% expm::sqrtm(sigmaF)
  svd_decomp = base::svd(dum1)
  R = expm::sqrtm(base::solve(sigmaF)) %*% svd_decomp$u
  loadings_no_intercept = loadings_no_intercept %*% R;
  scores_no_intercept = scores_no_intercept %*% t(base::solve(R))

  # Add intercept
  loadings = cbind(loadings[,1], loadings_no_intercept)
  scores = cbind(scores[,1], scores_no_intercept)
  results = list(loadings=loadings, scores=scores, normalization_matrix=R)
  return(results)
}

#' FQPCA explained variability
#'
#' Computes the percentage of explained variability based on the variance of the scores matrix
#'
#' @param scores Matrix of scores
#'
#' @return The percentage of variability each component is explaining.
compute_explained_variability = function(scores)
{
  if(base::ncol(scores) == 2)
  {
    scores = matrix(scores[,-1], ncol=1)
  } else{
    scores = scores[,-1]
  }
  variability = base::diag(stats::var(scores))
  percentage_variability = variability / base::sum(variability)
  return(percentage_variability)
}

# MAIN ------------------------------------------------------------------------

new_fqpca <- function(loadings, scores, explained_variability, objective_function_value,
                      list_objective_function_values, execution_time, function_warnings,
                      mean_scores_unnormalized, n_components, quantile_value, method,
                      alpha_ridge, periodic, splines_df, tol, n_iters, verbose, seed,
                      normalization_matrix, spline_coefficients, spline_basis)
{
  structure(list(
    loadings=loadings,
    scores=scores,
    explained_variability=explained_variability,
    objective_function_value=objective_function_value,
    list_objective_function_values=list_objective_function_values,
    execution_time=execution_time,
    function_warnings=function_warnings,
    mean_scores_unnormalized=mean_scores_unnormalized,
    n_components=n_components,
    quantile_value=quantile_value,
    method=method,
    alpha_ridge=alpha_ridge,
    periodic=periodic,
    splines_df=splines_df,
    tol=tol,
    n_iters=n_iters,
    verbose=verbose,
    seed=seed,
    normalization_matrix=normalization_matrix,
    spline_coefficients=spline_coefficients,
    spline_basis=spline_basis
  ),
  class = "fqpca_object")
}


#' FQPCA
#'
#' Solves the functional quantile principal component analysis methodology
#'
#' @param x The N by T matrix of observed time instants
#' @param n_components Default=2. The number of estimated components.
#' @param quantile_value Default=0.5. The quantile considered.
#' @param periodic Default=TRUE. Boolean indicating if the data is expected to
#'                  be periodic (start coincides with end) or not.
#' @param splines_df Default=10. Degrees of freedom for the splines.
#' @param method Default='fn'. Method used in the resolution of the quantile
#'                regression model. It currently accepts the methods
#'                c('br', 'fn', 'pfn', 'sfn') along with any available solver in
#'                CVXR package like the free 'SCS' or the comercial 'MOSEK'.
#' @param alpha_ridge  Default=0. Hyper parameter controlling the penalization
#'                      on the second derivative of the splines. It has effect
#'                      only with CVXR solvers.
#' @param tol Default=1e-3. Tolerance on the convergence of the algorithm.
#'             Smaller values can speed up computation but may affect the
#'             quality of the estimations.
#' @param n_iters Default=30. Maximum number of iterations.
#' @param verbose Default=FALSE. Boolean indicating the verbosity.
#' @param seed Default=NULL. Seed for the random generator number.
#'              Parameter included for reproducibility purposes.
#'
#' @return fqpca_object
#' @export
#'
#' @examples
#' # Generate fake dataset with 150 observations and 144 time points
#'
#' x = matrix(rep(sin(seq(0, 2*pi, length.out=144)), 150), byrow=TRUE, nrow=150)
#' x = x + matrix(rnorm(150*144, 0, 0.4), nrow=150)
#'
#' # Add missing observations
#' x[sample(150*144, as.integer(0.2*150*144))] = NA
#'
#' results = fqpca(x=x, n_components=2, quantile_value=0.5)
#'
#' loadings = results$loadings
#' scores = results$scores
fqpca = function(x, n_components=2,  quantile_value=0.5,  periodic=TRUE, splines_df=10, method='fn', alpha_ridge=0, tol=1e-3, n_iters=30, verbose=FALSE, seed=NULL)
{
  global_start_time = base::Sys.time()

  # BASIC COMPROBATIONS -------------------------------------------------------

  if(!is.data.frame(x) & !is.matrix(x)){stop('x is not of a valid type. Object provided: ', typeof(x))}
  if(!n_components == floor(n_components)){stop('n_components must be an integer number. Value provided: ', n_components)}
  if(quantile_value<0 | quantile_value>1){stop('quantile_value must be a value between 0 and 1. Value provided: ', quantile_value)}
  if(!splines_df == floor(splines_df)){stop('splines_df must be an integer number. Value provided: ', splines_df)}
  if(!n_iters == floor(n_iters)){stop('n_iters must be an integer number. Value provided: ', n_iters)}
  if(alpha_ridge<0){stop('alpha_ridge must be equal or larger than 0')}
  if(tol<0){stop('tol must be equal or larger than 0')}
  if(!n_iters == floor(n_iters)){stop('n_iters must be an integer number. Value provided: ', n_iters)}

  valid_in_quantreg = c('br', 'fn', 'sfn', 'pfn')
  valid_in_cvxr = CVXR::installed_solvers()
  if(method %in% valid_in_cvxr)
  {
    penalized = TRUE
  } else if(method %in% valid_in_quantreg)
  {
    penalized = FALSE
  } else{
    stop('Invalid method provided. Valid methods include the quantreg methods
         c("br", "fn", "sfn", "pfn") and any installed solver for CVXR.
         Value provided: ', method)
  }

  if(!base::is.null(seed)){base::set.seed(seed)}
  x = unname(base::as.matrix(x))
  # Step 1: Initialize values
  n_obs = base::dim(x)[1]
  n_time = base::dim(x)[2]
  x_axis = base::seq(0, 1, length.out = n_time)
  x_mask = !base::is.na(x)

  # INITIALIZATION OF SPLINE BASIS --------------------------------------------

  # Initialize splines and second derivative ---
  if(periodic)
  {
    knots = base::seq(0, 1, length.out = 2 + splines_df - 1)  # 2 boundary knots + df - intercept
  } else{
    knots = base::seq(0, 1, length.out = 2 + splines_df - 3 - 1) # 2 boundary knots + df - degree - intercept
  }
  knots = knots[2:(base::length(knots)-1)]
  spline_basis =  splines2::mSpline(x_axis, degree=3, intercept=TRUE, knots=knots, periodic=periodic, Boundary.knots=c(min(x_axis), max(x_axis)))

  if(penalized)
  {
    deriv_basis = stats::deriv(spline_basis, derivs=2)
    R = t(deriv_basis) %*% deriv_basis
    R_block = simex::diag.block(R, n_components + 1)
    # Ensure that R_block is semi definite positive: (all eigenvalues are greater or equal to 0)
    for(i in 1:5000)
    {
      R_block = (1 - 1e-3) * R_block + 1e-3 * base::diag(base::nrow(R_block))
      if(all(base::eigen(R_block)$values > 0)){break}
    }
    if(!all(base::eigen(R_block)$values > 0)) {stop('Splines penalization matrix is not positive definite')}
  } else{R_block = NULL}

  # INITIALIZATION OF SCORES AND LOADINGS -------------------------------------

  # Generate initial splines basis coefficients
  spline_coefficients_0 = base::matrix(stats::rnorm(base::dim(spline_basis)[2] * (n_components+1), mean=0, sd=1/n_time),
                            nrow=base::dim(spline_basis)[2], ncol=n_components+1)

  # Initialize scores and loadings with these spline coefficients
  loadings_0 = compute_loadings(spline_basis=spline_basis, spline_coefficients=spline_coefficients_0)
  scores_0 = compute_scores(x=x, x_mask=x_mask, loadings=loadings_0, quantile_value=quantile_value)

  # Compute objective value function
  objective_function_0 = compute_objective_value(quantile_value=quantile_value, x=x, scores=scores_0, loadings=loadings_0)

  # Initialize storage of final solution
  objective_function_array = objective_function_0
  best_results = list(loadings=loadings_0, scores=scores_0, objective_function=objective_function_0, spline_coefficients=spline_coefficients_0, iteration=0)

  # Iterate until convergence or until max_iters is reached
  convergence_criteria = c()
  function_warnings = list(splines=FALSE, scores=FALSE, diverging_loop=FALSE, normalization=FALSE)

  for(i in 1:n_iters)
  {
    loop_start_time = base::Sys.time()

    # COMPUTATION OF LOADINGS AND SCORES ---------------------------

    # Obtain splines coefficients
    spline_coefficients_1 = try(compute_spline_coefficients(x=x, x_mask=x_mask, scores=scores_0, spline_basis=spline_basis, quantile_value=quantile_value, method=method, alpha_ridge=alpha_ridge, R_block=R_block, valid_in_quantreg=valid_in_quantreg), silent=TRUE)
    if(!is.matrix(spline_coefficients_1))
    {
      warning('Computation of spline coefficients failed at iteration ', i, '. Providing results from previous iteration.')
      function_warnings$splines = TRUE
      break
    }

    # Compute loadings
    loadings_1 = compute_loadings(spline_basis=spline_basis, spline_coefficients=spline_coefficients_1)

    # Compute scores
    scores_1 = try(compute_scores(x=x, x_mask=x_mask, loadings=loadings_1, quantile_value=quantile_value), silent=FALSE)
    if(!is.matrix(scores_1))
    {
      warning('Computation of scores failed at iteration ', i, '. Providing results from previous iteration.')
      function_warnings$scores = TRUE
      break
    }

    # Compute objective function
    objective_function_1 = compute_objective_value(quantile_value=quantile_value, x=x, scores=scores_1, loadings=loadings_1)
    objective_function_array = c(objective_function_array, objective_function_1)

    # COMPROBATION OF CONVERGENCE ---------------------------------------------

    if(objective_function_1 < best_results$objective_function)
    {
      best_results = list(loadings=loadings_1, scores=scores_1, objective_function=objective_function_1, spline_coefficients=spline_coefficients_1, iteration=i)
    }

    convergence_criteria[i] = base::abs(objective_function_1 - objective_function_0)

    # Measure computation time
    loop_end_time = base::Sys.time()
    loop_execution_time = difftime(loop_end_time, loop_start_time, units='secs')

    if(verbose)
    {
      message('Iteration ', i, ' completed in ', base::round(loop_execution_time, 3), ' seconds',
              '\n', 'Convergence criteria value: ', base::round(convergence_criteria[i], 4),
              '\n', 'Objective function (i-1) value: ', base::round(objective_function_0, 4), ' Objective function i value: ', round(objective_function_1, 4))
      message('___________________________________________________')
    }

    if(convergence_criteria[i] < tol)
    {
      if(verbose){message('Convergence criteria is satisfied. Value: ', base::round(convergence_criteria[i], 4))}
      break
    }

    # If the objective value function (that should be decreasing) is
    # 100 times larger than the original objective value function, break
    if( objective_function_array[i+1] > 100 * objective_function_array[2])
    {
      function_warnings$diverging_loop = TRUE
      if(verbose)
      {
        message('Breaking diverging loop at iteration ', i,
                '\n', 'Objective function value: ', base::round(objective_function_array[i], 4),
                '\n', 'Initial objective value function: ', base::round(objective_function_array[2], 4))
        message('Results of the iteration with the smallest objective value: ', base::round(best_results$objective_function, 4), ' will be provided')
      }
      break
    }

    # UPDATE VALUES -----------------------------------------------------------

    objective_function_0 = objective_function_1
    scores_0 = scores_1
  }
  if(i == n_iters){ if(verbose) {warning('Algorithm reached maximum number of iterations without convergence: ', n_iters, ' iterations')}}

  # NORMALIZATION TO ENSURE IDENTIFIABILITY -----------------------------------

  mean_scores = colMeans(matrix(best_results$scores[,-1], ncol=n_components)) # Required in the prediction step
  normalization_result = try(algorithm_normalization(loadings=best_results$loadings, scores=best_results$scores), silent=TRUE)
  error_chr = class(normalization_result)
  if(error_chr=="try-error")
  {
    warning('Failed normalization step to ensure identifiability. Providing unnormalized results')
    function_warnings$normalization = TRUE
    normalization_result = list(loadings=best_results$loadings, scores=best_results$scores, normalization_matrix=diag(n_components))
  }

  # EXPLAINED VARIABILITY -----------------------------------------------------

  explained_variability = compute_explained_variability(normalization_result$scores)

  global_end_time = base::Sys.time()
  global_execution_time = difftime(global_end_time, global_start_time, units='secs')

  results = new_fqpca(
    loadings=normalization_result$loadings,
    scores=normalization_result$scores,
    explained_variability=explained_variability,
    objective_function_value=best_results$objective_function,
    list_objective_function_values=objective_function_array,
    execution_time=global_execution_time,
    function_warnings=function_warnings,
    mean_scores_unnormalized = mean_scores,
    n_components=n_components,
    quantile_value=quantile_value,
    method=method,
    alpha_ridge=alpha_ridge,
    periodic=periodic,
    splines_df=splines_df,
    tol=tol,
    n_iters=best_results$iteration,
    verbose=verbose,
    seed=seed,
    normalization_matrix=normalization_result$normalization_matrix,
    spline_coefficients=best_results$spline_coefficients,
    spline_basis=spline_basis
  )
  return(results)
}


# PREDICTIONS -----------------------------------------------------------------

#' Predict FQPCA
#'
#' ## S3 method for class 'fqpca_object'
#' Given a new matrix x, predicts the value of the scores associated to the given matrix
#'
#' @param object An object output of the fqpca function.
#' @param newdata The N by T matrix of observed time instants to be tested.
#' @param ... further arguments passed to or from other methods.
#'
#' @return The normalized matrix of scores.
#' @export
#'
#' @examples
#' # Generate fake dataset with 150 observations and 144 time points
#'
#' x = matrix(rep(sin(seq(0, 2*pi, length.out=144)), 150), byrow=TRUE, nrow=150)
#' x = x + matrix(rnorm(150*144, 0, 0.4), nrow=150)
#'
#' # Add missing observations
#' x[sample(150*144, as.integer(0.2*150*144))] = NA
#'
#' results = fqpca(x=x[1:100,], n_components=2, quantile_value=0.5)
#'
#' predictions = predict(object=results, newdata=x[101:150,])
predict.fqpca_object = function(object, newdata, ...)
{
  if (!inherits(object, "fqpca_object"))
  {
    stop('The object must be of class fqpca_object')
  }
  if(!is.data.frame(newdata) & !is.matrix(newdata)){stop('newdata is not of a valid type. Object provided: ', typeof(newdata), '\nValid types: data.frame or matrix')}
  newdata = unname(as.matrix(newdata))
  if(ncol(newdata) == 1){ newdata = t(newdata)}

  x_mask = !base::is.na(newdata)
  n_obs = base::nrow(newdata)
  n_time = base::ncol(newdata)
  if(sum(!x_mask) == n_obs * n_time){stop('newdata contains no information. Check that there are values different from NA')}

  # Estimation of unnormalized loadings
  loadings_1 = compute_loadings(spline_basis=object$spline_basis, spline_coefficients=object$spline_coefficients)

  # Estimation of scores
  scores_1 = try(compute_scores(x=newdata, x_mask=x_mask, loadings=loadings_1, quantile_value=object$quantile_value), silent=FALSE)
  if(!is.matrix(scores_1))
  {
    stop('Computation of scores failed.')
  }
  scores_1 = matrix(scores_1[,-1], ncol=object$n_components) # Remove intercept column
  scores_1 = scores_1 - matrix(rep(object$mean_scores_unnormalized, n_obs), nrow=n_obs, byrow=T)
  scores_1 = scores_1 %*% t(base::solve(object$normalization_matrix)) # Normalize
  scores_1 = base::cbind(1, scores_1) # Add intercept column
  return(scores_1)
}

# BASIC PLOT ------------------------------------------------------------------

#' Plot FQPCA loading functions
#'
#' ## S3 method for class 'fqpca'
#' Given a fqpca object, plot the loading functions
#'
#' @param x An object output of the fqpca function.
#' @param ... further arguments passed to or from other methods.
#'
#' @return The plot of loadings
#' @export
#'
#' @examples
#' # Generate fake dataset with 150 observations and 144 time points
#'
#' x = matrix(rep(sin(seq(0, 2*pi, length.out=144)), 150), byrow=TRUE, nrow=150)
#' x = x + matrix(rnorm(150*144, 0, 0.4), nrow=150)
#'
#' # Add missing observations
#' x[sample(150*144, as.integer(0.2*150*144))] = NA
#'
#' results = fqpca(x=x, n_components=2, quantile_value=0.5)
#'
#' plot(results)
plot.fqpca_object = function(x, ...)
{
  if (!inherits(x, "fqpca_object"))
  {
    stop('The object must be of class fqpca_object')
  }
  loadings = x$loadings
  graphics::par(mfrow = c(1, 2))
  graphics::plot(loadings[,1], type='l', lty=1, lwd=2, xlab='', ylab='')
  graphics::title('Intercept curve')

  no_intercept_loadings = loadings[,-1]
  if(base::ncol(loadings) == 2)
  {
    no_intercept_loadings = base::matrix(no_intercept_loadings, ncol=1)
  }
  graphics::matplot(no_intercept_loadings, type = "l", lty = 1, lwd = 2, xlab='', ylab='')
  graphics::title('Principal components functions')
}
