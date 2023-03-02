
#' FQPCA normalization
#'
#' Performs the normalization of matrices of loadings and scores in order to ensure the solution is unique.
#'
#' @param loadings Matrix of loadings
#' @param scores Matrix of scores
#'
#' @return The normalized matrices of loadings and scores and the normalization matrix.
algorithm_normalization_original = function(loadings, scores)
{
  # Normalization process to ensure an identifiable solution
  n_time = nrow(loadings)
  n_obs = nrow(scores)
  n_components = ncol(loadings)-1

  # Remove intercept from normalization
  loadings_no_intercept = matrix(loadings[, 2:ncol(loadings)], ncol=n_components)
  scores_no_intercept = matrix(scores[, 2:ncol(scores)], ncol=n_components)

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

new_fqpca_original <- function(loadings, scores, explained_variability, objective_function_value,
                      list_objective_function_values, execution_time, function_warnings,
                      n_components, quantile_value, alpha_ridge,
                      periodic, splines_df, tol, n_iters, solver, penalized, method, verbose, seed,
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
    n_components=n_components,
    quantile_value=quantile_value,
    alpha_ridge=alpha_ridge,
    periodic=periodic,
    splines_df=splines_df,
    tol=tol,
    n_iters=n_iters,
    solver=solver,
    penalized=penalized,
    method=method,
    verbose=verbose,
    seed=seed,
    normalization_matrix=normalization_matrix,
    spline_coefficients=spline_coefficients,
    spline_basis=spline_basis
  ),
  class = "fqpca_object_original")
}

#' FQPCA
#'
#' Solves the functional quantile principal component analysis methodology
#'
#' @param x The N by T matrix of observed time instants
#' @param n_components The number of estimated components. Default is 2.
#' @param quantile_value The quantile considered. Default is the median 0.5.
#' @param alpha_ridge The hyperparameter controlling the penalization on the splines. This parameter has no effect if parameter penalized is set to FALSE. Default is 1e-16.
#' @param periodic Boolean indicating if the data is expected to be periodic (start coincides with end) or not. Default is TRUE.
#' @param splines_df Degrees of freedom for the splines. It is recommended not to modify this parameter and control the smoothness via the hyperparameter alpha_ridge. Default is 30.
#' @param tol Tolerance on the convergence of the algorithm. Smaller values can spped up computation but may affect the quality of the estimations. Default is 1e-3.
#' @param n_iters Maximum number of iterations. Default is 30.
#' @param solver Solver to be used by the CVXR package in the resolution of the penalized quantile regression models. The speed of the methodology can vary depending on the solver used. There are free alternatives available like 'SCS', 'ECOS' or 'OSQP'. It is also possible to use licensed programs like 'MOSEK' or 'GUROBI'. Default is 'SCS'.
#' @param penalized Boolean indicating if the model should solve a penalized splines model or a non penelized splines model
#' @param method Method to be used by the quantreg package when solving quantile regression models. This parameter has no effect if penalized is TRUE (default)
#' @param verbose Boolean indicating the verbosity of the function.
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
#' # Add missing observations
#' x[sample(200*144, as.integer(0.2*200*144))] = NA
#'
#' results = fqpca_rank(x=x, n_components=2, quantile_value=0.5, alpha_ridge=1e-12)
#'
#' loadings = results$loadings
#' scores = results$scores
fqpca_rank = function(x, n_components=2,  quantile_value=0.5, alpha_ridge=1e-16,  periodic=TRUE, splines_df=30, tol=1e-3, n_iters=30, solver='SCS', penalized=TRUE, method='br', verbose=FALSE, seed=1)
{
  if(!is.data.frame(x) & !is.matrix(x)){stop('x is not of a valid type. Object provided: ', typeof(x))}
  if(!n_components == floor(n_components)){stop('n_components must be an integer number. Value provided: ', n_components)}
  if(quantile_value<0 | quantile_value>1){stop('quantile_value must be a value between 0 and 1. Value provided: ', quantile_value)}
  if(!n_iters == floor(n_iters)){stop('n_iters must be an integer number. Value provided: ', n_iters)}

  global_start_time = base::Sys.time()

  if(!base::is.null(seed)){base::set.seed(seed)}
  x = unname(base::as.matrix(x))
  # Step 1: Initialize values
  n_obs = base::dim(x)[1]
  n_time = base::dim(x)[2]
  x_axis = base::seq(0, 1, length.out = n_time)
  x_mask = !base::is.na(x)
  true_n_components = n_components

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

  # INITIALIZATION OF SCORES AND LOADINGS -------------------------------------

  # Generate initial splines basis coefficients
  spline_coefficients_0 = base::matrix(stats::rnorm(base::dim(spline_basis)[2] * (n_components+1), mean=0, sd=1/n_time),
                            nrow=base::dim(spline_basis)[2], ncol=n_components+1)

  # Initialize scores and loadings with these spline coefficients
  loadings_0 = compute_loadings(spline_basis=spline_basis, spline_coefficients=spline_coefficients_0)
  scores_0 = compute_scores(x=x, x_mask=x_mask, loadings=loadings_0, quantile_value=quantile_value, method=method)

  # Compute objective value function
  objective_function_0 = compute_objective_value(quantile_value=quantile_value, x=x, scores=scores_0, loadings=loadings_0)

  # Initialize storage of final solution
  objective_function_array = objective_function_0
  best_results = list(loadings=loadings_0, scores=scores_0, objective_function=objective_function_0, spline_coefficients=spline_coefficients_0, iteration=0)

  # Iterate until convergence or until max_iters is reached
  convergence_criteria = c()
  function_warnings = list(splines=FALSE, scores=FALSE, diverging_loop=FALSE, normalization=FALSE, loadings_singularity=FALSE)

  for(i in 1:n_iters)
  {
    loop_start_time = base::Sys.time()

    # COMPUTATION OF LOADINGS AND SCORES ---------------------------

    # Obtain splines coefficients
    spline_coefficients_1 = try(compute_spline_coefficients(x=x, x_mask=x_mask, scores=scores_0, spline_basis=spline_basis, quantile_value=quantile_value, alpha_ridge=alpha_ridge, R_block=R_block, solver=solver, penalized=penalized, method=method), silent=TRUE)
    if(!is.matrix(spline_coefficients_1))
    {
      warning('Computation of spline coefficients failed at iteration ', i, '. Providing results from previous iteration.')
      function_warnings$splines = TRUE
      break
    }

    # Compute loadings
    loadings_1 = compute_loadings(spline_basis=spline_basis, spline_coefficients=spline_coefficients_1)

    # Normalization result
    normalization_result = try(algorithm_normalization_original(loadings=loadings_1, scores=scores_0), silent=TRUE)
    error_chr = class(normalization_result)
    if(error_chr=="try-error")
    {
      warning('Failed normalization step to ensure identifiability at iteration ', i, ' Providing unnormalized results')
      normalization_result = list(loadings=loadings_1, scores=scores_0, normalization_matrix=diag(n_components))
    }

    loadings_1 = normalization_result$loadings

    rank_loadings = base::qr(loadings_1)$rank
    if(rank_loadings != n_components+1)
    {
      function_warnings$loadings_singularity = TRUE
      n_components = rank_loadings-1
      loadings_1 = base::matrix(loadings_1[, 1:rank_loadings], ncol=rank_loadings)
    }

    # Compute scores
    scores_1 = try(compute_scores(x=x, x_mask=x_mask, loadings=loadings_1, quantile_value=quantile_value, method=method), silent=FALSE)
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
      best_results = list(loadings=loadings_1, scores=scores_1, objective_function=objective_function_1, spline_coefficients=spline_coefficients_1, iteration=i, n_components=ncol(loadings_1)-1)
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
  if(true_n_components != best_results$n_components){warning('Loadings matrix is singular. Reducing the number of components from ',true_n_components, ' to ',best_results$n_components, ' to obtain a full-rank matrix')}

  # NORMALIZATION TO ENSURE IDENTIFIABILITY -----------------------------------

  normalization_result = try(algorithm_normalization_original(loadings=best_results$loadings, scores=best_results$scores), silent=TRUE)
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

  results = new_fqpca_original(
    loadings=normalization_result$loadings,
    scores=normalization_result$scores,
    explained_variability=explained_variability,
    objective_function_value=best_results$objective_function,
    list_objective_function_values=objective_function_array,
    execution_time=global_execution_time,
    function_warnings=function_warnings,
    n_components=best_results$n_components,
    quantile_value=quantile_value,
    alpha_ridge=alpha_ridge,
    periodic=periodic,
    splines_df=splines_df,
    tol=tol,
    n_iters=best_results$iteration,
    solver=solver,
    penalized=penalized,
    method=method,
    verbose=verbose,
    seed=seed,
    normalization_matrix=normalization_result$normalization_matrix,
    spline_coefficients=best_results$spline_coefficients,
    spline_basis=spline_basis
  )
  return(results)
}
