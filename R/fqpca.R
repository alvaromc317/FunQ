# FUNCTIONAL QUANTILE PCA ALGORITHM -------------------------------------------

# INNER LOOP FUNCTIONS --------------------------------------------------------

#' Compute objective value function
#'
#' Inner function to compute the objective value
#'
#' @param Y The N by T matrix of observed time instants.
#' @param quantile.value The quantile considered.
#' @param scores The matrix of estimated scores
#' @param loadings The matrix of estimated loadings
#'
#' @return The objective value function
compute_objective_value = function(Y, quantile.value, scores, loadings)
{
  Y_reconstruction = scores %*% t(loadings)
  objective_value = base::sum(quantile_function(quantile.value=quantile.value, x=(Y - Y_reconstruction)), na.rm=T)
  objective_value = objective_value / sum(!is.na(Y))
  return(objective_value)
}

#' Compute splines coefficients
#'
#' Inner function to compute the splines coefficients of the fqpca methodology.
#'
#' @param Y The N by T matrix of observed time instants.
#' @param Y_mask Mask matrix of the same dimensions as Y indicating wich observations in Y are known.
#' @param scores Initial matrix of scores.
#' @param spline.basis The spline basis matrix.
#' @param quantile.value The quantile considered. Default is the median 0.5.
#' @param method Method used in the resolution of the quantile regression model.
#'                It currently accepts the methods c('br', 'fn', 'pfn', 'sfn')
#'                from quantreg package along with any available solver in CVXR
#'                package like the free 'SCS' or the commercial 'MOSEK'.
#' @param alpha.ridge The hyper-parameter controlling the penalization on the splines. This parameter has no effect if parameter penalized is set to FALSE. Default is 1e-12.
#' @param R_block Block diagonal matrix made of quadratic matrices used to apply a ridge based penalty.  This only has effect if parameter 'penalized' is TRUE.
#' @param valid_in_quantreg Array of valid quantreg methods
#'
#' @return The matrix of spline coefficients splines coefficients
compute_spline_coefficients = function(Y, Y_mask, scores, spline.basis, quantile.value, method, alpha.ridge, R_block, valid_in_quantreg)
{
  # At each iteration, call this function to obtain the estimation of F and Lambda
  n_obs = base::nrow(Y)
  npc = base::ncol(scores)
  n_basis = base::ncol(spline.basis)
  # vectorize Y [Y11, ..., Y1T, Y21, ..., Y2T,...YNT]
  Y_vector = c()
  for(i in base::seq(n_obs))
  {
    Yi = Y[i, Y_mask[i,]]
    Y_vector = c(Y_vector, Yi)
  }
  # Compute tensor product of Lambda and spline basis
  tensor_matrix = base::matrix(0, nrow=base::length(Y_vector), ncol=npc * n_basis)
  row_idx = 1
  for(i in base::seq(n_obs))
  {
    lambda_i = base::matrix(scores[i,], nrow=1)
    tmp_splines = spline.basis[Y_mask[i,], ]
    n_time_i = base::nrow(tmp_splines)
    tensor_matrix[row_idx:(row_idx+n_time_i-1),] = base::kronecker(lambda_i, tmp_splines)
    row_idx = row_idx+n_time_i
  }

  if(method %in% valid_in_quantreg)
  {
    B_vector = quantreg::rq(Y_vector ~ -1+tensor_matrix, tau=quantile.value, method=method)$coefficients
  } else{
    B_vector = quantile_regression_ridge(x=tensor_matrix, y=Y_vector,  R=R_block, quantile.value=quantile.value, lambda=alpha.ridge, solver=method)
  }
  spline.coefficients = base::matrix(B_vector, byrow=FALSE, ncol=npc)
  return(spline.coefficients)
}

#' Compute loadings (aka principal components)
#'
#' Inner function to compute the loadings of the fqpca methodology.
#'
#' @param spline.basis The spline basis matrix.
#' @param spline.coefficients the matrix of spline coefficients
#'
#' @return The matrix of loadings
compute_loadings = function(spline.basis, spline.coefficients)
{
  loadings = spline.basis %*% spline.coefficients
  return(loadings)
}

#' Compute scores
#'
#' Inner function to compute the scores of the fqpca methodology.
#'
#' @param Y The N by T matrix of observed time instants.
#' @param Y_mask Mask matrix of the same dimensions as Y indicating wich observations in Y are known.
#' @param loadings Matrix of loading coefficients
#' @param quantile.value The quantile considered.
#'
#' @return The matrix of scores
compute_scores = function(Y, Y_mask, loadings, quantile.value)
{
  n_obs = base::nrow(Y)
  npc = base::ncol(loadings)
  # Initialize the matrix of scores
  scores =  base::matrix(0, n_obs, npc-1)
  for(i in base::seq(n_obs))
  {
    # Obtain the loadings associated to observation i with a specific time grid
    loadings_obs = loadings[Y_mask[i, ], ]
    # Obtain the observation removing missing data
    Yi = Y[i, Y_mask[i, ]]
    # Substract intercept
    Yi = Yi - loadings_obs[,1]

    scores[i,] =  quantreg::rq(Yi ~ -1 + loadings_obs[,2:npc], tau=quantile.value)$coefficients
  }
  # Add intercept column
  scores = base::cbind(1, scores)
  return(scores)
}

# OUTER LOOP FUNCTIONS --------------------------------------------------------

#' FQPCA rotation
#'
#' Performs the rotation of matrices of loadings and scores in order to ensure the solution is unique.
#'
#' @param loadings Matrix of loadings
#' @param scores Matrix of scores
#'
#' @return The rotated matrices of loadings and scores and the rotation matrix.
rotate_scores_and_loadings = function(loadings, scores)
{
  n_time = nrow(loadings)
  n_obs = nrow(scores)
  npc = ncol(loadings)-1

  # Remove intercept from rotation
  loadings_no_intercept = matrix(loadings[, 2:ncol(loadings)], ncol=npc)
  scores_no_intercept = matrix(scores[, 2:ncol(scores)], ncol=npc)

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
  results = list(loadings=loadings, scores=scores, rotation.matrix=R)
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

new_fqpca <- function(loadings, scores, pve, objective.function.value,
                      list.objective.function.values, execution.time, function.warnings,
                      mean_scores_unnormalized, npc, quantile.value, method,
                      alpha.ridge, periodic, splines.df, tol, n.iters, verbose, seed,
                      rotation.matrix, spline.coefficients, spline.basis)
{
  structure(list(
    loadings=loadings,
    scores=scores,
    pve=pve,
    objective.function.value=objective.function.value,
    list.objective.function.values=list.objective.function.values,
    execution.time=execution.time,
    function.warnings=function.warnings,
    mean_scores_unnormalized=mean_scores_unnormalized,
    npc=npc,
    quantile.value=quantile.value,
    method=method,
    alpha.ridge=alpha.ridge,
    periodic=periodic,
    splines.df=splines.df,
    tol=tol,
    n.iters=n.iters,
    verbose=verbose,
    seed=seed,
    rotation.matrix=rotation.matrix,
    spline.coefficients=spline.coefficients,
    spline.basis=spline.basis
  ),
  class = "fqpca_object")
}


#' FQPCA
#'
#' Solves the functional quantile principal component analysis methodology
#'
#' @param Y The object storing the functional data. This can be either a matrix
#'          data.frame or tibble of dimensions N by T, or a tf object from the
#'          tidyfun package
#' @param npc Default=2 The number of estimated components.
#' @param quantile.value Default=0.5. The quantile considered.
#' @param periodic Default=TRUE. Boolean indicating if the data is expected to
#'                  be periodic (start coincides with end) or not.
#' @param splines.df Default=10. Degrees of freedom for the splines.
#' @param method Default='fn'. Method used in the resolution of the quantile
#'                regression model. It currently accepts the methods
#'                c('br', 'fn', 'pfn', 'sfn') from quantreg package along with
#'                any available solver in CVXR package like the free 'SCS'
#'                or the commercial 'MOSEK'.
#' @param alpha.ridge  Default=0. Hyper parameter controlling the penalization
#'                      on the second derivative of the splines. It has effect
#'                      only with CVXR solvers. Large values are associated with
#'                      smoother results. This component is experimental and may
#'                      lead to computational issues.
#' @param tol Default=1e-3. Tolerance on the convergence of the algorithm.
#'             Smaller values can speed up computation but may affect the
#'             quality of the estimations.
#' @param n.iters Default=30. Maximum number of iterations.
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
#' Y = matrix(rep(sin(seq(0, 2*pi, length.out=144)), 150), byrow=TRUE, nrow=150)
#' Y = Y + matrix(rnorm(150*144, 0, 0.4), nrow=150)
#'
#' # Add missing observations
#' Y[sample(150*144, as.integer(0.2*150*144))] = NA
#'
#' results = fqpca(Y=Y, npc=2, quantile.value=0.5)
#'
#' loadings = results$loadings
#' scores = results$scores
fqpca = function(Y, npc=2,  quantile.value=0.5,  periodic=TRUE, splines.df=10, method='fn', alpha.ridge=0, tol=1e-3, n.iters=30, verbose=FALSE, seed=NULL)
{
  global_start_time = base::Sys.time()

  # BASIC COMPROBATIONS -------------------------------------------------------

  if(!npc == floor(npc)){stop('npc must be an integer number. Value provided: ', npc)}
  if(quantile.value<0 | quantile.value>1){stop('quantile.value must be a value between 0 and 1. Value provided: ', quantile.value)}
  if(!splines.df == floor(splines.df)){stop('splines.df must be an integer number. Value provided: ', splines.df)}
  if(!n.iters == floor(n.iters)){stop('n.iters must be an integer number. Value provided: ', n.iters)}
  if(alpha.ridge<0){stop('alpha.ridge must be greater or equal than 0')}
  if(tol<0){stop('tol must be greater or equal than 0')}
  if(!n.iters == floor(n.iters)){stop('n.iters must be an integer number. Value provided: ', n.iters)}

  if(!is.data.frame(Y) & !is.matrix(Y))
  {
    # This structure allows tf to be a suggested package rather than mandatory
    if(!tf::is_tf(Y))
    {
      stop('Y is not of a valid type. Object provided: ', class(Y))
    }
  }

  # DETERMINE IF THE ALGORITHM RUNS USING PENALIZED OR UNPENALIZED SPLINES ----

  valid_in_quantreg = c('br', 'fn', 'sfn', 'pfn')
  if(method %in% valid_in_quantreg)
  {
    penalized = FALSE
  } else{
    # Check available solvers in CVXR
    valid_in_cvxr = CVXR::installed_solvers()
    if(method %in% valid_in_cvxr)
    {
      penalized = TRUE
    } else{
      stop('Invalid method provided. Valid methods include the quantreg methods
         c("br", "fn", "sfn", "pfn") and any installed solver for CVXR.
         Value provided: ', method)
    }
  }

  if(!base::is.null(seed)){base::set.seed(seed)}

  # Step 1: Initialize values
  Y = base::unname(base::as.matrix(Y))
  n_obs = base::dim(Y)[1]
  n_time = base::dim(Y)[2]
  Y_axis = base::seq(0, 1, length.out = n_time)
  Y_mask = !base::is.na(Y)

  # INITIALIZATION OF SPLINE BASIS --------------------------------------------

  # Initialize splines and second derivative ---
  if(periodic)
  {
    knots = base::seq(0, 1, length.out = 2 + splines.df - 1)  # 2 boundary knots + df - intercept
  } else{
    knots = base::seq(0, 1, length.out = 2 + splines.df - 3 - 1) # 2 boundary knots + df - degree - intercept
  }
  knots = knots[2:(base::length(knots)-1)]
  spline.basis =  splines2::mSpline(Y_axis, degree=3, intercept=TRUE, knots=knots, periodic=periodic, Boundary.knots=c(min(Y_axis), max(Y_axis)))

  if(penalized)
  {
    deriv_basis = stats::deriv(spline.basis, derivs=2)
    R = t(deriv_basis) %*% deriv_basis
    R_block = simex::diag.block(R, npc + 1)
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
  spline_coefficients_0 = base::matrix(stats::rnorm(base::dim(spline.basis)[2] * (npc+1), mean=0, sd=1/n_time),
                            nrow=base::dim(spline.basis)[2], ncol=npc+1)

  # Initialize scores and loadings with these spline coefficients
  loadings_0 = compute_loadings(spline.basis=spline.basis, spline.coefficients=spline_coefficients_0)
  scores_0 = compute_scores(Y=Y, Y_mask=Y_mask, loadings=loadings_0, quantile.value=quantile.value)

  # Compute objective value function
  objective_function_0 = compute_objective_value(quantile.value=quantile.value, Y=Y, scores=scores_0, loadings=loadings_0)

  # Initialize storage of final solution
  objective_function_array = objective_function_0
  best_results = list(loadings=loadings_0, scores=scores_0, objective_function=objective_function_0, spline.coefficients=spline_coefficients_0, iteration=0)

  # Iterate until convergence or until max_iters is reached
  convergence_criteria = c()
  function.warnings = list(splines=FALSE, scores=FALSE, diverging.loop=FALSE, rotation=FALSE)

  for(i in 1:n.iters)
  {
    loop_start_time = base::Sys.time()

    # COMPUTATION OF LOADINGS AND SCORES ---------------------------

    # Obtain splines coefficients
    spline_coefficients_1 = try(compute_spline_coefficients(Y=Y, Y_mask=Y_mask, scores=scores_0, spline.basis=spline.basis, quantile.value=quantile.value, method=method, alpha.ridge=alpha.ridge, R_block=R_block, valid_in_quantreg=valid_in_quantreg), silent=TRUE)
    if(!is.matrix(spline_coefficients_1))
    {
      warning('Computation of spline coefficients failed at iteration ', i, '. Providing results from previous iteration.')
      function.warnings$splines = TRUE
      break
    }

    # Compute loadings
    loadings_1 = compute_loadings(spline.basis=spline.basis, spline.coefficients=spline_coefficients_1)

    # Compute scores
    scores_1 = try(compute_scores(Y=Y, Y_mask=Y_mask, loadings=loadings_1, quantile.value=quantile.value), silent=FALSE)
    if(!is.matrix(scores_1))
    {
      warning('Computation of scores failed at iteration ', i, '. Providing results from previous iteration.')
      function.warnings$scores = TRUE
      break
    }

    # Compute objective function
    objective_function_1 = compute_objective_value(quantile.value=quantile.value, Y=Y, scores=scores_1, loadings=loadings_1)
    objective_function_array = c(objective_function_array, objective_function_1)

    # COMPROBATION OF CONVERGENCE ---------------------------------------------

    if(objective_function_1 < best_results$objective_function)
    {
      best_results = list(loadings=loadings_1, scores=scores_1, objective_function=objective_function_1, spline.coefficients=spline_coefficients_1, iteration=i)
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
      function.warnings$diverging.loop = TRUE
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
  if(i == n.iters){warning('Algorithm reached maximum number of iterations without convergence: ', n.iters, ' iterations')}

  # ROTATION TO ENSURE IDENTIFIABILITY -----------------------------------

  mean_scores = colMeans(matrix(best_results$scores[,-1], ncol=npc)) # Required in the prediction step
  normalization_result = try(rotate_scores_and_loadings(loadings=best_results$loadings, scores=best_results$scores), silent=TRUE)
  error_chr = class(normalization_result)
  if(error_chr=="try-error")
  {
    warning('Failed rotation step to ensure identifiability. Providing un-rotated results. You can try changing the method, or reducing the number of components')
    function.warnings$rotation = TRUE
    normalization_result = list(loadings=best_results$loadings, scores=best_results$scores, rotation.matrix=diag(npc))
  }

  # EXPLAINED VARIABILITY -----------------------------------------------------

  pve = compute_explained_variability(normalization_result$scores)

  global_end_time = base::Sys.time()
  global_execution_time = difftime(global_end_time, global_start_time, units='secs')

  results = new_fqpca(
    loadings=normalization_result$loadings,
    scores=normalization_result$scores,
    pve=pve,
    objective.function.value=best_results$objective_function,
    list.objective.function.values=objective_function_array,
    execution.time=global_execution_time,
    function.warnings=function.warnings,
    mean_scores_unnormalized = mean_scores,
    npc=npc,
    quantile.value=quantile.value,
    method=method,
    alpha.ridge=alpha.ridge,
    periodic=periodic,
    splines.df=splines.df,
    tol=tol,
    n.iters=best_results$iteration,
    verbose=verbose,
    seed=seed,
    rotation.matrix=normalization_result$rotation.matrix,
    spline.coefficients=best_results$spline.coefficients,
    spline.basis=spline.basis
  )
  return(results)
}

# PREDICTIONS -----------------------------------------------------------------

#' Predict FQPCA
#'
#' ## S3 method for class 'fqpca_object'
#' Given a new matrix Y, predicts the value of the scores associated to the given matrix
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
#' Y = matrix(rep(sin(seq(0, 2*pi, length.out=144)), 150), byrow=TRUE, nrow=150)
#' Y = Y + matrix(rnorm(150*144, 0, 0.4), nrow=150)
#'
#' # Add missing observations
#' Y[sample(150*144, as.integer(0.2*150*144))] = NA
#'
#' results = fqpca(Y=Y[1:100,], npc=2, quantile.value=0.5)
#'
#' predictions = predict(object=results, newdata=Y[101:150,])
predict.fqpca_object = function(object, newdata, ...)
{
  if (!inherits(object, "fqpca_object"))
  {
    stop('The object must be of class fqpca_object')
  }

  if(!is.data.frame(newdata) & !is.matrix(newdata))
  {
    # This structure allows tf to be a suggested package rather than mandatory
    if(!tf::is_tf(newdata))
    {
      stop('newdata is not of a valid type. Object provided: ', class(newdata))
    }
  }

  newdata = unname(as.matrix(newdata))
  if(ncol(newdata) == 1){ newdata = t(newdata)}

  Y_mask = !base::is.na(newdata)
  n_obs = base::nrow(newdata)
  n_time = base::ncol(newdata)
  if(sum(!Y_mask) == n_obs * n_time){stop('newdata contains no information. Check that there are values different from NA')}

  # Estimation of unnormalized loadings
  loadings_1 = compute_loadings(spline.basis=object$spline.basis, spline.coefficients=object$spline.coefficients)

  # Estimation of scores
  scores_1 = try(compute_scores(Y=newdata, Y_mask=Y_mask, loadings=loadings_1, quantile.value=object$quantile.value), silent=FALSE)
  if(!is.matrix(scores_1))
  {
    stop('Computation of scores failed.')
  }
  scores_1 = matrix(scores_1[,-1], ncol=object$npc) # Remove intercept column
  scores_1 = scores_1 - matrix(rep(object$mean_scores_unnormalized, n_obs), nrow=n_obs, byrow=T)
  scores_1 = scores_1 %*% t(base::solve(object$rotation.matrix)) # Normalize
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
#' Y = matrix(rep(sin(seq(0, 2*pi, length.out=144)), 150), byrow=TRUE, nrow=150)
#' Y = Y + matrix(rnorm(150*144, 0, 0.4), nrow=150)
#'
#' # Add missing observations
#' Y[sample(150*144, as.integer(0.2*150*144))] = NA
#'
#' results = fqpca(Y=Y, npc=2, quantile.value=0.5)
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
