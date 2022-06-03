# FUNCTIONAL QUANTILE PCA ALGORITHM -------------------------------------------

#' QR objective value
#'
#' Objective value function of a quantile regression model
#'
#' @param quantile_value The quantile considered. Default is the median 0.5.
#' @param x The vector.
#'
#' @return Objective value function for a quantile regression model.
quantile_function = function(x, quantile_value=0.5)
{
  # Quantile check function
  return(0.5 * base::abs(x) + (quantile_value - 0.5) * x)
}

#' Penalized quantile regression
#'
#' Uses the CVXR library to solve a penalized quantile regression model
#'
#' @param x N times p matrix of predictive variables.
#' @param y N times 1 vector of response.
#' @param R Quadratic matrix used to apply a ridge based penalty.
#' @param quantile_value The quantile considered. Default is the median 0.5.
#' @param lambda The hyperparameter controlling the penalization. Default is 1e-3.
#' @param verbose Boolean indicating verbosity of the function. Default is FALSE.
#' @param solver Solver to be used by the CVXR package in the resolution of the penalized quantile regression models. The speed of the methodology can vary depending on the solver used. There are free alternatives available like 'SCS', 'ECOS' or 'OSQP'. It is also possible to use licensed programs like 'MOSEK' or 'GUROBI'. Default is 'SCS'.
#'
#' @return Beta coefficients of the penalized quantile regression model.
qr_ridge = function(x, y, R, quantile_value=0.5, lambda=1e-3, verbose=FALSE, solver='SCS')
{
  # Solve a quantile regression model with a ridge type penalty
  n = dim(x)[1]
  m = dim(x)[2]
  beta_var = CVXR::Variable(m)
  objective_function = (1.0 / n) * base::sum(quantile_function(quantile_value=quantile_value, x=(y - x %*% beta_var)))
  ridge_penalization = lambda * CVXR::quad_form(beta_var, R)
  objective = CVXR::Minimize(objective_function + ridge_penalization)
  problem = CVXR::Problem(objective)
  solution = suppressWarnings(try(CVXR::solve(problem, solver=solver), silent=TRUE))
  error_chr = class(solution)
  if(error_chr=="try-error" & verbose==TRUE){message('Error with solver: ', solver, '. Use an alternative solver like ECOS, OSQP, or licensed programs like MOSEK or GUROBI')}
  beta_sol = c(solution$getValue(beta_var))
  results = beta_sol
  return(results)
}

# NORMALIZATION PROCESS -------------------------------------------------------

#' FQPCA normalization
#'
#' Performs the normalization of matrices of loadings and scores in order to ensure the solution is unique.
#'
#' @param Fhat Matrix of loadings
#' @param Lhat Matrix of scores
#'
#' @return The normalized matrices of loadings and scores and the normalization matrix.
algorithm_normalization = function(Fhat, Lhat)
{
  # Normalization process to ensure identifiability of solution
  n_time = nrow(Fhat)
  n_obs = nrow(Lhat)
  Fhat_tmp = Fhat[, 2:dim(Fhat)[2]]
  Lhat_tmp = Lhat[, 2:dim(Lhat)[2]]
  # Step 4: Normalize variables (code extracted from IQR.m file (original))
  sigmaF = (1/n_time) * t(Fhat_tmp) %*% Fhat_tmp
  sigmaA = (1/n_obs) * t(Lhat_tmp) %*% Lhat_tmp
  dum1 = expm::sqrtm(sigmaF) %*% sigmaA %*% expm::sqrtm(sigmaF)
  svd_decomp = base::svd(dum1)
  R = expm::sqrtm(base::solve(sigmaF)) %*% svd_decomp$u
  Fhat_tmp = Fhat_tmp %*% R;
  Lhat_tmp = Lhat_tmp %*% t(base::solve(R))
  Fhat = cbind(Fhat[,1], Fhat_tmp)
  Lhat = cbind(Lhat[,1], Lhat_tmp)
  results = list(Fhat=Fhat, Lhat=Lhat, normalization_matrix=R)
  return(results)
}

# INNER FUNCTIONS -------------------------------------------------------------

#' Compute loadings and scores
#'
#' Given matrix x and spline information, computes the matrices of loadings and scores.
#'
#' @param x The N by T matrix of observed time instants.
#' @param x_mask Mask matrix of the same dimensions asx indicating wihch observations in x are known.
#' @param quantile_value The quantile considered. Default is the median 0.5.
#' @param basis_coef The matrix of spline basis coefficients.
#' @param spline_basis The matrix of spline basis.
#' @param spline_basis_i A list containing, for each observation xi, the matrix of associated spline basis.
#'
#' @return The unnormalized matrices of loadings and scores, along with the mean and standard deviations of the scores (these last are used in the prediction function).
compute_F_Lambda = function(x, x_mask, quantile_value, basis_coef, spline_basis, spline_basis_i)
{
  # Given the spline basis and the basis coefficients, obtain F and Lambda matrices
  n_components = ncol(basis_coef)
  n_obs = nrow(x)
  # Obtain the final matrix of F factors based on full spline basis
  F1 = spline_basis %*% basis_coef
  # Given F1 compute Lambda1 based on the intercept - centered data
  Lambda1 =  matrix(0, n_obs, n_components-1)
  for(i in seq(n_obs))
  {
    Fi = spline_basis_i[[i]] %*% basis_coef
    xi = x[i, x_mask[i,]]
    xi = xi - Fi[,1]
    Lambda1[i,] =  quantreg::rq(xi ~ -1 + Fi[,2:(n_components)], tau=quantile_value)$coefficients
  }
  Lambda1 = base::scale(Lambda1)
  mean_Lambda1 = attr(Lambda1, 'scaled:center')
  sd_Lambda1 = attr(Lambda1, 'scaled:scale')
  Lambda1 = cbind(1, Lambda1)
  results = list(F1 = F1, Lambda1 =Lambda1, mean_Lambda1=mean_Lambda1, sd_Lambda1=sd_Lambda1)
  return(results)
}

#' Inner loop function of the fqpca methodology.
#'
#' @param x The N by T matrix of observed time instants.
#' @param x_mask Mask matrix of the same dimensions asx indicating wihch observations in x are known.
#' @param Lambda0 Initial matrix of scores.
#' @param spline_basis_i A list containing, for each observation xi, the matrix of associated spline basis.
#' @param quantile_value The quantile considered. Default is the median 0.5.
#' @param lambda_ridge The hyperparameter controlling the penalization on the splines. Default is 1e-12.
#' @param R_block Block diagonal matrix made of quadratic matrices used to apply a ridge based penalty.
#' @param solver Solver to be used by the CVXR package in the resolution of the penalized quantile regression models. The speed of the methodology can vary depending on the solver used. There are free alternatives available like 'SCS', 'ECOS' or 'OSQP'. It is also possible to use licensed programs like 'MOSEK' or 'GUROBI'. Default is 'SCS'.
#'
#' @return The loadings and scores matrices, the matrix of spline coefficients and the objective function value at the current iteration.
inner_loop = function(x, x_mask, Lambda0, spline_basis_i, quantile_value, lambda_ridge, R_block, solver)
{
  # At each iteration, call this function to obtain the estimation of F and Lambda
  n_obs = base::nrow(x)
  n_components = base::ncol(Lambda0)
  n_basis = base::ncol(spline_basis_i[[1]])
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
    lambda_i = base::matrix(Lambda0[i,], nrow=1)
    n_time_i = base::nrow(spline_basis_i[[i]])
    tensor_matrix[row_idx:(row_idx+n_time_i-1),] = base::kronecker(lambda_i, spline_basis_i[[i]])
    row_idx = row_idx+n_time_i
  }
  # Obtain spline coefficients and reshape as matrix
  B_vector = qr_ridge(x=tensor_matrix, y=x_vector,  R=R_block, quantile_value=quantile_value, lambda=lambda_ridge, solver=solver)
  B = base::matrix(B_vector, byrow=FALSE, ncol=n_components)

  # Obtain estimation of F
  F1 = list()
  for(i in base::seq(n_obs))
  {
    F1[[i]] = spline_basis_i[[i]] %*% B
  }

  # Given F1 compute Lambda1 based on the intercept - centered data
  Lambda1 =  base::matrix(0, n_obs, n_components-1)
  for(i in base::seq(n_obs))
  {
    xi = x[i, x_mask[i,]]
    xi = xi - F1[[i]][,1]
    Lambda1[i,] =  quantreg::rq(xi ~ -1 + F1[[i]][,2:(n_components)], tau=quantile_value)$coefficients
  }
  Lambda1 = base::scale(Lambda1)
  Lambda1 = base::cbind(1, Lambda1)

  # Given Lambda1 and F1 compute of convergence criteria
  of_value1 = 0
  for(i in base::seq(n_obs))
  {
    xi = x[i, x_mask[i,]]
    of_value1 = of_value1 + base::mean(quantile_function(quantile_value=quantile_value, x=(xi - Lambda1[i,] %*% t(F1[[i]]))))
  }
  of_value1 = 1/n_obs * of_value1
  result = list(F1=F1, Lambda1=Lambda1, B=B, of_value1=of_value1)
  return(result)
}

# MAIN ------------------------------------------------------------------------

new_fqpca <- function(loadings, scores, unnormalized_loadings, normalization_matrix,
                      spline_coefficients, mean_unnormalized_scores, sd_unnormalized_scores,
                      objective_function_value, list_objective_function_values, splines_df,
                      periodic, quantile_value, n_components, n_iters, lambda_ridge, spline_basis, execution_time,
                      error_checker_normalization, error_checker_loop, diverging_loop)
{
  structure(list(
    loadings = loadings,
    scores = scores,
    unnormalized_loadings = unnormalized_loadings,
    normalization_matrix = normalization_matrix,
    spline_coefficients = spline_coefficients,
    mean_unnormalized_scores = mean_unnormalized_scores,
    sd_unnormalized_scores = sd_unnormalized_scores,
    objective_function_value = objective_function_value,
    list_objective_function_values = list_objective_function_values,
    splines_df = splines_df,
    periodic = periodic,
    quantile_value = quantile_value,
    n_iters = n_iters,
    lambda_ridge=lambda_ridge,
    spline_basis = spline_basis,
    execution_time = execution_time,
    error_checker_normalization = error_checker_normalization,
    error_checker_loop = error_checker_loop,
    diverging_loop = diverging_loop
  ),
  class = "fqpca")
}


#' FQPCA
#'
#' Solves the functional quantile principal component analysis methodology
#'
#' @param x The N by T matrix of observed time instants
#' @param n_components The number of estimated components. Default is 2.
#' @param quantile_value The quantile considered. Default is the median 0.5.
#' @param lambda_ridge The hyperparameter controlling the penalization on the splines. Default is 1e-12.
#' @param periodic Boolean indicating if the data is expected to be periodic (start coincides with end) or not. Default is TRUE.
#' @param splines_df Degrees of freedom for the splines. It is recommended not to modify this parameter and control the smoothness via the hyperparameter lambda_ridge. Default is 30.
#' @param tol Tolerance on the convergence of the algorithm. Smaller values can spped up computation but may affect the quality of the estimations. Default is 1e-3.
#' @param n_iters Maximum number of iterations. Default is 30.
#' @param solver Solver to be used by the CVXR package in the resolution of the penalized quantile regression models. The speed of the methodology can vary depending on the solver used. There are free alternatives available like 'SCS', 'ECOS' or 'OSQP'. It is also possible to use licensed programs like 'MOSEK' or 'GUROBI'. Default is 'SCS'.
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
#' # Add missing observations
#' x[sample(200*144, as.integer(0.2*200*144))] = NA
#'
#' results = fqpca(x=x, n_components=1, quantile_value=0.5, lambda_ridge=1e-12)
#'
#' loadings = results$loadings
#' scores = results$scores
fqpca = function(x, n_components=2,  quantile_value=0.5, lambda_ridge=1e-12,  periodic=TRUE, splines_df=30, tol=1e-3, n_iters=30, solver='SCS', verbose=FALSE, seed=NULL)
{
  global_start_time = base::Sys.time()
  if(!base::is.null(seed)){base::set.seed(seed)}
  x = base::as.matrix(x)
  # Step 1: Initialize values
  n_obs = base::dim(x)[1]
  n_time = base::dim(x)[2]
  x_axis = base::seq(0, 1, length.out = n_time)
  x_mask = !base::is.na(x)
  # Initialize splines and second derivative ---
  if(periodic)
  {
    knots = base::seq(0, 1, length.out = 2 + splines_df - 1)  # 2 boundary knots + df - intercept
  } else{
    knots = base::seq(0, 1, length.out = 2 + splines_df - 3 - 1) # 2 boundary knots + df - degree - intercept
  }
  knots = knots[2:(base::length(knots)-1)]
  spline_basis_i = list()
  for(i in base::seq(n_obs))
  {
    xi_axis = x_axis[x_mask[i,]]
    spline_basis_i[[i]] = splines2::mSpline(xi_axis, degree=3, intercept=TRUE, knots=knots, periodic=periodic, Boundary.knots=c(min(x_axis), max(x_axis)))
  }
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
  if(!all(base::eigen(R_block)$values > 0)) {message('Splines penalization matrix is not positive definite')}

  # Generate F0 and compute Lambda0 based on centered data ---
  basis_coef = base::matrix(stats::rnorm(base::dim(spline_basis)[2] * (n_components+1), mean=0, sd=1/n_time),
                            nrow=base::dim(spline_basis)[2], ncol=n_components+1)
  Lambda0 = base::matrix(0, n_obs, n_components)
  F0 = list()
  for(i in base::seq(n_obs))
  {
    xi = x[i, x_mask[i,]]
    F0[[i]] = spline_basis_i[[i]] %*% basis_coef
    xi = xi - F0[[i]][,1]
    Lambda0[i,] =  quantreg::rq(xi ~ -1 + F0[[i]][,2:(n_components+1)], tau=quantile_value)$coefficients
  }
  # Center and scale Lambda. Add a 1 as an intercept column
  Lambda0 = base::scale(Lambda0)
  Lambda0 = base::cbind(1, Lambda0)
  # Given Lambda0 and F0 compute convergence criteria (of==objective function)
  of_value0 = 0
  for(i in base::seq(n_obs))
  {
    xi = x[i, x_mask[i,]]
    of_value0 = of_value0 + base::mean(quantile_function(quantile_value=quantile_value, x=(xi - Lambda0[i,] %*% t(F0[[i]]))))
  }
  of_value0 = 1/n_obs * of_value0
  of_value = c(of_value0)

  # Initialize solution
  best_iteration = 0
  best_F = F0
  best_Lambda = Lambda0
  best_of = of_value0
  best_B = basis_coef
  loop_result_backup = list(F1=F0, Lambda1=Lambda0, B=basis_coef, of_value1=of_value0)

  # Iterate until convergence or until max_iters is reached
  convergence_criteria = c()
  error_checker_loop = FALSE
  error_checker_normalization = FALSE
  diverging_loop = FALSE
  for(i in 1:n_iters)
  {
    loop_start_time = base::Sys.time()
    loop_result = try(inner_loop(x=x, x_mask=x_mask, Lambda0=Lambda0, spline_basis_i=spline_basis_i, quantile_value=quantile_value, lambda_ridge=lambda_ridge, R_block=R_block, solver=solver), silent=TRUE)
    error_chr = class(loop_result)
    if(error_chr=="try-error")
    {
      # If an error occurred (aka, algorithm failed to converge) exit the loop
      message('Inner loop error')
      error_checker_loop = TRUE
      loop_result = loop_result_backup
      break
    } else {
      loop_result_backup = loop_result
    }
    # Extract results from inner loop function
    F1 = loop_result$F1
    Lambda1 = loop_result$Lambda1
    of_value1 = loop_result$of_value1
    of_value = c(of_value, of_value1)
    # Objective is to minimize objective value. Pick the solution associated to the minimum objective value function.
    if(of_value1 < best_of)
    {
      best_of = of_value1
      best_iteration = i
      best_F = F1
      best_Lambda = Lambda1
      best_B = loop_result$B
    }
    convergence_criteria[i] = base::abs(of_value1 - of_value0)
    # Measure computation time
    loop_end_time = base::Sys.time()
    loop_execution_time = loop_end_time - loop_start_time
    if(verbose)
    {
      message('Iteration ', i, ' completed in ', base::round(loop_execution_time, 3), ' seconds',
              '\n', 'Convergence criteria value: ', base::round(convergence_criteria[i], 4),
              '\n', 'OF0 value: ', base::round(of_value0, 4), ' OF1 value: ', base::round(of_value1, 4))
      message('___________________________________________________')
    }

    if(convergence_criteria[i] < tol)
    {
      if(verbose){message('Convergence criteria is satisfied. Value: ', base::round(convergence_criteria[i], 4))}
      break
    }
    # If the objective value function (that should be decreasing) is
    # 50 times larger than the original objective value function, break
    if( of_value[i+1] > 50 * of_value[2])
    {
      diverging_loop = TRUE
      if(verbose)
      {
        message('Breaking diverging loop at iteration ', i,
                '\n', 'Objective function value: ', base::round(of_value[i], 3),
                '\n', 'Initial objective value function: ', base::round(of_value[2], 3))
        message('Results of the iteration with the smallest objective value: ', base::round(best_of, 3), ' will be provided')
      }
      break
    }
    # Update values
    of_value0 = of_value1
    Lambda0 = Lambda1
    F0 = F1
  }
  if(i == n_iters){ if(verbose) {message('Algorithm reached max_iters: ', n_iters, ' iterations')}}

  best_results = compute_F_Lambda(x=x, x_mask=x_mask, quantile_value=quantile_value, basis_coef=best_B, spline_basis=spline_basis, spline_basis_i=spline_basis_i)
  # Normalization of best results
  normalization_result = try(algorithm_normalization(Fhat=best_results$F1, Lhat=best_results$Lambda1), silent=TRUE)
  error_chr = class(normalization_result)
  if(error_chr=="try-error")
  {
    message('Failed normalization of results. Providing unnormalized results')
    error_checker_normalization = TRUE
    normalization_result = list(Fhat=best_F, Lhat=best_Lambda, normalization_matrix=diag(n_components))
  }
  best_Fhat = normalization_result$Fhat
  best_Lhat = normalization_result$Lhat
  best_normalization_matrix = normalization_result$normalization_matrix

  global_end_time = base::Sys.time()
  global_execution_time = global_end_time - global_start_time

  results = new_fqpca(
    loadings = best_Fhat,
    scores = best_Lhat,
    unnormalized_loadings = best_results$F1,
    normalization_matrix = best_normalization_matrix,
    spline_coefficients = best_B,
    mean_unnormalized_scores = best_results$mean_Lambda1,
    sd_unnormalized_scores = best_results$sd_Lambda1,
    objective_function_value = best_of,
    list_objective_function_values = of_value,
    splines_df=splines_df,
    periodic = periodic,
    quantile_value = quantile_value,
    n_components = n_components,
    n_iters = best_iteration,
    lambda_ridge=lambda_ridge,
    spline_basis=spline_basis,
    execution_time=global_execution_time,
    error_checker_normalization=error_checker_normalization,
    error_checker_loop=error_checker_loop,
    diverging_loop=diverging_loop
  )

  return(results)
}

# PREDICTIONS -----------------------------------------------------------------

#' Predict FQPCA
#'
#' ## S3 method for class 'fqpca'
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
#' # Generate fake dataset with 200 observations and 144 time points
#'
#' x = matrix(rep(sin(seq(0, 2*pi, length.out=144)), 200), byrow=TRUE, nrow=200)
#' x = x + matrix(rnorm(200*144, 0, 0.4), nrow=200)
#'
#' # Add missing observations
#' x[sample(200*144, as.integer(0.2*200*144))] = NA
#'
#' results = fqpca(x=x[1:150,], n_components=1, quantile_value=0.5)
#'
#' predictions = predict(results, newdata=x[151:200,])
predict.fqpca = function(object, newdata, ...)
{
  if (!inherits(object, "fqpca"))
  {
    stop('The object must be of class fqpca')
  }
  # Unload required information from fqfa object
  spline_coefficients=object$spline_coefficients
  normalization_matrix=object$normalization_matrix
  mean_unnormalized_scores=object$mean_unnormalized_scores
  sd_unnormalized_scores=object$sd_unnormalized_scores
  splines_df =object$splines_df
  periodic = object$periodic
  quantile_value = object$quantile_value
  x = base::as.matrix(newdata)
  n_obs = base::nrow(x)
  n_time = base::ncol(x)
  n_components = base::ncol(spline_coefficients)
  x_axis = base::seq(0, 1, length.out = n_time)
  x_mask = !base::is.na(x)
  mean_Lhat_matrix = base::matrix(base::rep(mean_unnormalized_scores, n_obs), nrow=n_obs, byrow=TRUE)
  sd_Lhat_matrix = base::matrix(base::rep(sd_unnormalized_scores, n_obs), nrow=n_obs, byrow=TRUE)
  if(periodic)
  {
    knots = base::seq(0, 1, length.out = 2 + splines_df - 1)  # 2 boundary knots + df - intercept
  } else{
    knots = base::seq(0, 1, length.out = 2 + splines_df - 3 - 1) # 2 boundary knots + df - degree - intercept
  }
  knots = knots[2:(length(knots)-1)]
  spline_basis_i = list()
  for(i in base::seq(n_obs))
  {
    xi_axis = x_axis[x_mask[i,]]
    spline_basis_i[[i]] = splines2::mSpline(xi_axis, degree=3, intercept=TRUE, knots=knots, periodic=periodic, Boundary.knots=c(min(x_axis), max(x_axis)))
  }

  # Given F1 compute Lambda1 based on the intercept - centered data
  Lambda1 =  base::matrix(0, n_obs, n_components-1)
  for(i in base::seq(n_obs))
  {
    Fi = spline_basis_i[[i]] %*% spline_coefficients
    xi = x[i, x_mask[i,]]
    xi = xi - Fi[,1]
    Lambda1[i,] =  quantreg::rq(xi ~ -1 + Fi[,2:(n_components)], tau=quantile_value)$coefficients
  }
  Lambda1 = (Lambda1 - mean_Lhat_matrix) / sd_Lhat_matrix
  Lambda1 = Lambda1 %*% t(base::solve(normalization_matrix))
  Lambda1 = base::cbind(1, Lambda1)
  return(Lambda1)
}
