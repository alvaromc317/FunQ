
#' QR objective value
#'
#' Objective value function of a quantile regression model
#'
#' @param quantile.value The quantile considered. Default is the median 0.5.
#' @param x The vector.
#'
#' @return Objective value function for a quantile regression model.
quantile_function = function(x, quantile.value=0.5)
{
  # Quantile check function
  return(0.5 * base::abs(x) + (quantile.value - 0.5) * x)
}

#' Penalized quantile regression
#'
#' Uses the CVXR library to solve a penalized quantile regression model
#'
#' @param x N times p matrix of predictive variables.
#' @param y N times 1 vector of response.
#' @param R Quadratic matrix used to apply a ridge based penalty.
#' @param quantile.value The quantile considered. Default is the median 0.5.
#' @param lambda The hyperparameter controlling the penalization. Default is 1e-3.
#' @param verbose Boolean indicating verbosity of the function. Default is FALSE.
#' @param solver Solver to be used by the CVXR package in the resolution of the penalized quantile regression models. The speed of the methodology can vary depending on the solver used. There are free alternatives available like 'SCS', 'ECOS' or 'OSQP'. It is also possible to use licensed programs like 'MOSEK' or 'GUROBI'. Default is 'SCS'.
#'
#' @return Beta coefficients of the penalized quantile regression model.
qr_ridge = function(x, y, R, quantile.value=0.5, lambda=1e-3, verbose=FALSE, solver='SCS')
{
  # Solve a quantile regression model with a ridge type penalty
  n = dim(x)[1]
  m = dim(x)[2]
  beta_var = CVXR::Variable(m)
  objective_function = (1.0 / n) * base::sum(quantile_function(quantile.value=quantile.value, x=(y - x %*% beta_var)))
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

#' Compute loadings and scores
#'
#' Given matrix x and spline information, computes the matrices of loadings and scores.
#'
#' @param x The N by T matrix of observed time instants.
#' @param x_mask Mask matrix of the same dimensions asx indicating wihch observations in x are known.
#' @param quantile.value The quantile considered. Default is the median 0.5.
#' @param basis_coef The matrix of spline basis coefficients.
#' @param spline_basis The matrix of spline basis.
#' @param spline_basis_i A list containing, for each observation xi, the matrix of associated spline basis.
#'
#' @return The unnormalized matrices of loadings and scores, along with the mean and standard deviations of the scores (these last are used in the prediction function).
compute_F_Lambda = function(x, x_mask, quantile.value, basis_coef, spline_basis, spline_basis_i)
{
  # Given the spline basis and the basis coefficients, obtain F and Lambda matrices
  n.components = ncol(basis_coef)
  n_obs = nrow(x)
  # Obtain the final matrix of F factors based on full spline basis
  F1 = spline_basis %*% basis_coef
  # Given F1 compute Lambda1 based on the intercept - centered data
  Lambda1 =  matrix(0, n_obs, n.components-1)
  for(i in seq(n_obs))
  {
    Fi = spline_basis_i[[i]] %*% basis_coef
    xi = x[i, x_mask[i,]]
    xi = xi - Fi[,1]
    Lambda1[i,] =  quantreg::rq(xi ~ -1 + Fi[,2:(n.components)], tau=quantile.value)$coefficients
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
#' @param quantile.value The quantile considered. Default is the median 0.5.
#' @param lambda.ridge The hyperparameter controlling the penalization on the splines. Default is 1e-12.
#' @param R_block Block diagonal matrix made of quadratic matrices used to apply a ridge based penalty.
#' @param solver Solver to be used by the CVXR package in the resolution of the penalized quantile regression models. The speed of the methodology can vary depending on the solver used. There are free alternatives available like 'SCS', 'ECOS' or 'OSQP'. It is also possible to use licensed programs like 'MOSEK' or 'GUROBI'. Default is 'SCS'.
#'
#' @return The loadings and scores matrices, the matrix of spline coefficients and the objective function value at the current iteration.
inner_loop = function(x, x_mask, Lambda0, spline_basis_i, quantile.value, lambda.ridge, R_block, solver)
{
  # At each iteration, call this function to obtain the estimation of F and Lambda
  n_obs = base::nrow(x)
  n.components = base::ncol(Lambda0)
  n_basis = base::ncol(spline_basis_i[[1]])
  # vectorize X [x11, ..., x1T, x21, ..., x2T,...XNT]
  x_vector = c()
  for(i in base::seq(n_obs))
  {
    xi = x[i, x_mask[i,]]
    x_vector = c(x_vector, xi)
  }
  # Compute tensor product of Lambda and spline basis
  tensor_matrix = base::matrix(0, nrow=base::length(x_vector), ncol=n.components * n_basis)
  row_idx = 1
  for(i in base::seq(n_obs))
  {
    lambda_i = base::matrix(Lambda0[i,], nrow=1)
    n_time_i = base::nrow(spline_basis_i[[i]])
    tensor_matrix[row_idx:(row_idx+n_time_i-1),] = base::kronecker(lambda_i, spline_basis_i[[i]])
    row_idx = row_idx+n_time_i
  }
  # Obtain spline coefficients and reshape as matrix
  B_vector = qr_ridge(x=tensor_matrix, y=x_vector,  R=R_block, quantile.value=quantile.value, lambda=lambda.ridge, solver=solver)
  B = base::matrix(B_vector, byrow=FALSE, ncol=n.components)

  # Obtain estimation of F
  F1 = list()
  for(i in base::seq(n_obs))
  {
    F1[[i]] = spline_basis_i[[i]] %*% B
  }

  # Given F1 compute Lambda1 based on the intercept - centered data
  Lambda1 =  base::matrix(0, n_obs, n.components-1)
  for(i in base::seq(n_obs))
  {
    xi = x[i, x_mask[i,]]
    xi = xi - F1[[i]][,1]
    Lambda1[i,] =  quantreg::rq(xi ~ -1 + F1[[i]][,2:(n.components)], tau=quantile.value)$coefficients
  }
  Lambda1 = base::scale(Lambda1)
  Lambda1 = base::cbind(1, Lambda1)

  # Given Lambda1 and F1 compute of convergence criteria
  of_value1 = 0
  for(i in base::seq(n_obs))
  {
    xi = x[i, x_mask[i,]]
    of_value1 = of_value1 + base::mean(quantile_function(quantile.value=quantile.value, x=(xi - Lambda1[i,] %*% t(F1[[i]]))))
  }
  of_value1 = 1/n_obs * of_value1
  result = list(F1=F1, Lambda1=Lambda1, B=B, of_value1=of_value1)
  return(result)
}
