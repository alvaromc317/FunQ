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

#' Quantile regression model
#'
#' Uses the CVXR library to solve a quantile regression model
#'
#' @param x N times p matrix of predictive variables.
#' @param y N times 1 vector of response.
#' @param quantile_value The quantile considered. Default is the median 0.5.
#' @param solver Solver to be used by the CVXR package in the resolution of the penalized quantile regression models. The speed of the methodology can vary depending on the solver used. There are free alternatives available like 'SCS', 'ECOS' or 'OSQP'. It is also possible to use licensed programs like 'MOSEK' or 'GUROBI'. Default is 'SCS'.
#'
#' @return Beta coefficients of the penalized quantile regression model.
qr = function(x, y, quantile_value=0.5, solver='SCS')
{
  # Solve a quantile regression model with a ridge type penalty
  n = dim(x)[1]
  m = dim(x)[2]
  beta_var = CVXR::Variable(m)
  objective_function = (1.0 / n) * base::sum(quantile_function(quantile_value=quantile_value, x=(y - x %*% beta_var)))
  objective = CVXR::Minimize(objective_function)
  problem = CVXR::Problem(objective)
  solution = CVXR::solve(problem, solver=solver)
  beta_sol = c(solution$getValue(beta_var))
  results = beta_sol
  return(results)
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
#' @param solver Solver to be used by the CVXR package in the resolution of the penalized quantile regression models. The speed of the methodology can vary depending on the solver used. There are free alternatives available like 'SCS', 'ECOS' or 'OSQP'. It is also possible to use licensed programs like 'MOSEK' or 'GUROBI'. Default is 'SCS'.
#'
#' @return Beta coefficients of the penalized quantile regression model.
qr_ridge = function(x, y, R, quantile_value=0.5, lambda=1e-3, solver='SCS')
{
  # Solve a quantile regression model with a ridge type penalty
  n = dim(x)[1]
  m = dim(x)[2]
  beta_var = CVXR::Variable(m)
  objective_function = (1.0 / n) * base::sum(quantile_function(quantile_value=quantile_value, x=(y - x %*% beta_var)))
  ridge_penalization = lambda * CVXR::quad_form(beta_var, R)
  objective = CVXR::Minimize(objective_function + ridge_penalization)
  problem = CVXR::Problem(objective)
  solution = CVXR::solve(problem, solver=solver)
  beta_sol = c(solution$getValue(beta_var))
  results = beta_sol
  return(results)
}

