
# QUANTILE BASED FUNCTIONS ----------------------------------------------------

#' @title Quantile check (pinball) loss
#' @description Evaluates the quantile check loss elementwise.
#' @param x The vector.
#' @param quantile.value The quantile considered. Defaults to 0.5.
#' @return Objective value function for a quantile regression model.
check_loss <- function(x, quantile.value = 0.5)
{
  return(0.5 * abs(x) + (quantile.value - 0.5) * x)
}

#' @title Quantile error computation
#' @description Quantile error metric computation
#' @param Y \eqn{(N \times T)} matrix of observed time instants.
#' @param Y.pred \eqn{(N \times T)} matrix of predicted time instants.
#' @param quantile.value The quantile considered.
#' @return The quantile error between the two matrices
#' @export
quantile_error <- function(Y, Y.pred, quantile.value)
{
  mean(check_loss(quantile.value = quantile.value, x = (Y - Y.pred)), na.rm = TRUE)
}

#' @title Proportion of points under quantile
#' @description Proportion of points under quantile.
#' @param Y \eqn{(N \times T)} matrix of observed time instants.
#' @param Y.pred \eqn{(N \times T)} matrix of predicted time instants.
#' @return The proportion of points under the estimated quantile function
#' @export
proportion_under_quantile <- function(Y, Y.pred)
{
  mean(Y <= Y.pred, na.rm = TRUE)
}

#' @title Compute objective value function.
#' @description Inner function to compute the objective value of the alternating quantile algorithms at each iteration.
#' @param Y \eqn{(N \times T)} matrix of observed time instants.
#' @param quantile.value The quantile considered.
#' @param scores The matrix of estimated scores.
#' @param intercept population intercept
#' @param loadings The matrix of estimated loadings.
#' @return The objective value function.
compute_objective_value <- function(Y, quantile.value, scores, intercept, loadings)
{
  Y.pred <- sweep(
    scores %*% t(loadings),
    MARGIN = 2,
    STATS = intercept,
    FUN = "+")
  objective.value <- quantile_error(Y=Y, Y.pred=Y.pred, quantile.value=quantile.value)
  return(objective.value)
}
