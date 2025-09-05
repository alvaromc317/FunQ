
# SCORES ----------------------------------------------------------------------

#' @title Compute scores
#' @description Inner function to compute the scores of the fosqr-fqpca methodology.
#' @param Y The \eqn{(N \times T)} matrix of observed time instants.
#' @param Y.mask Mask matrix of the same dimensions as Y indicating which observations in Y are known.
#' @param Y.fosr.hat fosqr fitted values
#' @param fqpca.loadings Matrix of loading coefficients.
#' @param quantile.value The quantile considered.
#' @return The matrix of scores.
compute_scores_fosr_fqpca <- function(
    Y,
    Y.mask,
    Y.fosr.hat,
    fqpca.loadings,
    quantile.value)
{
  n.obs <- base::nrow(Y)
  npc <- base::ncol(fqpca.loadings)
  scores <-  base::matrix(0, n.obs, npc)
  for(i in base::seq(n.obs))
  {
    loadings.i <- fqpca.loadings[Y.mask[i, ], , drop=FALSE]
    Yi <- Y[i, Y.mask[i, ]] - Y.fosr.hat[i, Y.mask[i, ]]
    if(length(Yi) < ncol(loadings.i)){stop('Observation ', i, ' must have a larger number of timepoints than npcs.')}
    scores[i,] <-  quantreg::rq.fit.br(y=Yi, x=loadings.i, tau=quantile.value)$coefficients
  }
  return(scores)
}

# SPLINES ---------------------------------------------------------------------

#' @title Compute fosqr spline coefficients
#' @description Inner function to compute the spline coefficients of the fosqr-fqpca methodology
#' @param Y.vector vectorised version of Y, the \eqn{(N \times T)} matrix of observed time instants.
#' @param Y.mask Mask matrix of the same dimensions as Y indicating which observations in Y are known.
#' @param regressors Matrix of covariates.
#' @param spline.basis The spline basis matrix.
#' @param quantile.value The quantile considered.
#' @param method Method used in the resolution of the quantile regression model. It currently accepts the methods \code{c('conquer', 'quantreg')}.
#' @return The matrix of spline coefficients.
fosqrfqpca_compute_fosr_spline_coefficients <- function(
    Y.vector,
    Y.mask,
    regressors,
    spline.basis,
    quantile.value,
    method)
{
  n.obs <- base::nrow(regressors)
  n.regressors <- base::ncol(regressors)
  n.basis <- base::ncol(spline.basis)
  intercept.spline.basis <-  spline.basis[, -1]
  tensor.matrix <- base::matrix(0, nrow=base::length(Y.vector), ncol=n.regressors*n.basis-1)
  row.idx <- 1
  for(i in base::seq(n.obs))
  {
    regressors.i <- regressors[i, -1, drop=FALSE]
    tmp.splines <- spline.basis[Y.mask[i, ], , drop = FALSE]
    tmp.intercept.splines <- intercept.spline.basis[Y.mask[i, ], , drop = FALSE]
    n.time.i <- base::nrow(tmp.splines)
    tensor.matrix[row.idx:(row.idx + n.time.i - 1), ] <- base::cbind(tmp.intercept.splines, base::kronecker(regressors.i, tmp.splines))
    row.idx <- row.idx+n.time.i
  }
  if(method == 'conquer'){
    B.vector <- conquer::conquer(Y=Y.vector, X=tensor.matrix, tau=quantile.value)$coeff
  }else if(method == 'quantreg'){
    B.vector <- quantreg::rq(Y.vector~tensor.matrix, tau=quantile.value, method='fnb')$coefficients
  }
  spline.coefficients <- base::matrix(B.vector, byrow=FALSE, ncol=n.regressors)
  return(spline.coefficients)
}

#' @title Compute fosqr and fqpca spline coefficients
#' @description Inner function to compute the spline coefficients of the fosqr-fqpca methodology
#' @param Y.vector vectorised version of Y, the \eqn{(N \times T)} matrix of observed time instants.
#' @param Y.mask Mask matrix of the same dimensions as Y indicating which observations in Y are known.
#' @param fqpca.scores Matrix of fqpca scores
#' @param regressors Matrix of covariates.
#' @param spline.basis The spline basis matrix.
#' @param quantile.value The quantile considered.
#' @param method Method used in the resolution of the quantile regression model. It currently accepts the methods \code{c('conquer', 'quantreg')}.
#' @return The matrix of spline coefficients.
fosqrfqpca_compute_spline_coefficients <- function(
    Y.vector,
    Y.mask,
    fqpca.scores,
    regressors,
    spline.basis,
    quantile.value,
    method)
{
  n.obs <- base::nrow(Y.mask)
  npc <- base::ncol(fqpca.scores)
  n.regressors <- base::ncol(regressors)
  n.basis <- base::ncol(spline.basis)
  intercept.spline.basis <-  spline.basis[, -1]
  tensor.matrix <- base::matrix(0, nrow=base::length(Y.vector), ncol=(n.regressors+npc)*n.basis-1)
  row.idx <- 1
  for(i in base::seq(n.obs))
  {
    regressors.i <- regressors[i, -1, drop=FALSE]
    scores.i <- fqpca.scores[i, , drop=FALSE]
    tmp.splines <- spline.basis[Y.mask[i, ], , drop=FALSE]
    tmp.intercept.splines <- intercept.spline.basis[Y.mask[i, ], , drop = FALSE]
    n.time.i <- base::nrow(tmp.splines)
    tensor.matrix[row.idx:(row.idx + n.time.i - 1), ] <- base::cbind(tmp.intercept.splines, base::kronecker(matrix(c(regressors.i, scores.i), nrow=1), tmp.splines))
    row.idx <- row.idx+n.time.i
  }
  if(method == 'conquer'){
    B.vector <- conquer::conquer(Y=Y.vector, X=tensor.matrix, tau=quantile.value)$coeff
  }else if(method == 'quantreg'){
    B.vector <- quantreg::rq(Y.vector~tensor.matrix, tau=quantile.value, method='fnb')$coefficients
  }
  spline.coefficients <- base::matrix(B.vector, byrow=FALSE, ncol=(n.regressors+npc))
  fosr.spline.coef <- spline.coefficients[, 1:n.regressors, drop=FALSE]
  fqpca.spline.coef <- spline.coefficients[, (n.regressors+1):ncol(spline.coefficients), drop=FALSE]
  results <- list(fosr.spline.coef=fosr.spline.coef, fqpca.spline.coef=fqpca.spline.coef)
  return(results)
}

# OTHER FUNCTIONS -------------------------------------------------------------

# MAIN ------------------------------------------------------------------------

# PREDICTIONS -----------------------------------------------------------------

# BASIC PLOT ------------------------------------------------------------------
