
# SCORES -----------------------------------------------------------------------

#' @title Compute scores
#' @description Inner function to compute the per-subject scores via independent quantile regressions on the loading functions.
#' @param Y The \eqn{(N \times T)} matrix of observed time instants.
#' @param Y.mask Mask matrix of the same dimensions as Y indicating which observations in Y are known.
#' @param loadings Matrix of loading coefficients.
#' @param quantile.value The quantile considered.
#' @param offset Optional offset subtracted from Y before fitting: NULL (no offset), a length-T vector (e.g. a population intercept), or an \eqn{(N \times T)} matrix (e.g. fitted values from a regression part).
#' @return The matrix of scores.
compute_scores <- function(
    Y,
    Y.mask,
    loadings,
    quantile.value,
    offset = NULL)
{
  # Initialize the matrix of scores
  n.obs <- nrow(Y)
  npc <- ncol(loadings)
  # A per-observation offset must match the dimensions of Y; anything else
  # (including a T x 1 matrix intercept) is treated as a length-T vector
  offset.per.obs <- is.matrix(offset) && nrow(offset) == nrow(Y) && ncol(offset) == ncol(Y)
  scores <-  matrix(0, n.obs, npc)
  for(i in seq_len(n.obs))
  {
    mask.i <- Y.mask[i, ]
    # Obtain the loadings associated to observation i with a specific time grid
    loadings.obs <- loadings[mask.i, , drop = FALSE]
    # Obtain the observation removing missing data and subtract the offset
    Yi <- Y[i, mask.i]
    if(!is.null(offset))
    {
      Yi <- Yi - if(offset.per.obs){offset[i, mask.i]}else{offset[mask.i]}
    }
    if(length(Yi) < ncol(loadings.obs)){stop('Observation ', i, ' must have a larger number of timepoints than npcs.')}
    scores[i,] <-  quantreg::rq.fit.br(y=Yi, x=loadings.obs, tau=quantile.value)$coefficients
  }
  return(scores)
}

# SPLINE COEFFICIENTS ----------------------------------------------------------

#' @title Extract the coefficient vector from a quantile regression model object
#' @description Returns the flat coefficient vector regardless of whether the model was fitted with \code{conquer::conquer} (\code{$coeff}) or \code{quantreg::rq} (\code{$coefficients}).
#' @param model A fitted model object.
#' @return The coefficient vector.
extract_coefficients <- function(model)
{
  if(!is.null(model$coeff)){return(model$coeff)}
  if(!is.null(model$coefficients)){return(model$coefficients)}
  stop("Model object contains no 'coeff' or 'coefficients' element.")
}

#' @title Fit the tensor-design quantile regression
#' @description Inner function to fit the spline coefficients of the alternating quantile algorithms: one large quantile regression of the vectorised data on the tensor-product design built from the current scores and the spline basis.
#' @param Y.vector vectorised version of Y, the \eqn{(N \times T)} matrix of observed time instants.
#' @param Y.mask Mask matrix of the same dimensions as Y indicating which observations in Y are known.
#' @param scores Matrix of scores (and / or regressors) scaling the spline blocks.
#' @param spline.basis The spline basis matrix.
#' @param quantile.value The quantile considered.
#' @param method Method used in the resolution of the quantile regression model. It currently accepts the methods \code{c('conquer', 'quantreg')}.
#' @param intercept Boolean. If TRUE (default) the model contains a functional intercept; if FALSE the fit has no intercept of any kind (used by the mfqpca within level).
#' @param penalized Boolean indicating if the smoothness should be controlled using a second derivative penalty (conquer only).
#' @param lambda.ridge Hyper parameter controlling the penalization on the second derivative of the splines. Only used when \code{penalized=TRUE}.
#' @param return.model Boolean. If TRUE, returns a list with the fitted model object, the spline coefficient matrix, the flat coefficient vector and the dense tensor design matrix (needed for the fosqr variance estimation).
#' @return The matrix of spline coefficients (with \code{npc + 1} columns when \code{intercept=TRUE}, \code{npc} otherwise), or the extended list when \code{return.model=TRUE}.
fit_tensor_quantile_regression <- function(
    Y.vector,
    Y.mask,
    scores,
    spline.basis,
    quantile.value,
    method,
    intercept = TRUE,
    penalized = FALSE,
    lambda.ridge = 0,
    return.model = FALSE)
{
  npc <- ncol(scores)
  n.basis <- ncol(spline.basis)

  # No-intercept fit (mfqpca within level)
  if(!intercept)
  {
    tensor.matrix <- build_tensor_matrix(
      scores=scores,
      Y.mask=Y.mask,
      spline.basis=spline.basis,
      intercept=FALSE)
    if(method == 'conquer'){
      # conquer2 is a fork of conquer that allows fitting the model without an intercept
      B.vector <- conquer2::conquer(X=tensor.matrix, Y=Y.vector, tau=quantile.value)$coeff
    }else if(method == 'quantreg'){
      B.vector <- quantreg::rq.fit.fnb(y=Y.vector, x=tensor.matrix, tau=quantile.value)$coefficients
    }
    return(matrix(B.vector, byrow=FALSE, ncol=npc))
  }

  # Penalized fit (conquer only; callers enforce method='conquer')
  if(penalized)
  {
    tensor.matrix <- build_tensor_matrix(scores=scores, Y.mask=Y.mask, spline.basis=spline.basis)
    if(lambda.ridge == 0)
    {
      B.vector <- conquer::conquer(Y=Y.vector, X=tensor.matrix, tau=quantile.value)$coeff
    }else{
      B.vector <- conquer::conquer.reg(
        Y=Y.vector,
        X=tensor.matrix,
        tau=quantile.value,
        penalty='elastic',
        lambda=lambda.ridge,
        para.elastic=0)$coeff
    }
    return(matrix(B.vector, byrow=FALSE, ncol=(npc+1)))
  }

  # Model-returning fit (fosqr family: the model object and the dense design
  # matrix are needed downstream for the sandwich variance estimation)
  if(return.model)
  {
    tensor.matrix <- build_tensor_matrix(scores=scores, Y.mask=Y.mask, spline.basis=spline.basis)
    if(method == 'conquer'){
      model <- conquer::conquer(Y=Y.vector, X=tensor.matrix, tau=quantile.value)
    }else if(method == 'quantreg'){
      model <- quantreg::rq(Y.vector~tensor.matrix, tau=quantile.value, method='fnb')
    }
    B.vector <- extract_coefficients(model)
    spline.coefficients <- matrix(B.vector, byrow=FALSE, ncol=(npc+1))
    return(list(
      model=model,
      spline.coefficients=spline.coefficients,
      B.vector=B.vector,
      tensor.matrix=tensor.matrix))
  }

  # Standard intercept model. For large problems the quantreg path uses a sparse
  # design and the sparse Frisch-Newton solver, which is substantially faster;
  # for small designs the dense solver has lower fixed overhead
  # (see additional_32_33_tests.md).
  use.sparse <- method == 'quantreg' && length(Y.vector) * ((npc+1) * n.basis) > 5e6
  if(use.sparse){
    sparse.design <- build_tensor_matrix_sparse(scores=scores, Y.mask=Y.mask, spline.basis=spline.basis)
    B.vector <- quantreg::rq.fit.sfn(sparse.design, Y.vector, tau=quantile.value)$coef
  }else{
    tensor.matrix <- build_tensor_matrix(scores=scores, Y.mask=Y.mask, spline.basis=spline.basis)
    if(method == 'conquer'){
      B.vector <- conquer::conquer(Y=Y.vector, X=tensor.matrix, tau=quantile.value)$coeff
    }else if(method == 'quantreg'){
      B.vector <- quantreg::rq(Y.vector~tensor.matrix, tau=quantile.value, method='fnb')$coefficients
    }
  }
  spline.coefficients <- matrix(B.vector, byrow=FALSE, ncol=(npc+1))
  return(spline.coefficients)
}

# LOADINGS ----------------------------------------------------------------------

#' @title Compute loadings (aka principal components)
#' @description Inner function to compute the intercept and loading functions from the spline coefficient matrix.
#' @param spline.basis The spline basis matrix.
#' @param spline.coefficients the matrix of spline coefficients.
#' @return The matrix of loadings
compute_loadings <- function(
    spline.basis,
    spline.coefficients)
{
  intercept.spline.basis <-  spline.basis[, -1]
  intercept.part <- c(cbind(1, intercept.spline.basis) %*% matrix(spline.coefficients[,1], ncol=1))
  fqpc.part <- spline.basis %*% spline.coefficients[, -1, drop=FALSE]
  return(list(intercept=intercept.part, loadings=fqpc.part))
}

# ROTATION ----------------------------------------------------------------------

#' @title Rotation of loadings and scores
#' @description Performs the SVD-based rotation of loadings and scores (with sign alignment) in order to ensure the solution is unique.
#' @param scores Matrix of scores.
#' @param loadings Matrix of loadings.
#' @param intercept Optional population intercept. When provided, the scores are mean-centered and the effect of their mean is moved into the intercept; when NULL (e.g. the mfqpca within level) the scores are rotated as they are.
#' @return The rotated matrices of loadings and scores, the (possibly updated) intercept and the rotation matrix.
rotate_factors <- function(
    scores,
    loadings,
    intercept = NULL)
{
  npc <- ncol(loadings)

  if(!is.null(intercept))
  {
    # Ensure scores are mean-centered and move the effect of the mean scores to the intercept
    means.scores <- matrix(colMeans(scores), ncol=1)
    scores <- scale(scores, center=TRUE, scale=FALSE)
    intercept <- intercept + loadings %*% means.scores
  }

  # Make covariance matrix positive definite using eigen-decomposition
  cov.score <- stats::cov(scores)
  eig <- eigen(cov.score, symmetric=TRUE)
  eig$values[eig$values < 1e-8] <- 1e-8
  if (npc == 1) {
    eig_vals_diag <- eig$values
  } else {
    eig_vals_diag <- diag(eig$values)
  }

  cov.score.pd <- tcrossprod(eig$vectors %*% eig_vals_diag, eig$vectors)
  chol_cov <- chol(cov.score.pd)
  svd.decomp <- svd(tcrossprod(chol_cov, loadings))
  rotation.matrix <- backsolve(chol_cov, sweep(
      svd.decomp$u,
      2,
      svd.decomp$d,
      "*"))
  scores <- scores %*% rotation.matrix
  loadings <- svd.decomp$v

  # Sign alignment: make the peak of each loading function positive so the solution is unique
  max.idx <- apply(abs(loadings), 2, which.max)
  peak.vals <- loadings[cbind(max.idx, seq_len(ncol(loadings)))]
  peak.vals[peak.vals==0] <- 1
  signs <- sign(peak.vals)
  loadings  <- sweep(
    loadings,
    2,
    signs,
    "*")
  scores    <- sweep(
    scores,
    2,
    signs,
    "*")
  rotation.matrix <- sweep(
    rotation.matrix,
    2,
    signs,
    "*")

  results <- list(
    scores=scores,
    intercept=intercept,
    loadings=loadings,
    rotation.matrix=rotation.matrix)
  return(results)
}

# RECONSTRUCTION AND COMPONENT SELECTION ----------------------------------------

#' @title Truncated reconstruction from scores and loadings
#' @description Helper function to compute fitted values from the first \code{n.comp} components.
#' @param scores matrix of scores
#' @param loadings matrix of loadings
#' @param n.comp number of components used. If 0 then a matrix of 0s is returned.
reconstruct_scores_loadings <- function(
    scores,
    loadings,
    n.comp)
{
  if (n.comp > 0)
  {
    r = scores[, seq_len(n.comp), drop = FALSE] %*% t(loadings[, seq_len(n.comp), drop = FALSE])
  }else{
    r = matrix(0, nrow = nrow(scores), ncol = nrow(loadings))
  }
  return(r)
}

#' @title Explained variance ratio of the scores
#' @description Computes the percentage of explained variability based on the variance of the scores matrix
#' @param scores Matrix of scores.
#' @return The percentage of variability each component is explaining.
explained_variance_ratio <- function(scores)
{
  variability <- diag(stats::var(scores))
  pve <- variability / sum(variability)
  return(pve)
}

#' @title Select the number of components
#' @description Computes the required number of components to achieve the desired pve value
#' @param scores matrix of scores
#' @param pve Percentage of explained variability, a number between 0 and 1. Set to NULL to use all components.
#' @param npc Number of components. If provided (non NULL), it supersedes \code{pve}. Must be an integer between 0 and the number of components in \code{scores}.
select_npc <- function(scores, pve, npc = NULL)
{
  if(!is.null(npc)){
    if(!is.numeric(npc) || length(npc) != 1 || npc %% 1 != 0 || npc < 0 || npc > ncol(scores)){
      stop('npc must be an integer number between 0 and the number of components (', ncol(scores),'). Value provided: ', npc)
    }
    return(npc)
  }
  if(is.null(pve)){
    n.components <- ncol(scores)
  }else if(pve == 0){
    n.components <- 0
  }else if(pve > 0 && pve < 1){
    score.variance <- cumsum(explained_variance_ratio(scores))
    n.components <- min(which(score.variance >= pve))
  }else if(pve == 1){
    n.components <- ncol(scores)
  }else{
    stop('pve must be a number between 0 and 1 (percentage of explained variability). To select a number of components directly, use the npc argument. Value provided: ', pve)
  }
  return(n.components)
}
