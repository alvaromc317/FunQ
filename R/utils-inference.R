
# INFERENCE UTILITIES ----------------------------------------------------------

#' @title Compute variance
#' @description Generic function to compute variance.
#' @param object An object for which to compute variance.
#' @param ... further arguments passed to or from other methods.
#' @export
compute_variance <- function(object, ...) {
  UseMethod("compute_variance")
}

#' @title Compute covariance matrix of quantile regression coefficients using kernel estimation
#' @description Replicates the logic of \code{summary.rq(..., se = 'ker')$cov} from a design matrix, response, coefficient vector, and quantile probability, without requiring a \code{quantreg::rq} model object. This lets the sandwich covariance be computed identically regardless of whether the coefficients were estimated with \code{quantreg::rq} or \code{conquer::conquer}.
#' @param x The design matrix (n x p). Should include an intercept column if the fit included one.
#' @param y The response vector (n).
#' @param coefficients The vector of estimated regression coefficients (p).
#' @param tau The quantile probability (scalar, e.g., 0.5 for the median).
#' @param hs Logical, whether to use Hall-Sheather bandwidth selection (default TRUE).
#' @return A p x p covariance matrix of the coefficients.
get_kernel_cov <- function(x, y, coefficients, tau = 0.5, hs = TRUE)
{
  x <- as.matrix(x)
  y <- as.vector(y)
  p <- ncol(x)
  n <- length(y)

  if (length(coefficients) != p) {
    stop("Length of coefficients does not match number of columns in x.")
  }
  if (nrow(x) != n) {
    stop("Number of rows in x does not match length of y.")
  }

  # Calculate residuals
  uhat <- c(y - x %*% coefficients)

  # Initial bandwidth calculation
  h <- quantreg::bandwidth.rq(tau, n, hs = hs)

  # Adjust bandwidth if it extends beyond [0, 1]
  while ((tau - h < 0) || (tau + h > 1)) {
    h <- h / 2
  }

  # Re-scale bandwidth based on residual distribution
  # This matches the implementation in summary.rq for se='ker'
  scale_factor <- min(sqrt(stats::var(uhat)), (stats::quantile(uhat, .75) - stats::quantile(uhat, .25)) / 1.34)
  h_scaled <- (stats::qnorm(tau + h) - stats::qnorm(tau - h)) * scale_factor

  # Compute kernel weights using Gaussian kernel
  f <- stats::dnorm(uhat / h_scaled) / h_scaled
  f <- pmax(f, 1e-8)

  # Compute sandwich components
  # The original implementation uses QR decomposition of weighted X for stability
  qq <- qr(sqrt(f) * x)
  fxxinv <- backsolve(qq$qr[seq_len(p), seq_len(p), drop = FALSE], diag(p))
  fxxinv <- fxxinv %*% t(fxxinv)
  inv_piv <- order(qq$pivot)
  fxxinv <- fxxinv[inv_piv, inv_piv]

  # Compute final covariance matrix
  # cov = tau(1-tau) * fxxinv * (X'X) * fxxinv
  cov_matrix <- tau * (1 - tau) * fxxinv %*% crossprod(x) %*% fxxinv
  cov_matrix <- (cov_matrix + t(cov_matrix)) / 2

  return(cov_matrix)
}

#' @title Compute sandwich estimator variance associated to fosqr loadings
#' @description Given a fosqr-fqpca splines coefficients model, estimates variances associated to fosqr loadings.
#' @param cov.model Covariance matrix of the spline coefficients, as computed by \code{get_kernel_cov}.
#' @param n.regressors number of covariates
#' @param spline.basis spline basis matrix
#' @return A matrix of variances associated to the fosqr loadings
compute_var_sandwich <- function(
    cov.model,
    n.regressors,
    spline.basis)
{
  # Dimensions from the basis
  n.time  <- nrow(spline.basis)
  n.basis <- ncol(spline.basis)

  # Basis used by the intercept block (scalar intercept + non-constant spline cols)
  B_int <- cbind(1, spline.basis[, -1, drop = FALSE])   # n.time x n.basis

  # Block starts for the covariance (each block has size n.basis)
  block_starts <- seq.int(1L, by = n.basis, length.out = (n.regressors + 1L))

  # Allocate output: columns = [intercept, slopes...]
  var.matrix <- matrix(NA_real_, nrow = n.time, ncol = (n.regressors + 1L))

  # Intercept block variance
  idx <- block_starts[1]:(block_starts[1] + n.basis - 1L)
  cov.block <- cov.model[idx, idx, drop = FALSE]
  var.matrix[, 1] <- rowSums((B_int %*% cov.block) * B_int)

  # Slope blocks variance
  for (k in seq_len(n.regressors)) {
    idxk <- block_starts[k + 1L]:(block_starts[k + 1L] + n.basis - 1L)
    cov.block.k <- cov.model[idxk, idxk, drop = FALSE]
    var.matrix[, k + 1L] <- rowSums((spline.basis %*% cov.block.k) * spline.basis)
  }
  var.matrix <- pmax(var.matrix, 0)
  return(var.matrix)
}

#' @title Compute  correction variance associated to fosqr loadings
#' @description Given a fosqr-fqpca splines coefficients model, estimates variances associated to fosqr loadings.
#' @param n.bootstrap number of bootstraps
#' @param fqpca.scores matrix of fqpca scores
#' @param regressors matrix of covariates
#' @param fqpca.loadings matrix of fqpca loadings
#' @return A matrix of variances associated to the fosqr loadings
compute_var_correction <- function(
    n.bootstrap,
    fqpca.scores,
    regressors,
    fqpca.loadings)
{
  regressors <- cbind(1, regressors)
  n.obs <- nrow(fqpca.scores)
  npc <- ncol(fqpca.scores)
  n.regressors <- ncol(regressors)
  # Compute bootstrap coefficients
  bootstrap_coefs_t <- sapply(seq_len(n.bootstrap), function(b.idx){
    idx <- sample(seq_len(n.obs), size = nrow(fqpca.scores), replace = TRUE)
    scores.b = fqpca.scores[idx, , drop=FALSE]
    regressors.b = regressors[idx, , drop=FALSE]
    coefs <- stats::lm.fit(x = regressors.b, y = scores.b)$coefficients
    as.vector(coefs)
  })
  bootstrap.coefs.matrix <- t(bootstrap_coefs_t)
  coef_names <- paste0(
    "B",
    0:(n.regressors-1),
    "_PC",
    rep(seq_len(npc), each = n.regressors))
  colnames(bootstrap.coefs.matrix) <- coef_names
  # Compute variance-covariance matrix
  var_cov_matrix <- stats::cov(bootstrap.coefs.matrix)
  # Compute variance for each regressor
  variance.list <- lapply(seq_len(n.regressors), function(j) {
    regressor_j_coef_names <- paste0(
      "B",
      j-1,
      "_PC",
      seq_len(npc))
    cov_j <- var_cov_matrix[regressor_j_coef_names, regressor_j_coef_names, drop = FALSE]
    rowSums((fqpca.loadings %*% cov_j) * fqpca.loadings)
  })
  var.matrix <- matrix(unlist(variance.list), ncol = n.regressors)
  return(var.matrix)
}
