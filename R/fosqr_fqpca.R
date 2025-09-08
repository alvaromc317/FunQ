# SCORES ----------------------------------------------------------------------

#' @title Compute scores
#' @description Inner function to compute the scores of the fosqr-fqpca methodology.
#' @param Y The \eqn{(N \times T)} matrix of observed time instants.
#' @param Y.mask Mask matrix of the same dimensions as Y indicating which observations in Y are known.
#' @param Y.fosqr.hat fosqr fitted values
#' @param fqpca.loadings Matrix of loading coefficients.
#' @param quantile.value The quantile considered.
#' @return The matrix of scores.
fosqrfqpca_compute_scores <- function(
    Y,
    Y.mask,
    Y.fosqr.hat,
    fqpca.loadings,
    quantile.value)
{
  n.obs <- base::nrow(Y)
  npc <- base::ncol(fqpca.loadings)
  scores <-  base::matrix(0, n.obs, npc)
  for(i in base::seq(n.obs))
  {
    loadings.i <- fqpca.loadings[Y.mask[i, ], , drop=FALSE]
    Yi <- Y[i, Y.mask[i, ]] - Y.fosqr.hat[i, Y.mask[i, ]]
    if(length(Yi) < ncol(loadings.i)){stop('Observation ', i, ' must have a larger number of timepoints than npcs.')}
    scores[i,] <-  quantreg::rq.fit.br(y=Yi, x=loadings.i, tau=quantile.value)$coefficients
  }
  return(scores)
}

# SPLINES ---------------------------------------------------------------------

#' @title Compute spline coefficients
#' @description Inner function to compute the spline coefficients of the fosqr - fqpca methodology.
#' @param Y.vector vectorised version of Y, the \eqn{(N \times T)} matrix of observed time instants.
#' @param Y.mask Mask matrix of the same dimensions as Y indicating which observations in Y are known.
#' @param scores Initial matrix of regressors and / or scores
#' @param spline.basis The spline basis matrix.
#' @param quantile.value The quantile considered.
#' @param method Method used in the resolution of the quantile regression model. It currently accepts the methods \code{c('conquer', 'quantreg')}.
#' @return The matrix of spline coefficients.
fosqrfqpca_compute_spline_coefficients <- function(
    Y.vector,
    Y.mask,
    scores,
    spline.basis,
    quantile.value,
    method)
{
  n.obs <- base::nrow(scores)
  npc <- base::ncol(scores)
  n.basis <- base::ncol(spline.basis)
  intercept.spline.basis <-  spline.basis[, -1]
  tensor.matrix <- base::matrix(0, nrow=base::length(Y.vector), ncol=(npc+1)*n.basis-1)
  row.idx <- 1
  for(i in base::seq(n.obs))
  {
    scores.i <- scores[i, , drop=FALSE]
    tmp.splines <- spline.basis[Y.mask[i, ], , drop = FALSE]
    tmp.intercept.splines <- intercept.spline.basis[Y.mask[i, ], , drop = FALSE]
    n.time.i <- base::nrow(tmp.splines)
    tensor.matrix[row.idx:(row.idx + n.time.i - 1), ] <- base::cbind(tmp.intercept.splines, base::kronecker(scores.i, tmp.splines))
    row.idx <- row.idx+n.time.i
  }
  if(method == 'conquer'){
    B.vector <- conquer::conquer(Y=Y.vector, X=tensor.matrix, tau=quantile.value)$coeff
  }else if(method == 'quantreg'){
    B.vector <- quantreg::rq(Y.vector~tensor.matrix, tau=quantile.value, method='fnb')$coefficients
  }
  spline.coefficients <- base::matrix(B.vector, byrow=FALSE, ncol=(npc+1))
  return(spline.coefficients)
}

#' @title Compute spline coefficients
#' @description Inner function to compute the spline coefficients of the fosqr - fqpca methodology using quantreg to estimate the variance for the sandwich estimator.
#' @param Y.vector vectorised version of Y, the \eqn{(N \times T)} matrix of observed time instants.
#' @param Y.mask Mask matrix of the same dimensions as Y indicating which observations in Y are known.
#' @param scores Initial matrix of regressors and / or scores
#' @param spline.basis The spline basis matrix.
#' @param quantile.value The quantile considered.
#' @return The matrix of spline coefficients.
fosqrfqpca_compute_spline_coefficients_variance_estimation <- function(
    Y.vector,
    Y.mask,
    scores,
    spline.basis,
    quantile.value)
{
  n.obs <- base::nrow(scores)
  npc <- base::ncol(scores)
  n.basis <- base::ncol(spline.basis)
  intercept.spline.basis <-  spline.basis[, -1]
  tensor.matrix <- base::matrix(0, nrow=base::length(Y.vector), ncol=(npc+1)*n.basis-1)
  row.idx <- 1
  for(i in base::seq(n.obs))
  {
    scores.i <- scores[i, , drop=FALSE]
    tmp.splines <- spline.basis[Y.mask[i, ], , drop = FALSE]
    tmp.intercept.splines <- intercept.spline.basis[Y.mask[i, ], , drop = FALSE]
    n.time.i <- base::nrow(tmp.splines)
    tensor.matrix[row.idx:(row.idx + n.time.i - 1), ] <- base::cbind(tmp.intercept.splines, base::kronecker(scores.i, tmp.splines))
    row.idx <- row.idx+n.time.i
  }
  model <- quantreg::rq(Y.vector~tensor.matrix, tau=quantile.value, method='fn')
  B.vector <- model$coefficients
  spline.coefficients <- base::matrix(B.vector, byrow=FALSE, ncol=(npc+1))
  results <- list(spline.coefficients=spline.coefficients, model=model)
  return(results)
}

# OTHER FUNCTIONS -------------------------------------------------------------

#' @title Compute loadings (aka principal components)
#' @description Inner function to compute the loadings of the fqpca methodology.
#' @param spline.basis The spline basis matrix.
#' @param spline.coefficients the matrix of spline coefficients.
#' @return The matrix of loadings
fosqrfqpca_compute_loadings_fosqr <- function(
    spline.basis,
    spline.coefficients)
{
  intercept.spline.basis <-  spline.basis[, -1]
  intercept.part <- c(cbind(1, intercept.spline.basis) %*% matrix(spline.coefficients[,1], ncol=1))
  fqpc.part <- spline.basis %*% spline.coefficients[, -1, drop=FALSE]
  return(list(intercept=intercept.part, loadings=fqpc.part))
}

#' @title Rotation of fqpca loadings and scores
#' @description Performs the rotation of loadings and scores in order to ensure the solution is unique.
#' @param fqpca.loadings Matrix of fqpca loadings.
#' @param fosqr.intercept population intercept
#' @param fqpca.scores Matrix of fqpca scores.
#' @return The rotated matrices of loadings and scores and the rotation matrix.
fosqrfqpca_rotate_scores_and_loadings <- function(
    fqpca.loadings,
    fosqr.intercept,
    fqpca.scores)
{
  npc <- base::ncol(fqpca.loadings)
  # Ensure scores are mean-centered and move the effect of the mean scores to the fosqr intercept
  means.scores <- base::matrix(colMeans(fqpca.scores), ncol=1)
  fqpca.scores <- base::scale(fqpca.scores, center=TRUE, scale=FALSE)
  fosqr.intercept <- fosqr.intercept + fqpca.loadings %*% means.scores
  # Make covariance matrix positive definite using eigen-decomposition
  cov.score <- stats::cov(fqpca.scores)
  eig <- eigen(cov.score, symmetric=TRUE)
  eig$values[eig$values < 1e-8] <- 1e-8
  if (npc == 1) {
    eig_vals_diag <- eig$values
  } else {
    eig_vals_diag <- diag(eig$values)
  }
  cov.score.pd <- eig$vectors %*% eig_vals_diag %*% t(eig$vectors)
  chol_cov <- chol(cov.score.pd)
  svd.decomp <- base::svd((chol_cov)  %*% t(fqpca.loadings))
  rotation.matrix <- base::solve(chol_cov) %*% svd.decomp$u %*% base::diag(svd.decomp$d, nrow=length(svd.decomp$d))
  fqpca.scores <- fqpca.scores %*% rotation.matrix
  fqpca.loadings <- svd.decomp$v
  results <- list(fqpca.loadings=fqpca.loadings, fqpca.scores=fqpca.scores, fosqr.intercept=fosqr.intercept, rotation.matrix=rotation.matrix)
  return(results)
}

#' @title Checks input data
#' @description Checks that Y is valid and formats it.
#' @param Y An \eqn{(N \times T)} matrix or a tf object from the tidyfun package
#' @return The prepared Y matrix
check_Y <- function(Y)
{
  # Process Y: must be a matrix or tf object
  if (is.matrix(Y)) {
    Y_out <- base::unname(Y)
  } else if (tf::is_tf(Y)) {
    Y_out <- base::unname(as.matrix(Y))
  } else {
    stop("Y must be a matrix or a 'tf' object.")
  }
  return(Y_out)
}

#' @title Checks input data
#' @description Checks that regressors is valid and formats it.
#' @param regressors An \eqn{(N \times P)} matrix
#' @return The prepared regressors matrix
check_regressors <- function(regressors)
{
  if (!is.matrix(regressors)) {
    stop("regressors must be a matrix.")
  }
  if (!is.numeric(regressors)) {
    stop("regressors must be a numeric matrix.")
  }
  if (any(is.na(regressors))) {
    stop("regressors cannot contain NA values.")
  }

  col_variances <- apply(regressors, 2, stats::var)
  zero_var_count <- sum(col_variances <= 1e-8)

  # Stop if any column is constant
  if (zero_var_count >= 1) {
    stop("There are constant regressor columns.")
  }
  names.regressors <- colnames(regressors)
  regressors <- base::unname(regressors)
  return(list(regressors=regressors, names.regressors=names.regressors))
}

#' @title Check input parameters
#' @description Check input parameters
#' @param npc The number of estimated components.
#' @param quantile.value The quantile considered.
#' @param periodic Boolean indicating if the data is expected to be periodic (start coincides with end) or not.
#' @param splines.df Degrees of freedom for the splines.
#' @param splines.method Method used in the resolution of the splines quantile regression model. It currently accepts the methods \code{c('conquer', 'quantreg')}.
#' @param tol Tolerance on the convergence of the algorithm.
#' @param max.iters Maximum number of iterations.
#' @param verbose Boolean indicating the verbosity.
#' @param seed Seed for the random generator number.
#' @return No return
check_fosqr_fqpca_params <- function(
    npc,
    quantile.value,
    periodic,
    splines.df,
    splines.method,
    tol,
    max.iters,
    verbose,
    seed)
{
  # Check 'npc': integer number, positive
  if (!is.numeric(npc) || length(npc) != 1 || npc %% 1 != 0 || npc <= 0) {
    stop("Invalid input for 'npc': ", npc,
         ". Expected a positive integer number.")
  }

  # Check 'quantile.value': float number in (0, 1)
  if (!is.numeric(quantile.value) || length(quantile.value) != 1 ||
      quantile.value <= 0 || quantile.value >= 1) {
    stop("Invalid input for 'quantile.value': ", quantile.value,
         ". Expected a float number in (0, 1).")
  }

  # Check 'periodic': Boolean
  if (!is.logical(periodic) || length(periodic) != 1) {
    stop("Invalid input for 'periodic': ", periodic,
         ". Expected a Boolean value (TRUE or FALSE).")
  }

  # Check 'splines.df': integer positive number, larger than 'npc'
  if (!is.numeric(splines.df) || length(splines.df) != 1 ||
      splines.df %% 1 != 0 || splines.df < npc) {
    stop("Invalid input for 'splines.df': ", splines.df,
         ". Expected an integer number larger or equal than 'npc' (", npc, ").")
  }

  # Check 'splines.method': must be either 'conquer' or 'quantreg'
  if (!is.character(splines.method) || length(splines.method) != 1 ||
      !splines.method %in% c("conquer", "quantreg")) {
    stop("Invalid input for 'splines.method': '", splines.method,
         "'. Expected 'conquer' or 'quantreg'.")
  }

  # Check 'tol': positive float number
  if (!is.numeric(tol) || length(tol) != 1 || tol <= 0) {
    stop("Invalid input for 'tol': ", tol,
         ". Expected a positive float number.")
  }

  # Check 'max.iters': positive integer number
  if (!is.numeric(max.iters) || length(max.iters) != 1 ||
      max.iters %% 1 != 0 || max.iters <= 0) {
    stop("Invalid input for 'max.iters': ", max.iters,
         ". Expected a positive integer number.")
  }

  # Check 'verbose': Boolean
  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("Invalid input for 'verbose': ", verbose,
         ". Expected a Boolean value (TRUE or FALSE).")
  }

  # Check 'seed': positive integer number or NULL
  if (!(is.null(seed) || (is.numeric(seed) && length(seed) == 1 &&
                          seed %% 1 == 0 && seed > 0))) {
    stop("Invalid input for 'seed': ", seed,
         ". Expected a positive integer number or NULL.")
  }
}

# MAIN ------------------------------------------------------------------------

#' @title FOSQR - FQPCA main function
#' @description Solves the Function on Scalar Quantile Regression (FOSQR) model
#' accounting for residual correlation based on Functional Principal Component
#' Analysis (FQPCA)
#' @param Y An \eqn{(N \times T)} matrix or a tf object from the tidyfun package.
#' @param regressors An \eqn{(N \times P)} matrix of covariates.
#' @param npc The number of estimated components on the FQPCA part (used for modeling residual correlation)
#' @param quantile.value The quantile considered.
#' @param periodic Boolean indicating if the data is expected to be periodic (start coincides with end) or not.
#' @param splines.df Degrees of freedom for the splines.
#' @param splines.method Method used in the resolution of the splines quantile regression model. It currently accepts the methods \code{c('conquer', 'quantreg')}.
#' @param tol Tolerance on the convergence of the algorithm.
#' @param max.iters Maximum number of iterations.
#' @param verbose Boolean indicating the verbosity.
#' @param seed Seed for the random generator number.
#' @return fosqr_fqpca_object
#' @export
#' @examples
#'
#' n.obs = 150
#' n.time = 144
#' time.grid <- seq(0, 2 * pi, length.out = n.time)
#'
#' # Generate regressor and principal component functions
#' b1 = sin(time.grid)
#' pc1 = sin(2 * time.grid)
#'
#' # Generate covariates and scores
#' regressors =  5 * matrix(rnorm(n.obs), ncol=1)
#' scores = matrix(10 * rnorm(n.obs), ncol = 1)
#' epsilon.matrix <- 0.1 * matrix(rnorm(n.obs * n.time), nrow = n.obs)
#'
#' # Build dataset
#' Y = 100 + regressors[, 1] %o% b1 + scores[, 1] %o% pc1 + epsilon.matrix
#'
#' # Add missing observations
#' Y[sample(n.obs*n.time, as.integer(0.2*n.obs*n.time))] <- NA
#'
#' results <- fosqr_fqpca(Y = Y, regressors = regressors, npc = 1, quantile.value = 0.5, seed=1)
fosqr_fqpca <- function(
    Y = NULL,
    regressors = NULL,
    npc = 2,
    quantile.value = 0.5,
    periodic = TRUE,
    splines.df = 10,
    splines.method = 'conquer',
    tol = 1e-3,
    max.iters = 20,
    verbose = FALSE,
    seed = NULL)
{
  global.start.time <- base::Sys.time()

  # Get the input parameters
  formal_names <- base::setdiff(names(formals(sys.function())), "...")
  inputs <- base::mget(formal_names, envir = environment())
  inputs$Y = NULL
  inputs$regressors = NULL

  # Check the input parameters except Y and colname
  check_fosqr_fqpca_params(npc=npc, quantile.value=quantile.value, periodic=periodic,
                          splines.df=splines.df, splines.method=splines.method,
                          tol=tol, max.iters=max.iters,verbose=verbose, seed=seed)

  # Check Y and colname and return an unnamed matrix
  Y <- check_Y(Y)
  regressors_checked <- check_regressors(regressors)
  regressors <- regressors_checked$regressors

  if(nrow(Y) != nrow(regressors)){stop('nrow(Y) and nrow(regressors) must be the same.')}

  # If seed is provided, set seed for computations
  if(!base::is.null(seed)){base::set.seed(seed)}

  n.obs <- base::dim(Y)[1]
  n.time <- base::dim(Y)[2]
  n.regressors <- base::ncol(regressors)
  Y.axis <- base::seq(0, 1, length.out = n.time)

  Y.mask <- !base::is.na(Y)
  Y.list <- lapply(base::seq(n.obs), function(j) Y[j, Y.mask[j, ]])
  Y.vector <- unlist(Y.list, use.names = FALSE)

  function.warnings <- list(splines = FALSE, scores = FALSE, diverging.loop = FALSE, rotation = FALSE)

  spline.basis <- pbs::pbs(Y.axis, degree = 3, df = splines.df, intercept = TRUE, periodic=periodic, Boundary.knots = c(min(Y.axis), max(Y.axis)))
  n.basis <- base::ncol(spline.basis)

  i <- 1
  loop.start.time <- base::Sys.time()

  # fosqr LEVEL
  fosqr.spline.coef <- try(fosqrfqpca_compute_spline_coefficients(Y.vector = Y.vector, Y.mask = Y.mask, scores = regressors, spline.basis = spline.basis, quantile.value = quantile.value, method=splines.method), silent = FALSE)
  if(!is.matrix(fosqr.spline.coef)){stop('Iteration 1. Failed computation of spline coefficients.')}
  fosqr.loadings.list <- fosqrfqpca_compute_loadings_fosqr(spline.basis = spline.basis, spline.coefficients = fosqr.spline.coef)
  fosqr.intercept <- fosqr.loadings.list$intercept
  fosqr.loadings <- fosqr.loadings.list$loadings
  Y.fosqr.hat <- sweep(regressors %*% t(fosqr.loadings), MARGIN = 2, STATS = fosqr.intercept, FUN = "+")

  # FQPCA level
  fqpca.spline.coef <- base::matrix(stats::rnorm(base::dim(spline.basis)[2] * npc, mean = 0, sd = 1), nrow = base::dim(spline.basis)[2], ncol = npc)
  fqpca.loadings <- spline.basis %*% fqpca.spline.coef
  fqpca.scores <- try(fosqrfqpca_compute_scores(Y = Y, Y.mask = Y.mask, Y.fosqr.hat = Y.fosqr.hat, fqpca.loadings = fqpca.loadings, quantile.value = quantile.value), silent = FALSE)
  if(!is.matrix(fqpca.scores)){stop('Iteration 1. Failed computation of scores')}

  # ROTATION PROCESS
  rotation.result <- try(fosqrfqpca_rotate_scores_and_loadings(fqpca.loadings = fqpca.loadings, fosqr.intercept = fosqr.intercept, fqpca.scores = fqpca.scores), silent = FALSE)
  if(inherits(rotation.result, "try-error"))
  {
    warning('Iteration 1. Failed rotation process. Skipping rotation on this iteration.')
    function.warnings$rotation <- TRUE
    rotation.result <- list(fqpca.loadings = fqpca.loadings, fqpca.scores = fqpca.scores, fosqr.intercept = fosqr.intercept, rotation.matrix = diag(npc))
  }
  fqpca.loadings <- rotation.result$fqpca.loadings
  fqpca.scores <- rotation.result$fqpca.scores
  fosqr.intercept <- rotation.result$fosqr.intercept
  rotation.matrix = rotation.result$rotation.matrix

  # Compute objective value function
  objective.function <- compute_objective_value(quantile.value=quantile.value, Y=Y, scores=cbind(regressors, fqpca.scores), intercept=fosqr.intercept, loadings=cbind(fosqr.loadings, fqpca.loadings))
  objective.function.array <- objective.function

  loop.execution.time <- difftime(base::Sys.time(), loop.start.time, units = 'secs')

  if(verbose)
  {
    message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '. Iteration ', 1, ' completed in ', base::round(loop.execution.time, 3), ' seconds',
            '\n', 'Objective function   i   value: ', round(objective.function.array[1], 4),
            '\n', '____________________________________________________________')
  }

  convergence <- FALSE
  for(i in 2:max.iters)
  {
    loop.start.time <- base::Sys.time()

    # Obtain splines coefficients
    spline.coefficients <- try(fosqrfqpca_compute_spline_coefficients(Y.vector = Y.vector, Y.mask = Y.mask, scores = cbind(regressors, fqpca.scores),  spline.basis = spline.basis, quantile.value = quantile.value, method=splines.method), silent = FALSE)
    if(inherits(spline.coefficients, "try-error"))
    {
      warning('Iteration ', i, '. Failed computation of spline coefficients. Providing results from previous iteration.')
      function.warnings$splines <- TRUE
      break
    }
    fqpca.spline.coef <- spline.coefficients[, (n.regressors+2):ncol(spline.coefficients), drop=FALSE]
    fosqr.spline.coef <- spline.coefficients[, 1:(n.regressors+1), drop=FALSE]

    # Compute loadings
    fqpca.loadings <- spline.basis %*% fqpca.spline.coef
    fosqr.loadings.list <- fosqrfqpca_compute_loadings_fosqr(spline.basis = spline.basis, spline.coefficients = fosqr.spline.coef)
    fosqr.intercept <- fosqr.loadings.list$intercept
    fosqr.loadings <- fosqr.loadings.list$loadings
    Y.fosqr.hat <- sweep(regressors %*% t(fosqr.loadings), MARGIN = 2, STATS = fosqr.intercept, FUN = "+")

    # Compute scores
    fqpca.scores <- try(fosqrfqpca_compute_scores(Y = Y, Y.mask = Y.mask, Y.fosqr.hat = Y.fosqr.hat, fqpca.loadings = fqpca.loadings, quantile.value = quantile.value), silent = FALSE)
    if(!is.matrix(fqpca.scores))
    {
      warning('Iteration ', i, '. Failed computation of scores. Providing results from previous iteration.')
      function.warnings$scores <- TRUE
      break
    }

    # Rotate loadings and scores
    rotation.result <- try(fosqrfqpca_rotate_scores_and_loadings(fqpca.loadings = fqpca.loadings, fosqr.intercept = fosqr.intercept, fqpca.scores = fqpca.scores), silent = FALSE)
    if(inherits(rotation.result, "try-error"))
    {
      warning('Iteration ', i, '. Failed rotation process. Skipping rotation on this iteration.')
      function.warnings$rotation <- TRUE
      rotation.result <- list(fqpca.loadings = fqpca.loadings, fqpca.scores = fqpca.scores, fosqr.intercept = fosqr.intercept, rotation.matrix = diag(npc))
    }
    fqpca.loadings <- rotation.result$fqpca.loadings
    fqpca.scores <- rotation.result$fqpca.scores
    fosqr.intercept <- rotation.result$fosqr.intercept
    rotation.matrix = rotation.result$rotation.matrix

    # Compute objective value function
    objective.function <- compute_objective_value(quantile.value=quantile.value, Y=Y, scores=cbind(regressors, fqpca.scores), intercept=fosqr.intercept, loadings=cbind(fosqr.loadings, fqpca.loadings))
    objective.function.array <- c(objective.function.array, objective.function)

    convergence.criteria <- base::abs(objective.function.array[i] - objective.function.array[i-1])
    if(convergence.criteria < tol){convergence = TRUE}

    # Measure computation time
    loop.execution.time <- difftime(base::Sys.time(), loop.start.time, units = 'secs')

    if(verbose)
    {
      message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '. Iteration ', i, ' completed in ', base::round(loop.execution.time, 3), ' seconds',
              '\n', 'Objective function (i-1) value: ', base::round(objective.function.array[i-1], 4),
              '\n', 'Objective function   i   value: ', round(objective.function.array[i], 4),
              '\n', 'Convergence criteria value    : ', base::round(convergence.criteria, 4),
              '\n', '____________________________________________________________')
    }

    if(convergence){break}

    # Avoid computational issues and do early stops if reaching diverging loops.
    if(i>=3 && (objective.function.array[i] > 10 * objective.function.array[2]))
    {
      function.warnings$diverging.loop <- TRUE
      message('Iteration ', i, '. Breaking diverging loop')
      break
    }
  }

  if(verbose & convergence){message("\u2705 Algorithm converged successfully")}
  if(i == max.iters & !convergence){warning('\u274C Algorithm reached maximum number of iterations without convergence. Consider increasing the value of max.iters')}

  # EXPLAINED VARIABILITY
  pve <- compute_explained_variability(fqpca.scores)

  # FINAL RESULTS
  global.execution.time <- difftime(Sys.time(), global.start.time, units = 'secs')

  results <- structure(list(
    fosqr.intercept = fosqr.intercept,
    fosqr.loadings = fosqr.loadings,
    regressors = regressors,
    fqpca.loadings = fqpca.loadings,
    fqpca.scores = fqpca.scores,
    pve = pve,
    objective.function.array = objective.function.array,
    execution.time = global.execution.time,
    function.warnings = function.warnings,
    rotation.matrix = rotation.matrix,
    spline.basis = spline.basis,
    names.regressors = regressors_checked$names.regressors,
    inputs = inputs),
    class = "fosqr_fqpca_object")
  return(results)
}

# PREDICTIONS -----------------------------------------------------------------

#' @title Predict scores
#' @description S3 method for class 'fosqr_fqpca_object' Given a new matrix Y, predicts the value of the scores associated to the given matrix.
#' @param object An object output of the fqpca function.
#' @param newdata The N by T matrix of observed time instants to be tested, or the tf object.
#' @param newregressors The N by p matrix of regressors to be tested
#' @param ... further arguments passed to or from other methods.
#' @return The predicted matrix of scores.
#' @export
#' @examples
#'
#' n.obs = 150
#' n.time = 144
#' time.grid <- seq(0, 2 * pi, length.out = n.time)
#'
#' # Generate regressor and principal component functions
#' b1 = sin(time.grid)
#' pc1 = sin(2 * time.grid)
#'
#' # Generate covariates and scores
#' regressors =  5 * matrix(rnorm(n.obs), ncol=1)
#' scores = matrix(10 * rnorm(n.obs), ncol = 1)
#' epsilon.matrix <- 0.1 * matrix(rnorm(n.obs * n.time), nrow = n.obs)
#'
#' # Build dataset
#' Y = 100 + regressors[, 1] %o% b1 + scores[, 1] %o% pc1 + epsilon.matrix
#'
#' # Add missing observations
#' Y[sample(n.obs*n.time, as.integer(0.2*n.obs*n.time))] <- NA
#' results <- fosqr_fqpca(
#'     Y = Y[1:100, ],
#'     regressors = regressors[1:100, , drop=FALSE],
#'     npc = 1, quantile.value = 0.5, seed=1)
#'
#' predictions <- predict(object = results, newdata = Y[101:150,],
#'     newregressors = regressors[101:150, , drop=FALSE])
#'
predict.fosqr_fqpca_object <- function(object, newdata, newregressors, ...)
{
  if (!base::inherits(object, "fosqr_fqpca_object")){stop('The object must be of class fosqr_fqpca_object')}

  Y <- check_Y(newdata)
  regressors_checked <- check_regressors(newregressors)
  regressors <- regressors_checked$regressors
  Y.fosqr.hat <- sweep(regressors %*% t(object$fosqr.loadings), MARGIN = 2, STATS = object$fosqr.intercept, FUN = "+")

  Y.mask <- !base::is.na(Y)
  n.obs <- base::nrow(Y)
  n.time <- base::ncol(Y)
  if(sum(!Y.mask) == n.obs * n.time){stop('newdata contains no information. Check that there are values different from NA')}

  # Estimation of scores
  scores <- try(fosqrfqpca_compute_scores(Y = Y, Y.mask = Y.mask, Y.fosqr.hat = Y.fosqr.hat, fqpca.loadings = object$fqpca.loadings, quantile.value = object$inputs$quantile.value), silent = FALSE)
  if(!is.matrix(scores))
  {
    stop('Computation of scores failed.')
  }
  return(scores)
}

#' @title Fit Yhat
#' @description S3 method for class 'fosqr_fqpca_object'. Given a fosqr_fqpca_object model, estimates Yhat for different pve values.
#' @param object An object output of the fosqr_fqpca function.
#' @param pve If smaller than 1, taken as percentage of explained variability used in Yhat estimation. If greater than 1, taken as number of components  used in Yhat estimation.
#' @param ... further arguments passed to or from other methods.
#' @return The matrix of fitted values
#' @export
#' @examples
#'
#' n.obs = 150
#' n.time = 144
#' time.grid <- seq(0, 2 * pi, length.out = n.time)
#'
#' # Generate regressor and principal component functions
#' b1 = sin(time.grid)
#' pc1 = sin(2 * time.grid)
#'
#' # Generate covariates and scores
#' regressors =  5 * matrix(rnorm(n.obs), ncol=1)
#' scores = matrix(10 * rnorm(n.obs), ncol = 1)
#' epsilon.matrix <- 0.1 * matrix(rnorm(n.obs * n.time), nrow = n.obs)
#'
#' # Build dataset
#' Y = 100 + regressors[, 1] %o% b1 + scores[, 1] %o% pc1 + epsilon.matrix
#'
#' # Add missing observations
#' Y[sample(n.obs*n.time, as.integer(0.2*n.obs*n.time))] <- NA
#'
#' results <- fosqr_fqpca(Y = Y, regressors=regressors, npc = 1, quantile.value = 0.5, seed=1)
#' Yhat <- fitted(object = results, pve=0.95)
fitted.fosqr_fqpca_object <- function(object, pve=0.95, ...)
{
  if (!base::inherits(object, "fosqr_fqpca_object")){stop('The object must be of class fosqr_fqpca_object')}
  n.components <- obtain_npc(scores=object$fqpca.scores, pve=pve)

  # Build fosqr fitted value
  fosqr.fitted <- sweep(object$regressors %*% t(object$fosqr.loadings), MARGIN = 2, STATS = object$fosqr.intercept, FUN = "+")

  # Build FQPCA fitted value
  if (n.components > 0){
    fqpca.fitted <- object$fqpca.scores[, seq_len(n.components), drop = FALSE] %*% t(object$fqpca.loadings[, seq_len(n.components), drop = FALSE])
  }else{
    fqpca.fitted <- matrix(0, nrow = nrow(object$fqpca.scores), ncol = nrow(object$fqpca.loadings))
  }
  Y.pred <- fosqr.fitted + fqpca.fitted
  return(Y.pred)
}

# BASIC PLOT ------------------------------------------------------------------

#' @title Plot fosqr and fqpca loading functions with ggplot2
#' @description S3 method for class 'fosqr_fqpca_object'. Given a fosqr_fqpca_object model, plots the principal component and fosqr loadings.
#' @param x An object output of the fosqr_fqpca function.
#' @param pve If smaller than 1, taken as percentage of explained variability used in Yhat estimation. If greater than 1, taken as number of components  used in Yhat estimation.
#' @param ... further arguments passed to or from other methods.
#' @return A ggplot object plotting the intercept and FQPC curves.
#' @export
#' @examples
#'
#' n.obs = 150
#' n.time = 144
#' time.grid <- seq(0, 2 * pi, length.out = n.time)
#'
#' # Generate regressor and principal component functions
#' b1 = sin(time.grid)
#' pc1 = sin(2 * time.grid)
#'
#' # Generate covariates and scores
#' regressors = 5 * matrix(rnorm(n.obs), ncol=1)
#' scores = matrix(10 * rnorm(n.obs), ncol = 1)
#' epsilon.matrix <- 0.1 * matrix(rnorm(n.obs * n.time), nrow = n.obs)
#'
#' # Build dataset
#' Y = 100 + regressors[, 1] %o% b1 + scores[, 1] %o% pc1 + epsilon.matrix
#'
#' # Add missing observations
#' Y[sample(n.obs*n.time, as.integer(0.2*n.obs*n.time))] <- NA
#'
#' results <- fosqr_fqpca(Y = Y, regressors=regressors, npc = 1, quantile.value = 0.5, seed=1)
#' plot(x = results, pve=0.95)
plot.fosqr_fqpca_object <- function(x, pve = 0.95, ...)
{
  if (!base::inherits(x, "fosqr_fqpca_object")){stop('The object x must be of class fosqr_fqpca_object')}

  # Ensure zero components when pve == 0
  n.components <- obtain_npc(scores = x$fqpca.scores, pve = pve)

  intercept <- x$fosqr.intercept

  # Intercept data frame
  intercept_df <- data.frame(
    Time = seq_along(intercept),
    Loading = intercept,
    Component = "Intercept",
    stringsAsFactors = FALSE
  )

  # Between components
  fosqr.loadings <- x$fosqr.loadings
  names.regressors = x$names.regressors
  if(is.null(names.regressors)){names.regressors = paste0("B ", seq_len(ncol(fosqr.loadings)))}
  fosqr_df <- data.frame(
    Time = rep(seq_len(nrow(fosqr.loadings)), times = ncol(fosqr.loadings)),
    Loading = as.vector(fosqr.loadings),
    Component = rep(names.regressors, each = nrow(fosqr.loadings)),
    stringsAsFactors = FALSE
  )

  # FQPCA components
  fqpca.loadings <- x$fqpca.loadings
  if (is.null(fqpca.loadings) || n.components <= 0) {
    fqpca_df <- data.frame(Time = integer(0), Loading = numeric(0), Component = character(0), stringsAsFactors = FALSE)
    n.components <- 0
  } else {
    fqpca.loadings.filter <- fqpca.loadings[, seq_len(n.components), drop = FALSE]
    fqpca_df <- data.frame(
      Time = rep(seq_len(nrow(fqpca.loadings.filter)), times = ncol(fqpca.loadings.filter)),
      Loading = as.vector(fqpca.loadings.filter),
      Component = rep(paste0("FQPC ", seq_len(ncol(fqpca.loadings.filter))), each = nrow(fqpca.loadings.filter)),
      stringsAsFactors = FALSE
    )
  }

  # Combine data
  plot_data <- rbind(intercept_df, fosqr_df, fqpca_df)

  # Safe factor levels (avoid 1:0 issue by using seq_len)
  fosqr_levels <- names.regressors
  fqpca_levels  <- if (n.components  > 0) paste0("FQPC ",  seq_len(n.components))  else character(0)
  plot_data$Component <- factor(plot_data$Component, levels = c("Intercept", fosqr_levels, fqpca_levels))

  # Plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$Time, y = .data$Loading)) +
    ggplot2::geom_line(color = "steelblue") +
    ggplot2::facet_wrap(~ Component, scales = "free_y", ncol = 4) +
    ggplot2::labs(title = "Loading Functions", x = "Time", y = "Loading") +
    ggplot2::theme_bw()
  return(p)
}

# INFERENCE SECTION -----------------------------------------------------------

#' @title Compute variance associated to fosqr loadings
#' @description S3 method for class 'fosqr_fqpca_object'. Given a fosqr_fqpca_object model, estimates variances associated to fosqr loadings.
#' @param object An object output of the fosqr_fqpca function.
#' @param Y An \eqn{(N \times T)} matrix or a tf object from the tidyfun package.
#' @param pve If smaller than 1, taken as percentage of explained variability used in Yhat estimation. If greater than 1, taken as number of components  used in Yhat estimation.
#' @param n.bootstrap An integer number used in the bootstrapping part of the variance estimation.
#' @return A list containing the two parts of the variance computation and the final variances.
#' @export
#' @examples
#' n.obs = 150
#' n.time = 144
#' time.grid <- seq(0, 2 * pi, length.out = n.time)
#'
#' # Generate regressor and principal component functions
#' b1 = sin(time.grid)
#' pc1 = sin(2 * time.grid)
#'
#' # Generate covariates and scores
#' regressors = 5 * matrix(rnorm(n.obs), ncol=1)
#' scores = matrix(10 * rnorm(n.obs), ncol = 1)
#' epsilon.matrix <- 0.1 * matrix(rnorm(n.obs * n.time), nrow = n.obs)
#'
#' # Build dataset
#' Y = 100 + regressors[, 1] %o% b1 + scores[, 1] %o% pc1 + epsilon.matrix
#'
#' # Add missing observations
#' Y[sample(n.obs*n.time, as.integer(0.2*n.obs*n.time))] <- NA
#'
#' results <- fosqr_fqpca(Y = Y, regressors=regressors, npc = 1, quantile.value = 0.5, seed=1)
#' variances <- compute_variance_fosqr(object=results, Y=Y, pve=0.95)
compute_variance_fosqr <- function(object, Y, pve=0.95, n.bootstrap=1000)
{
  if (!base::inherits(object, "fosqr_fqpca_object")){stop('The object must be of class fosqr_fqpca_object')}
  Y <- check_Y(Y)
  if(nrow(Y) != nrow(object$regressors)){stop('nrow(Y) does not match nrow(regressors).')}
  if(ncol(Y) != nrow(object$fosqr.loadings)){stop('ncol(Y) does not match nrow(fosqr.loadings)')}

  n.obs <- base::dim(Y)[1]
  n.regressors <- base::ncol(object$regressors)
  Y.mask <- !base::is.na(Y)
  Y.list <- lapply(base::seq(n.obs), function(j) Y[j, Y.mask[j, ]])
  Y.vector <- unlist(Y.list, use.names = FALSE)
  n.components <- obtain_npc(scores=object$fqpca.scores, pve=pve)
  spline.coefficients.results <- try(fosqrfqpca_compute_spline_coefficients_variance_estimation(
    Y.vector = Y.vector,
    Y.mask = Y.mask,
    scores = cbind(object$regressors, object$fqpca.scores[, seq_len(n.components), drop = FALSE]),
    spline.basis = object$spline.basis,
    quantile.value = object$inputs$quantile.value), silent = FALSE)

  spline.coefficients <- spline.coefficients.results$spline.coefficients
  model <- spline.coefficients.results$model

  fqpca.spline.coef <- spline.coefficients[, (n.regressors+2):ncol(spline.coefficients), drop=FALSE]
  fqpca.loadings <- (object$spline.basis %*% fqpca.spline.coef)[, seq_len(n.components), drop = FALSE]


  var.analytical <- compute_var_sandwich(
    n.regressors=n.regressors,
    model=model,
    spline.basis=object$spline.basis,
    se.method='ker')

  var.correction <- compute_var_correction(
    n.bootstrap = n.bootstrap,
    fqpca.scores = object$fqpca.scores[, seq_len(n.components), drop = FALSE],
    regressors = object$regressors,
    fqpca.loadings = fqpca.loadings)
  results <- list(variance = var.analytical+var.correction, var.analytical=var.analytical, var.correction=var.correction)
  return(results)
}

#' @title Compute sandwich estimator variance associated to fosqr loadings
#' @description Given a fosqr-fqpca splines coefficients model, estimates variances associated to fosqr loadings.
#' @param n.regressors number of covariates
#' @param model fosqr-fqpca splines quantile regression model
#' @param spline.basis spline basis matrix
#' @param se.method standard error method used by quantreg to obtain the sandwich estimator
#' @return A matrix of variances associated to the fosqr loadings
compute_var_sandwich <- function(
    n.regressors,
    model,
    spline.basis,
    se.method='ker')
{
  splines.df <- ncol(spline.basis)
  cov.model <- summary(model, covariance = TRUE, se = se.method)$cov
  var.matrix <- sapply(1:(n.regressors+1), function(idx) {
    indices <- ((idx - 1) * splines.df + 1):(idx * splines.df)
    cov.block <- cov.model[indices, indices, drop = FALSE]
    var.loading <- spline.basis %*% cov.block %*% t(spline.basis)
    diag(var.loading)
  })
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
    idx <- sample(1:n.obs, size = nrow(fqpca.scores), replace = TRUE)
    scores.b = fqpca.scores[idx, , drop=FALSE]
    regressors.b = regressors[idx, , drop=FALSE]
    coefs <- stats::lm.fit(x = regressors.b, y = scores.b)$coefficients
    as.vector(coefs)
  })
  bootstrap.coefs.matrix <- t(bootstrap_coefs_t)
  coef_names <- paste0("B", 0:(n.regressors-1), "_PC", rep(1:npc, each = n.regressors))
  colnames(bootstrap.coefs.matrix) <- coef_names
  # Compute variance-covariance matrix
  var_cov_matrix <- stats::cov(bootstrap.coefs.matrix)
  # Compute variance for each regressor
  variance.list <- lapply(1:n.regressors, function(j) {
    regressor_j_coef_names <- paste0("B", j-1, "_PC", 1:npc)
    cov_j <- var_cov_matrix[regressor_j_coef_names, regressor_j_coef_names, drop = FALSE]
    rowSums((fqpca.loadings %*% cov_j) * fqpca.loadings)
  })
  var.matrix <- matrix(unlist(variance.list), ncol = n.regressors)
  return(var.matrix)
}
