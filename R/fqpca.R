# SCORES ----------------------------------------------------------------------

#' @title Compute scores
#' @description Inner function to compute the scores of the fqpca methodology.
#' @param Y The \eqn{(N \times T)} matrix of observed time instants.
#' @param Y.mask Mask matrix of the same dimensions as Y indicating which observations in Y are known.
#' @param intercept population intercept
#' @param loadings Matrix of loading coefficients.
#' @param quantile.value The quantile considered.
#' @return The matrix of scores.
compute_scores <- function(Y, Y.mask, intercept, loadings, quantile.value)
{
  # Initialize the matrix of scores
  n.obs <- base::nrow(Y)
  npc <- base::ncol(loadings)
  scores <-  base::matrix(0, n.obs, npc)
  for(i in base::seq(n.obs))
  {
    # Obtain the loadings associated to observation i with a specific time grid
    loadings.obs <- loadings[Y.mask[i, ], , drop = FALSE]
    intercept.obs <- intercept[Y.mask[i, ]]
    # Obtain the observation removing missing data and subtract intercept
    Yi <- Y[i, Y.mask[i, ]] - intercept.obs
    scores[i,] <-  quantreg::rq.fit.br(y=Yi, x=loadings.obs, tau=quantile.value)$coefficients
  }
  return(scores)
}

# SPLINES ---------------------------------------------------------------------

#' @title Compute spline coefficients
#' @description Inner function to compute the spline coefficients of the fqpca methodology.
#' @param Y.vector vectorised version of Y, the \eqn{(N \times T)} matrix of observed time instants.
#' @param Y.mask Mask matrix of the same dimensions as Y indicating which observations in Y are known.
#' @param scores Initial matrix of scores.
#' @param spline.basis The spline basis matrix.
#' @param quantile.value The quantile considered.
#' @param method Method used in the resolution of the quantile regression model. It currently accepts the methods \code{c('conquer', 'quantreg')}.
#' @return The matrix of spline coefficients.
compute_spline_coefficients_unpenalized <- function(Y.vector, Y.mask, scores, spline.basis, quantile.value, method)
{
  n.obs <- base::nrow(scores)
  npc <- base::ncol(scores)
  n.basis <- base::ncol(spline.basis)
  intercept.spline.basis <-  spline.basis[, -1]
  tensor.matrix <- base::matrix(0, nrow=base::length(Y.vector), ncol=(npc+1)*n.basis-1)
  row.idx <- 1
  for(i in base::seq(n.obs))
  {
    scores.i <- scores[i, , drop=F]
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
#' @description Inner function to compute the spline coefficients of the fqpca methodology.
#' @param Y.vector vectorised version of Y, the \eqn{(N \times T)} matrix of observed time instants.
#' @param Y.mask Mask matrix of the same dimensions as Y indicating which observations in Y are known.
#' @param scores Initial matrix of scores.
#' @param spline.basis The spline basis matrix.
#' @param quantile.value The quantile considered.
#' @param lambda.ridge  Hyper parameter controlling the penalization on the second derivative of the splines. It has effect only with \code{penalized=TRUE} and \code{method='conquer'}.
#' @return The matrix of spline coefficients splines coefficients.
compute_spline_coefficients_penalized <- function(Y.vector, Y.mask, scores, spline.basis, quantile.value, lambda.ridge)
{
  n.obs <- base::nrow(scores)
  npc <- base::ncol(scores)
  n.basis <- base::ncol(spline.basis)
  intercept.spline.basis <-  spline.basis[, -1]
  tensor.matrix <- base::matrix(0, nrow=base::length(Y.vector), ncol=(npc+1)*n.basis-1)
  row.idx <- 1
  for(i in base::seq(n.obs))
  {
    scores.i <- scores[i, , drop=F]
    tmp.splines <- spline.basis[Y.mask[i, ], , drop = FALSE]
    tmp.intercept.splines <- intercept.spline.basis[Y.mask[i, ], , drop = FALSE]
    n.time.i <- base::nrow(tmp.splines)
    tensor.matrix[row.idx:(row.idx + n.time.i - 1), ] <- base::cbind(tmp.intercept.splines, base::kronecker(scores.i, tmp.splines))
    row.idx <- row.idx+n.time.i
  }
  if(lambda.ridge == 0)
  {
    B.vector <- conquer::conquer(Y=Y.vector, X=tensor.matrix, tau=quantile.value)$coeff
  }else{
    B.vector <- conquer::conquer.reg(Y=Y.vector, X=tensor.matrix, tau=quantile.value, penalty='elastic', lambda=lambda.ridge, para.elastic=0)$coeff
  }
  spline.coefficients <- base::matrix(B.vector, byrow=FALSE, ncol=(npc+1))
  return(spline.coefficients)
}

#' @title Compute spline coefficients
#' @description Inner function to compute the spline coefficients of the fqpca methodology.
#' @param penalized Boolean indicating if the smoothness should be controlled using a second derivative penalty.
#' @param Y.vector vectorised version of Y, the \eqn{(N \times T)} matrix of observed time instants.
#' @param Y.mask Mask matrix of the same dimensions as Y indicating which observations in Y are known.
#' @param scores Initial matrix of scores.
#' @param spline.basis The spline basis matrix.
#' @param quantile.value The quantile considered.
#' @param method Method used in the resolution of the quantile regression model. It currently accepts the methods \code{c('conquer', 'quantreg')}.
#' @param lambda.ridge  Hyper parameter controlling the penalization on the second derivative of the splines. It has effect only with \code{penalized=TRUE} and \code{method='conquer'}.
#' @return The matrix of spline coefficients.
compute_spline_coefficients <- function(penalized, Y.vector, Y.mask, scores, spline.basis, quantile.value, method, lambda.ridge)
{
  if(penalized)
  {
    spline.coefficients <- compute_spline_coefficients_penalized(Y.vector=Y.vector, Y.mask=Y.mask, scores=scores, spline.basis=spline.basis, quantile.value=quantile.value, lambda.ridge=lambda.ridge)
  }else{
    spline.coefficients <- compute_spline_coefficients_unpenalized(Y.vector=Y.vector, Y.mask=Y.mask, scores=scores, spline.basis=spline.basis, quantile.value=quantile.value, method=method)
  }
  return(spline.coefficients)
}

# OTHER FUNCTIONS -------------------------------------------------------------

#' @title Compute objective value function.
#' @description Inner function to compute the objective value of the fqpca methodology at each iteration.
#' @param Y \eqn{(N \times T)} matrix of observed time instants.
#' @param quantile.value The quantile considered.
#' @param scores The matrix of estimated scores.
#' @param intercept population intercept
#' @param loadings The matrix of estimated loadings.
#' @return The objective value function.
compute_objective_value <- function(Y, quantile.value, scores, intercept, loadings)
{
  Y.pred <- sweep(scores %*% t(loadings), MARGIN = 2, STATS = intercept, FUN = "+")
  objective.value <- quantile_error(Y=Y, Y.pred=Y.pred, quantile.value=quantile.value)
  return(objective.value)
}

#' @title Compute loadings (aka principal components)
#' @description Inner function to compute the loadings of the fqpca methodology.
#' @param spline.basis The spline basis matrix.
#' @param spline.coefficients the matrix of spline coefficients.
#' @return The matrix of loadings
compute_loadings <- function(spline.basis, spline.coefficients)
{
  intercept.spline.basis <-  spline.basis[, -1]
  intercept.part <- c(cbind(1, intercept.spline.basis) %*% matrix(spline.coefficients[,1], ncol=1))
  fqpc.part <- spline.basis %*% spline.coefficients[, -1, drop=FALSE]
  return(list(intercept=intercept.part, loadings=fqpc.part))
}

#' @title Rotation of fqpca loadings and scores
#' @description Performs the rotation of loadings and scores in order to ensure the solution is unique.
#' @param scores Matrix of scores.
#' @param intercept population intercept
#' @param loadings Matrix of loadings.
#' @return The rotated matrices of loadings and scores and the rotation matrix.
rotate_scores_and_loadings <- function(scores, intercept, loadings)
{
  npc <- base::ncol(loadings)

  # Ensure scores are mean-centered
  means.scores <- base::matrix(colMeans(scores), ncol=1)
  scores <- base::scale(scores, center=TRUE, scale=FALSE)

  # Move the effect of the mean scores to the intercept
  intercept <- intercept + loadings %*% means.scores

  # Make covariance matrix positive definite using eigen-decomposition
  cov.score <- stats::cov(scores)
  eig <- eigen(cov.score, symmetric=TRUE)
  eig$values[eig$values < 1e-8] <- 1e-8
  if (npc == 1) {
    eig_vals_diag <- eig$values
  } else {
    eig_vals_diag <- diag(eig$values)
  }
  cov.score.pd <- eig$vectors %*% eig_vals_diag %*% t(eig$vectors)

  # Perform rotation
  svd.decomp <- svd(chol(cov.score.pd) %*% t(loadings))
  rotation.matrix <- solve(chol(cov.score.pd)) %*% svd.decomp$u %*%
    diag(svd.decomp$d, nrow=length(svd.decomp$d))
  scores <- scores %*% rotation.matrix
  loadings <- svd.decomp$v

  results <- list(scores=scores, intercept=intercept, loadings=loadings, rotation.matrix=rotation.matrix)
  return(results)
}

#' @title Check input parameters
#' @description Check input parameters
#' @param npc The number of estimated components.
#' @param quantile.value The quantile considered.
#' @param periodic Boolean indicating if the data is expected to be periodic (start coincides with end) or not.
#' @param splines.df Degrees of freedom for the splines.
#' @param splines.method Method used in the resolution of the splines quantile regression model. It currently accepts the methods \code{c('conquer', 'quantreg')}.
#' @param penalized Boolean indicating if the smoothness should be controlled using a second derivative penalty. This functionality is experimental and is much slower than the control of the smoothness using the degrees of freedom.
#' @param lambda.ridge  Hyper parameter controlling the penalization on the second derivative of the splines. It has effect only with \code{penalized=TRUE}.
#' @param tol Tolerance on the convergence of the algorithm.
#' @param max.iters Maximum number of iterations.
#' @param verbose Boolean indicating the verbosity.
#' @param seed Seed for the random generator number.
#' @return No return
check_fqpca_params <- function(npc, quantile.value, periodic, splines.df, splines.method, penalized, lambda.ridge, tol, max.iters, verbose, seed)
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

  # Check 'penalized': Boolean
  if (!is.logical(penalized) || length(penalized) != 1) {
    stop("Invalid input for 'penalized': ", penalized,
         ". Expected a Boolean value (TRUE or FALSE).")
  }

  # Check 'lambda.ridge': positive number or zero
  if (!is.numeric(lambda.ridge) || length(lambda.ridge) != 1 ||
      lambda.ridge < 0) {
    stop("Invalid input for 'lambda.ridge': ", lambda.ridge,
         ". Expected a positive number or zero.")
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

#' @title Check input data
#' @description Check input data
#' @param data An \eqn{(N \times T)} matrix, a tf object from the tidyfun package or a data.frame containing the functional data as a tf column.
#' @param colname The name of the column containing the functional data. Use only if data is a dataframe and colname is a column in the dataframe.
#' @return an unnamed matrix
check_input_data <- function(data, colname)
{
  # Case 1: If data is a matrix, unname and return
  if(is.matrix(data)){return(base::unname(data))}

  # Case 2: If data is a 'tf' object
  if(tf::is_tf(data)){return(base::unname(base::as.matrix(data)))}

  # Case 3: If data is a data.frame
  if(is.data.frame(data))
  {
    # Check if 'colname' is provided and is a string
    if(missing(colname) || !is.character(colname) || length(colname) != 1)
    {
      stop("Data is a data.frame, but 'colname' is not provided as a string corresponding to a column name.")
    }
    # Check if 'colname' exists in the data.frame
    if(!(colname %in% colnames(data)))
    {
      stop("Data is a data.frame, but 'colname' ('", colname, "') does not correspond to any column in the data.")
    }
    # Extract the specified column
    column_data <- data[[colname]]
    # Check if the extracted column is a 'tf' object
    if(!tf::is_tf(column_data))
    {
      stop("The column '", colname, "' in the data.frame is not a 'tf' object.")
    }
    return(unname(as.matrix(column_data)))
  }
  # If data is none of the above types, stop with an error
  stop("Data is not a matrix, 'tf' object, or data.frame. Cannot process the input data.")
}

# MAIN ------------------------------------------------------------------------

#' @title FQPCA (Functional Quantile Principal Component Analysis)
#' @description Solves the functional quantile principal component analysis methodology
#' @param data An \eqn{(N \times T)} matrix, a tf object from the tidyfun package or a data.frame containing the functional data as a tf column.
#' @param colname The name of the column containing the functional data. Use only if data is a dataframe and colname is a column in the dataframe.
#' @param npc The number of estimated components.
#' @param quantile.value The quantile considered.
#' @param periodic Boolean indicating if the data is expected to be periodic (start coincides with end) or not.
#' @param splines.df Degrees of freedom for the splines.
#' @param splines.method Method used in the resolution of the splines quantile regression model. It currently accepts the methods \code{c('conquer', 'quantreg')}.
#' @param penalized Boolean indicating if the smoothness should be controlled using a second derivative penalty. This functionality is experimental.
#' @param lambda.ridge  Hyper parameter controlling the penalization on the second derivative of the splines. It has effect only with \code{penalized=TRUE} and \code{method='conquer'}.
#' @param tol Tolerance on the convergence of the algorithm.
#' @param max.iters Maximum number of iterations.
#' @param verbose Boolean indicating the verbosity.
#' @param seed Seed for the random generator number.
#' @return fqpca_object
#' @export
#' @examples
#'
#' n.obs = 150
#' n.time = 144
#'
#' # Generate scores
#' c1.vals = rnorm(n.obs)
#' c2.vals = rnorm(n.obs)
#'
#' # Generate pc's
#' pc1 = sin(seq(0, 2*pi, length.out = n.time))
#' pc2 = cos(seq(0, 2*pi, length.out = n.time))
#'
#' # Generate data
#' Y <- c1.vals * matrix(pc1, nrow = n.obs, ncol=n.time, byrow = TRUE) +
#' c2.vals * matrix(pc2, nrow = n.obs, ncol=n.time, byrow = TRUE)
#'
#' # Add noise
#' Y <- Y + matrix(rnorm(n.obs * n.time, 0, 0.4), nrow = n.obs)
#'
#' # Add missing observations
#' Y[sample(n.obs*n.time, as.integer(0.2*n.obs*n.time))] <- NA
#'
#' results <- fqpca(data = Y, npc = 2, quantile.value = 0.5, seed=1)
#'
#' intercept <- results$intercept
#' loadings <- results$loadings
#' scores <- results$scores
fqpca <- function(
    data, colname = NULL, npc = 2,  quantile.value = 0.5,  periodic = TRUE,
    splines.df = 10, splines.method = 'conquer', penalized = FALSE,
    lambda.ridge = 0, tol = 1e-3, max.iters = 20, verbose = FALSE, seed = NULL)
{
  global.start.time <- base::Sys.time()

  # Get the input parameters
  formal_names <- base::setdiff(names(formals(sys.function())), "...")
  inputs <- base::mget(formal_names, envir = environment())
  inputs$data = NULL
  inputs$colname = NULL

  # Check the input parameters except Y and colname
  check_fqpca_params(npc=npc, quantile.value=quantile.value, periodic=periodic,
                     splines.df=splines.df, splines.method=splines.method, penalized=penalized,
                     lambda.ridge=lambda.ridge, tol=tol, max.iters=max.iters,
                     verbose=verbose, seed=seed)

  # Check Y and colname and return an unnamed matrix
  Y <- check_input_data(data=data, colname=colname)

  # Check if splines.method is valid
  if(penalized && splines.method != 'conquer')
  {
    stop("Invalid input for 'splines.method and/or penalized':
         If penalized is set to  TRUE then splines.method must be 'conquer'")
  }

  # If seed is provided, set seed for computations
  if(!base::is.null(seed)){base::set.seed(seed)}

  # Step 1: Initialize values
  n.obs <- base::nrow(Y)
  n.time <- base::dim(Y)[2]
  Y.axis <- base::seq(0, 1, length.out = n.time)
  Y.mask <- !base::is.na(Y)
  # vectorize Y [Y11, ..., Y1T, Y21, ..., Y2T,...YNT]
  Y.list <- lapply(base::seq(n.obs), function(j) Y[j, Y.mask[j, ]])
  Y.vector <- unlist(Y.list, use.names = FALSE)
  function.warnings <- list(splines = FALSE, scores = FALSE, rotation = FALSE, diverging.loop = FALSE)

  # Initialize splines knots
  if(periodic)
  {
    knots <- base::seq(0, 1, length.out = 2 + splines.df - 1)  # 2 boundary knots + df - intercept
  } else{
    knots <- base::seq(0, 1, length.out = 2 + splines.df - 3 - 1) # 2 boundary knots + df - degree - intercept
  }
  knots <- knots[2:(base::length(knots)-1)]

  # Initialize spline basis
  spline.basis <- pbs::pbs(Y.axis, degree = 3, knots = knots, intercept = T, periodic=periodic, Boundary.knots = c(min(Y.axis), max(Y.axis)))
  n.basis <- base::ncol(spline.basis)

  # If penalized, compute different basis with an identity second derivative matrix with 0s in first and second positions
  # Based on Garcia de la Garza et al. (2024) (Adaptive Functional Principal Component Analysis)
  if(penalized)
  {
    eigen.decomp <- base::eigen(base::crossprod(base::diff(spline.basis, differences = 2)))
    first.part <- eigen.decomp$vectors[, (n.basis-1):n.basis]
    safe_eigen_values <- pmax(eigen.decomp$values[1:(n.basis-2)], 1e-8) # or a suitable small constant
    second.part <- eigen.decomp$vectors[, 1:(n.basis-2)] %*% base::diag(1 / base::sqrt(safe_eigen_values))
    U <- cbind(first.part, second.part)
    spline.basis <- spline.basis %*% U
  }

  i <- 1
  loop.start.time <- base::Sys.time()

  # Randomly generate splines coefficients
  spline.coefficients <- base::matrix(stats::rnorm(n.basis * (npc+1), mean = 0, sd = 1), nrow=n.basis, ncol=npc+1)

  loadings_list <- spline.basis %*% spline.coefficients
  intercept <- loadings_list[, 1]
  loadings <- loadings_list[, -1, drop=FALSE]

  scores <- try(compute_scores(Y = Y, Y.mask = Y.mask, intercept = intercept, loadings = loadings, quantile.value = quantile.value), silent = FALSE)
  if(!is.matrix(scores)){stop('Iteration 1. Failed computation of scores')}

  # Rotate loadings and scores
  rotation.result <- try(rotate_scores_and_loadings(scores = scores, intercept = intercept, loadings = loadings), silent = FALSE)
  if(inherits(rotation.result, "try-error"))
  {
    warning('Iteration 1. Failed rotation process. Skipping rotation on this iteration.')
    function.warnings$rotation <- TRUE
    rotation.result <- list(scores = scores, intercept = intercept, loadings = loadings, rotation.matrix = diag(npc))
  }
  scores <- rotation.result$scores
  intercept <- rotation.result$intercept
  loadings <- rotation.result$loadings
  rotation.matrix = rotation.result$rotation.matrix

  # Compute objective value function
  objective.function.array <- compute_objective_value(quantile.value = quantile.value, Y = Y, scores = scores, intercept = intercept,  loadings = loadings)

  last.iteration <- list(scores = scores, intercept = intercept, loadings = loadings, spline.coefficients = spline.coefficients, objective.function=objective.function.array[1], rotation.matrix = rotation.matrix, iteration = 1)

  loop.execution.time <- difftime(base::Sys.time(), loop.start.time, units = 'secs')

  if(verbose)
  {
    message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '. Iteration ', 1, ' completed in ', base::round(loop.execution.time, 3), ' seconds',
            '\n', 'Objective function   i   value: ', round(objective.function.array[1], 4),
            '\n', '____________________________________________________________')
  }

  # Iterate until convergence or until max.iters is reached
  convergence <- F
  for(i in 2:max.iters)
  {
    loop.start.time <- base::Sys.time()

    # OBTAIN SPLINE COEFFICIENTS
    spline.coefficients <- try(compute_spline_coefficients(penalized = penalized, Y.vector = Y.vector, Y.mask = Y.mask, scores = scores, spline.basis = spline.basis, quantile.value = quantile.value, method = splines.method, lambda.ridge = lambda.ridge), silent = FALSE)
    if(!is.matrix(spline.coefficients))
    {
      warning('Iteration ', i, '. Failed computation of spline coefficients. Providing results from previous iteration.')
      function.warnings$splines <- TRUE
      break
    }

    loadings_list <- compute_loadings(spline.basis=spline.basis, spline.coefficients=spline.coefficients)
    intercept <- loadings_list$intercept
    loadings <- loadings_list$loadings
    if(is_rank_deficient(loadings)){warning('Loadings matrix is singular.')}

    # OBTAIN SCORES
    scores <- try(compute_scores(Y = Y, Y.mask = Y.mask, intercept = intercept, loadings = loadings, quantile.value = quantile.value), silent = FALSE)
    if(!is.matrix(scores))
    {
      warning('Iteration ', i, '. Failed computation of scores. Providing results from previous iteration.')
      function.warnings$scores <- TRUE
      break
    }

    # Rotate loadings and scores
    rotation.result <- try(rotate_scores_and_loadings(scores = scores, intercept = intercept, loadings = loadings), silent = FALSE)
    if(inherits(rotation.result, "try-error"))
    {
      warning('Iteration ', i, '. Failed rotation process. Skipping rotation on this iteration.')
      function.warnings$rotation <- TRUE
      rotation.result <- list(scores = scores, intercept = intercept, loadings = loadings, rotation.matrix = diag(npc))
    }
    scores <- rotation.result$scores
    intercept <- rotation.result$intercept
    loadings <- rotation.result$loadings
    rotation.matrix = rotation.result$rotation.matrix

    # Compute objective function
    objective.function <- compute_objective_value(quantile.value = quantile.value, Y = Y, scores = scores, intercept = intercept, loadings = loadings)
    objective.function.array <- c(objective.function.array, objective.function)

    last.iteration <- list(scores = scores, intercept = intercept, loadings = loadings, spline.coefficients = spline.coefficients, objective.function=objective.function.array[i], rotation.matrix = rotation.matrix, iteration = i)

    convergence.criteria <- base::abs(objective.function.array[i] - objective.function.array[i-1])
    if(convergence.criteria < tol){convergence = T}

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
    if(i>=3 & (objective.function.array[i] > 10 * objective.function.array[2]))
    {
      function.warnings$diverging.loop <- TRUE
      message('Iteration ', i, '. Breaking diverging loop')
      break
    }
  }

  if(verbose & convergence){message("\u2705 Algorithm converged succesfully")}

  if(i == max.iters & !convergence){warning('\u274C Algorithm reached maximum number of iterations without convergence. Consider increasing the value of max.iters parameter')}

  # EXPLAINED VARIABILITY -----------------------------------------------------

  pve <- compute_explained_variability(last.iteration$scores)

  global.execution.time <- difftime(base::Sys.time(), global.start.time, units = 'secs')

  results <- structure(list(
    scores = last.iteration$scores,
    intercept = last.iteration$intercept,
    loadings = last.iteration$loadings,
    pve = pve,
    objective.function.value = last.iteration$objective.function,
    objective.function.array = objective.function.array,
    function.warnings = function.warnings,
    execution.time = global.execution.time,
    rotation.matrix = last.iteration$rotation.matrix,
    spline.coefficients = last.iteration$spline.coefficients,
    spline.basis = spline.basis,
    n.iters = i,
    inputs = inputs),
    class = "fqpca_object")

  return(results)
}

# PREDICTIONS -----------------------------------------------------------------

#' @title Predict fqpca scores
#' @description S3 method for class 'fqpca_object' Given a new matrix Y, predicts the value of the scores associated to the given matrix.
#' @param object An object output of the fqpca function.
#' @param newdata The N by T matrix of observed time instants to be tested, or the dataframe storing the tf functional vector using the same colname as the one used in the fqpca function.
#' @param ... further arguments passed to or from other methods.
#' @return The normalized matrix of scores.
#' @export
#' @examples
#'
#' n.obs = 150
#' n.time = 144
#'
#' # Generate scores
#' c1.vals = rnorm(n.obs)
#' c2.vals = rnorm(n.obs)
#'
#' # Generate pc's
#' pc1 = sin(seq(0, 2*pi, length.out = n.time))
#' pc2 = cos(seq(0, 2*pi, length.out = n.time))
#'
#' # Generate data
#' Y <- c1.vals * matrix(pc1, nrow = n.obs, ncol=n.time, byrow = TRUE) +
#' c2.vals * matrix(pc2, nrow = n.obs, ncol=n.time, byrow = TRUE)
#'
#' # Add noise
#' Y <- Y + matrix(rnorm(n.obs * n.time, 0, 0.4), nrow = n.obs)
#'
#' # Add missing observations
#' Y[sample(n.obs*n.time, as.integer(0.2*n.obs*n.time))] <- NA
#'
#' results <- fqpca(data = Y, npc = 2, quantile.value = 0.5, seed=1)
#'
#' predictions <- predict(object = results, newdata = Y[101:150,])
#'
predict.fqpca_object <- function(object, newdata, ...)
{
  if (!base::inherits(object, "fqpca_object")){stop('The object must be of class fqpca_object')}

  Y <- check_input_data(data=newdata, colname=object$inputs$colname)

  Y.mask <- !base::is.na(Y)
  n.obs <- base::nrow(Y)
  n.time <- base::ncol(Y)
  if(sum(!Y.mask) == n.obs * n.time){stop('newdata contains no information. Check that there are values different from NA')}

  # Estimation of scores
  scores <- try(compute_scores(Y = Y, Y.mask = Y.mask, intercept = object$intercept, loadings = object$loadings, quantile.value = object$inputs$quantile.value), silent = FALSE)
  if(!is.matrix(scores))
  {
    stop('Computation of scores failed.')
  }
  return(scores)
}

#' @title Fit Yhat
#' @description S3 method for class 'fqpca_object'. Given an fqpca_object model, estimates Yhat for different pve values.
#' @param object An object output of the fqpca function.
#' @param pve If smaller than 1, taken as percentage of explained variability used in Yhat estimation. If greater than 1, taken as number of components  used in Yhat estimation.
#' @param ... further arguments passed to or from other methods.
#' @return The normalized matrix of scores.
#' @export
#' @examples
#'
#' n.obs = 150
#' n.time = 144
#'
#' # Generate scores
#' c1.vals = rnorm(n.obs)
#' c2.vals = rnorm(n.obs)
#'
#' # Generate pc's
#' pc1 = sin(seq(0, 2*pi, length.out = n.time))
#' pc2 = cos(seq(0, 2*pi, length.out = n.time))
#'
#' # Generate data
#' Y <- c1.vals * matrix(pc1, nrow = n.obs, ncol=n.time, byrow = TRUE) +
#' c2.vals * matrix(pc2, nrow = n.obs, ncol=n.time, byrow = TRUE)
#'
#' # Add noise
#' Y <- Y + matrix(rnorm(n.obs * n.time, 0, 0.4), nrow = n.obs)
#'
#' # Add missing observations
#' Y[sample(n.obs*n.time, as.integer(0.2*n.obs*n.time))] <- NA
#'
#' results <- fqpca(data = Y, npc = 2, quantile.value = 0.5, seed=1)
#'
#' Yhat <- fitted(object = results, pve=0.99)
#'
fitted.fqpca_object <- function(object, pve=0.95, ...)
{
  if (!base::inherits(object, "fqpca_object")){stop('The object must be of class fqpca_object')}
  n.components <- obtain_npc(scores=object$scores, pve=pve)
  Y.pred <- sweep(object$scores[, 1:n.components, drop=F] %*% t(object$loadings[, 1:n.components, drop=F]), MARGIN = 2, STATS = object$intercept, FUN = "+")
  return(Y.pred)
}

# BASIC PLOT ------------------------------------------------------------------

#' @title Plot fqpca loading functions with ggplot2
#' @description S3 method for class 'fqpca_object'. Given a fqpca object, plot
#'              the loading functions using ggplot2 facets.
#' @param x An object output of the fqpca function.
#' @param pve Percentage of explained variability to determine the number of components.
#' @param ... Further arguments passed to or from other methods.
#' @return A ggplot object plotting the intercept and FQPC curves.
#' @importFrom rlang .data
#' @export
#' @examples
#'
#' n.obs = 150
#' n.time = 144
#'
#' # Generate scores
#' c1.vals = rnorm(n.obs)
#' c2.vals = rnorm(n.obs)
#'
#' # Generate pc's
#' pc1 = sin(seq(0, 2*pi, length.out = n.time))
#' pc2 = cos(seq(0, 2*pi, length.out = n.time))
#'
#' # Generate data
#' Y <- c1.vals * matrix(pc1, nrow = n.obs, ncol=n.time, byrow = TRUE) +
#' c2.vals * matrix(pc2, nrow = n.obs, ncol=n.time, byrow = TRUE) +
#' matrix(seq(1, 10, length.out=n.time), nrow = n.obs, ncol=n.time, byrow = TRUE)
#'
#' # Add noise
#' Y <- Y + matrix(rnorm(n.obs * n.time, 0, 0.4), nrow = n.obs)
#'
#' # Add missing observations
#' Y[sample(n.obs*n.time, as.integer(0.2*n.obs*n.time))] <- NA
#'
#' results <- fqpca(data = Y, npc = 2, quantile.value = 0.5, periodic = FALSE, seed=1)
#'
#' plot(results)
plot.fqpca_object <- function(x, pve = 0.99, ...)
{
  if (!inherits(x, "fqpca_object")) {stop("The x must be of class fqpca_object")}

  # Determine the number of components required based on cumulative PVE
  n.components <- obtain_npc(scores=x$scores, pve=pve)

  intercept <- x$intercept
  loadings <- x$loadings

  # Create a data frame for the intercept curve
  intercept_df <- data.frame(
    Time = seq_len(length(intercept)),
    Loading = intercept,
    Component = "Intercept")

  non_int <- loadings[, 1:n.components, drop = FALSE]

  # Reshape non-intercept loadings into a long format data frame
  fqpc_df <- data.frame(
    Time = rep(seq_len(nrow(non_int)), times = ncol(non_int)),
    Loading = as.vector(non_int),
    Component = rep(paste0("FQPC ", seq_len(ncol(non_int))), each = nrow(non_int)))

  # Combine the intercept and FQPC data frames
  plot_data <- rbind(intercept_df, fqpc_df)
  plot_data$Component <- factor(plot_data$Component, levels = c('Intercept', paste0('FQPC ', 1:ncol(non_int))))

  # Create the GGPlot object with facets displaying each component separately.
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$Time, y = .data$Loading)) +
    ggplot2::geom_line(color = "steelblue", linewidth = 1) +
    ggplot2::facet_wrap(~ Component, scales = "free_y", ncol = 3) +
    ggplot2::labs(title = "Loading Functions",
                  x = "Time",
                  y = "Loading") +
    ggplot2::theme_bw()
  return(p)
}
