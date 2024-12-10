
# SCORES ----------------------------------------------------------------------

#' @title Compute scores
#' @description Inner function to compute the scores of the fqpca methodology.
#' @param Y The \eqn{(N \times T)} matrix of observed time instants.
#' @param Y.mask Mask matrix of the same dimensions as Y indicating which observations in Y are known.
#' @param loadings Matrix of loading coefficients.
#' @param quantile.value The quantile considered.
#' @return The matrix of scores.
compute_scores_sequential <- function(Y, Y.mask, loadings, quantile.value)
{
  # Initialize the matrix of scores
  n.obs <- base::nrow(Y)
  npc1 <- base::ncol(loadings)
  scores <-  base::matrix(0, n.obs, npc1-1)
  for(i in base::seq(n.obs))
  {
    # Obtain the loadings associated to observation i with a specific time grid
    loadings.obs <- loadings[Y.mask[i, ], ]
    # Obtain the observation removing missing data and subtract intercept
    Yi <- Y[i, Y.mask[i, ]] - loadings.obs[,1]
    loadings.matrix = matrix(loadings.obs[,2:npc1], ncol=npc1-1)
    scores[i,] <-  quantreg::rq.fit.br(y=Yi, x=loadings.matrix, tau=quantile.value)$coefficients
  }
  # Add intercept column
  scores <- base::cbind(1, scores)
  return(scores)
}

#' @title Compute scores
#' @description Inner function to compute the scores of the fqpca methodology.
#' @param Y The \eqn{(N \times T)} matrix of observed time instants.
#' @param Y.mask Mask matrix of the same dimensions as Y indicating which observations in Y are known.
#' @param loadings Matrix of loading coefficients.
#' @param quantile.value The quantile considered.
#' @param num.cores Number of cores used in parallel executions.
#' @return The matrix of scores.
#' @importFrom foreach %dopar%
compute_scores_parallel <- function(Y, Y.mask, loadings, quantile.value, num.cores)
{
  i = 1
  # Initialize the matrix of scores
  n.obs <- base::nrow(Y)
  npc1 <- base::ncol(loadings)
  scores <-  base::matrix(0, n.obs, npc1-1)
  # Initialize the cluster each time the scores are computed
  cl = parallel::makeCluster(num.cores)
  doParallel::registerDoParallel(cl)
  scores <-
    foreach::foreach(i = base::seq(n.obs), .combine='rbind') %dopar%
    {
      loadings.obs <- loadings[Y.mask[i, ], ]
      Yi <- Y[i, Y.mask[i, ]] - loadings.obs[,1]
      loadings.matrix = matrix(loadings.obs[,2:npc1], ncol=npc1-1)
      quantreg::rq.fit.br(y=Yi, x=loadings.matrix, tau=quantile.value)$coefficients
    }
  parallel::stopCluster(cl)
  foreach::registerDoSEQ()
  rownames(scores) = colnames(scores) = NULL
  # Add intercept column
  scores <- base::cbind(1, scores)
  return(scores)
}

#' @title Compute scores
#' @description Inner function to compute the scores of the fqpca methodology.
#' @param Y The \eqn{(N \times T)} matrix of observed time instants.
#' @param Y.mask Mask matrix of the same dimensions as Y indicating which observations in Y are known.
#' @param loadings Matrix of loading coefficients.
#' @param quantile.value The quantile considered.
#' @param parallelized.scores Should the scores be computed in parallel?
#' @param num.cores Number of cores used in parallel executions.
#' @return The matrix of scores.
#' @importFrom foreach %dopar%
compute_scores <- function(Y, Y.mask, loadings, quantile.value, parallelized.scores, num.cores)
{
  if(parallelized.scores)
  {
    scores <- compute_scores_parallel(Y, Y.mask, loadings, quantile.value, num.cores)
  }else{
    scores <- compute_scores_sequential(Y, Y.mask, loadings, quantile.value)
  }
  return(scores)
}

# SPLINES ---------------------------------------------------------------------

#' @title Compute spline coefficients
#' @description Inner function to compute the spline coefficients of the fqpca methodology.
#' @param Y \eqn{(N \times T)} matrix of observed time instants.
#' @param Y.mask Mask matrix of the same dimensions as Y indicating wich observations in Y are known.
#' @param scores Initial matrix of scores.
#' @param spline.basis The spline basis matrix.
#' @param quantile.value The quantile considered.
#' @param method Method used in the resolution of the quantile regression model. It currently accepts the methods \code{c('conquer', 'quantreg')}.
#' @return The matrix of spline coefficients splines coefficients.
compute_spline_coefficients_unpenalized <- function(Y, Y.mask, scores, spline.basis, quantile.value, method)
{
  n.obs <- base::nrow(Y)
  npc1 <- base::ncol(scores)
  n.basis <- base::ncol(spline.basis)
  intercept.spline.basis <-  spline.basis[, -1]
  # vectorize Y [Y11, ..., Y1T, Y21, ..., Y2T,...YNT]
  Y.vector <- c()
  for(i in base::seq(n.obs)){Y.vector <- c(Y.vector, Y[i, Y.mask[i,]])}

  # Compute tensor product of scores and spline basis
  # Solve the large quantile regression model
  if(method == 'conquer')
  {
    tensor.matrix <- base::matrix(0, nrow=base::length(Y.vector), ncol=npc1*n.basis-1)
    row.idx <- 1
    for(i in base::seq(n.obs))
    {
      scores.i <- base::matrix(scores[i, -1], nrow=1)
      tmp.splines <- spline.basis[Y.mask[i, ], ]
      tmp.intercept.splines <- intercept.spline.basis[Y.mask[i, ], ]
      n.time.i <- base::nrow(tmp.splines)
      tensor.matrix[row.idx:(row.idx + n.time.i - 1), ] <- base::cbind(tmp.intercept.splines, base::kronecker(scores.i, tmp.splines))
      row.idx <- row.idx+n.time.i
    }
    B.vector <- conquer::conquer(Y=Y.vector, X=tensor.matrix, tau=quantile.value)$coeff
  }
  if(method == 'quantreg')
  {
    tensor.matrix <- base::matrix(0, nrow=base::length(Y.vector), ncol=npc1*n.basis)
    row.idx <- 1
    for(i in base::seq(n.obs))
    {
      scores.i <- base::matrix(scores[i, ], nrow=1)
      tmp.splines <- spline.basis[Y.mask[i, ], ]
      n.time.i <- base::nrow(tmp.splines)
      tensor.matrix[row.idx:(row.idx + n.time.i - 1), ] <- base::kronecker(scores.i, tmp.splines)
      row.idx <- row.idx+n.time.i
    }
    B.vector <- quantreg::rq.fit.fnb(y=Y.vector, x=tensor.matrix, tau=quantile.value)$coefficients
  }
  spline.coefficients <- base::matrix(B.vector, byrow=FALSE, ncol=npc1)
  return(spline.coefficients)
}

#' @title Compute spline coefficients
#' @description Inner function to compute the spline coefficients of the fqpca methodology.
#' @param Y \eqn{(N \times T)} matrix of observed time instants.
#' @param Y.mask Mask matrix of the same dimensions as Y indicating wich observations in Y are known.
#' @param scores Initial matrix of scores.
#' @param spline.basis The spline basis matrix.
#' @param quantile.value The quantile considered.
#' @param lambda.ridge  Hyper parameter controlling the penalization on the second derivative of the splines. It has effect only with \code{penalized=TRUE} and \code{method='conquer'}.
#' @return The matrix of spline coefficients splines coefficients.
compute_spline_coefficients_penalized <- function(Y, Y.mask, scores, spline.basis, quantile.value, lambda.ridge)
{
  n.obs <- base::nrow(Y)
  npc1 <- base::ncol(scores)
  n.basis <- base::ncol(spline.basis)
  intercept.spline.basis <-  spline.basis[, -1]
  Y.vector <- c()
  for(i in base::seq(n.obs)){Y.vector <- c(Y.vector, Y[i, Y.mask[i,]])}
  tensor.matrix <- base::matrix(0, nrow=base::length(Y.vector), ncol=npc1*n.basis-1)
  row.idx <- 1
  for(i in base::seq(n.obs))
  {
    scores.i <- base::matrix(scores[i, -1], nrow=1)
    tmp.splines <- spline.basis[Y.mask[i, ], ]
    tmp.intercept.splines <- intercept.spline.basis[Y.mask[i, ], ]
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
  spline.coefficients <- base::matrix(B.vector, byrow=FALSE, ncol=npc1)
  return(spline.coefficients)
}

#' @title Compute spline coefficients
#' @description Inner function to compute the spline coefficients of the fqpca methodology.
#' @param penalized Boolean indicating if the smoothness should be controlled using a second derivative penalty.
#' @param Y \eqn{(N \times T)} matrix of observed time instants.
#' @param Y.mask Mask matrix of the same dimensions as Y indicating wich observations in Y are known.
#' @param scores Initial matrix of scores.
#' @param spline.basis The spline basis matrix.
#' @param quantile.value The quantile considered.
#' @param method Method used in the resolution of the quantile regression model. It currently accepts the methods \code{c('conquer', 'quantreg')}.
#' @param lambda.ridge  Hyper parameter controlling the penalization on the second derivative of the splines. It has effect only with \code{penalized=TRUE} and \code{method='conquer'}.
#' @return The matrix of spline coefficients splines coefficients.
compute_spline_coefficients <- function(penalized, Y, Y.mask, scores, spline.basis, quantile.value, method, lambda.ridge)
{
  if(penalized)
  {
    spline.coefficients <- compute_spline_coefficients_penalized(Y=Y, Y.mask=Y.mask, scores=scores, spline.basis=spline.basis, quantile.value=quantile.value, lambda.ridge=lambda.ridge)
  }else{
    spline.coefficients <- compute_spline_coefficients_unpenalized(Y=Y, Y.mask=Y.mask, scores=scores, spline.basis=spline.basis, quantile.value=quantile.value, method=method)
  }
  return(spline.coefficients)
}

# OTHER FUNCTIONS -------------------------------------------------------------

#' @title Compute objective value function.
#' @description Inner function to compute the objective value of the fqpca methodology at each iteration.
#' @param Y \eqn{(N \times T)} matrix of observed time instants.
#' @param quantile.value The quantile considered.
#' @param scores The matrix of estimated scores.
#' @param loadings The matrix of estimated loadings.
#' @return The objective value function.
compute_objective_value <- function(Y, quantile.value, scores, loadings)
{
  Y.pred <- scores %*% t(loadings)
  objective.value <- quantile_error(Y=Y, Y.pred=Y.pred, quantile.value=quantile.value)
  return(objective.value)
}

#' @title Compute loadings (aka principal components)
#' @description Inner function to compute the loadings of the fqpca methodology.
#' @param spline.basis The spline basis matrix.
#' @param spline.coefficients the matrix of spline coefficients.
#' @param method Method used in the resolution of the quantile regression model. It currently accepts the methods \code{c('conquer', 'quantreg')} along with any available solver in \code{CVXR} package.
#' @return The matrix of loadings
compute_loadings <- function(spline.basis, spline.coefficients, method)
{
  intercept.spline.basis <-  spline.basis[, -1]
  if(method == 'conquer')
  {
    intercept.part <- cbind(1, intercept.spline.basis) %*% matrix(spline.coefficients[,1], ncol=1)
    fqpc.part <- spline.basis %*% spline.coefficients[,-1]
    loadings <- cbind(intercept.part, fqpc.part)
  }else{
    loadings <- spline.basis %*% spline.coefficients
  }
  return(loadings)
}

#' @title Rotation of fqpca loadings and scores
#' @description Performs the rotation of loadings and scores in order to ensure the solution is unique.
#' @param loadings Matrix of loadings.
#' @param scores Matrix of scores.
#' @return The rotated matrices of loadings and scores and the rotation matrix.
rotate_scores_and_loadings <- function(loadings, scores)
{
  npc <- base::ncol(loadings)-1

  # Remove intercept from rotation
  loadings.no.intercept <- base::matrix(loadings[, 2:ncol(loadings)], ncol=npc)
  scores.no.intercept <- base::matrix(scores[, 2:ncol(scores)], ncol=npc)

  # Ensure scores are mean-centered and move the effect of the mean scores to the intercept
  means.scores <- base::matrix(colMeans(scores.no.intercept), ncol=1)
  scores.no.intercept <- base::scale(scores.no.intercept, center=TRUE, scale=FALSE)
  loadings[, 1] <- loadings[, 1] + loadings.no.intercept %*% means.scores

  # Perform rotation
  cov.score <- stats::cov(scores.no.intercept)
  svd.decomp <- base::svd((chol(cov.score))  %*% t(loadings.no.intercept))
  rotation.matrix <- base::solve(chol(cov.score)) %*% svd.decomp$u %*% base::diag(svd.decomp$d, nrow=length(svd.decomp$d))

  scores.no.intercept <- scores.no.intercept %*% rotation.matrix
  loadings.no.intercept <- svd.decomp$v

  loadings <- base::cbind(loadings[,1], loadings.no.intercept)
  scores <- base::cbind(scores[,1], scores.no.intercept)
  results <- list(loadings=loadings, scores=scores, rotation.matrix=rotation.matrix)
  return(results)
}

#' @title Explained variability of the fqpca scores computation
#' @description Computes the percentage of explained variability based on the variance of the scores matrix
#' @param scores Matrix of scores.
#' @return The percentage of variability each component is explaining.
compute_explained_variability <- function(scores)
{
  npc <- base::ncol(scores)-1
  scores <- base::matrix(scores[,-1], ncol=npc)
  variability <- base::diag(stats::var(scores))
  pve <- variability / base::sum(variability)
  return(pve)
}

#' @title Check input parameters
#' @description Check input parameters
#' @param npc The number of estimated components.
#' @param quantile.value The quantile considered.
#' @param periodic Boolean indicating if the data is expected to be periodic (start coincides with end) or not.
#' @param splines.df Degrees of freedom for the splines.
#' @param method Method used in the resolution of the quantile regression model. It currently accepts the methods \code{c('conquer', 'quantreg')}.
#' @param penalized Boolean indicating if the smoothness should be controlled using a second derivative penalty. This functionality is experimental and is much slower than the control of the smoothness using the degrees of freedom.
#' @param lambda.ridge  Hyper parameter controlling the penalization on the second derivative of the splines. It has effect only with \code{penalized=TRUE}.
#' @param tol Tolerance on the convergence of the algorithm.
#' @param n.iters Maximum number of iterations.
#' @param verbose Boolean indicating the verbosity.
#' @param seed Seed for the random generator number.
#' @param parallelized.scores Should the scores be computed in parallel? Experimantal component.
#' @param num.cores Number of cores to use in parallelized executions.
#' @return No return
check_fqpca_params <- function(npc, quantile.value, periodic, splines.df, method, penalized, lambda.ridge, tol, n.iters, verbose, seed, parallelized.scores, num.cores)
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

  # Check 'method': must be either 'conquer' or 'quantreg'
  if (!is.character(method) || length(method) != 1 ||
      !method %in% c("conquer", "quantreg")) {
    stop("Invalid input for 'method': '", method,
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

  # Check 'n.iters': positive integer number
  if (!is.numeric(n.iters) || length(n.iters) != 1 ||
      n.iters %% 1 != 0 || n.iters <= 0) {
    stop("Invalid input for 'n.iters': ", n.iters,
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

  # Check 'parallelized.scores': Boolean
  if (!is.logical(parallelized.scores) || length(parallelized.scores) != 1) {
    stop("Invalid input for 'parallelized.scores': ", parallelized.scores,
         ". Expected a Boolean value (TRUE or FALSE).")
  }

  # Check 'num.cores': positive integer number or NULL
  if (!(is.null(num.cores) || (is.numeric(num.cores) && length(num.cores) == 1 &&
                               num.cores %% 1 == 0 && num.cores > 0))) {
    stop("Invalid input for 'num.cores': ", num.cores,
         ". Expected a positive integer number or NULL.")
  }
}

#' @title Check input data
#' @description Check input data
#' @title Explained variability of the fqpca scores computation
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

fqpca_structure <- function(loadings, scores, pve, objective.function.value,
                            objective.function.array, function.warnings,
                            execution.time, colname, npc, quantile.value, method,
                            penalized, lambda.ridge, periodic, splines.df, tol, n.iters, verbose,
                            seed, parallelized.scores, num.cores, rotation.matrix,
                            spline.coefficients, spline.basis)
{
  structure(list(
    loadings = loadings,
    scores = scores,
    pve = pve,
    objective.function.value = objective.function.value,
    objective.function.array = objective.function.array,
    function.warnings = function.warnings,
    execution.time = execution.time,
    colname = colname,
    npc = npc,
    quantile.value = quantile.value,
    method = method,
    penalized = penalized,
    lambda.ridge = lambda.ridge,
    periodic = periodic,
    splines.df = splines.df,
    tol = tol,
    n.iters = n.iters,
    verbose = verbose,
    seed = seed,
    parallelized.scores = parallelized.scores,
    num.cores = num.cores,
    rotation.matrix = rotation.matrix,
    spline.coefficients = spline.coefficients,
    spline.basis = spline.basis
  ),
  class = "fqpca_object")
}

#' @title FQPCA (Functional Quantile Principal Component Analysis)
#' @description Solves the functional quantile principal component analysis methodology
#' @param data An \eqn{(N \times T)} matrix, a tf object from the tidyfun package or a data.frame containing the functional data as a tf column.
#' @param colname The name of the column containing the functional data. Use only if data is a dataframe and colname is a column in the dataframe.
#' @param npc The number of estimated components.
#' @param quantile.value The quantile considered.
#' @param periodic Boolean indicating if the data is expected to be periodic (start coincides with end) or not.
#' @param splines.df Degrees of freedom for the splines.
#' @param method Method used in the resolution of the quantile regression model. It currently accepts the methods \code{c('conquer', 'quantreg')}.
#' @param penalized Boolean indicating if the smoothness should be controlled using a second derivative penalty. This functionality is experimental.
#' @param lambda.ridge  Hyper parameter controlling the penalization on the second derivative of the splines. It has effect only with \code{penalized=TRUE} and \code{method='conquer'}.
#' @param tol Tolerance on the convergence of the algorithm.
#' @param n.iters Maximum number of iterations.
#' @param verbose Boolean indicating the verbosity.
#' @param seed Seed for the random generator number.
#' @param parallelized.scores Should the scores be computed in parallel? Experimantal component.
#' @param num.cores Number of cores to use in parallelized executions.
#' @return fqpca_object
#' @export
#' @examples
#' # Generate fake dataset with 150 observations and 144 time points
#'
#' Y <- matrix(rep(sin(seq(0, 2*pi, length.out = 144)), 150), byrow = TRUE, nrow = 150)
#' Y <- Y + matrix(rnorm(150*144, 0, 0.4), nrow = 150)
#'
#' # Add missing observations
#' Y[sample(150*144, as.integer(0.2*150*144))] <- NA
#'
#' results <- fqpca(data = Y, npc = 2, quantile.value = 0.5)
#'
#' loadings <- results$loadings
#' scores <- results$scores
fqpca <- function(data, colname = NULL, npc = 2,  quantile.value = 0.5,  periodic = TRUE, splines.df = 10, method = 'conquer', penalized = FALSE, lambda.ridge = 0, tol = 1e-3, n.iters = 20, verbose = FALSE, seed = NULL, parallelized.scores=FALSE, num.cores=NULL)
{
  global.start.time <- base::Sys.time()

   # Check the input parameters except Y and colname
  check_fqpca_params(npc=npc, quantile.value=quantile.value, periodic=periodic,
                     splines.df=splines.df, method=method, penalized=penalized,
                     lambda.ridge=lambda.ridge, tol=tol, n.iters=n.iters,
                     verbose=verbose, seed=seed, parallelized.scores=parallelized.scores,
                     num.cores=num.cores)

  # Check Y and colname and return an unnamed matrix
  Y <- check_input_data(data=data, colname=colname)

  # Check if method is valid
  if(penalized && method != 'conquer')
  {
    stop("Invalid input for 'method and/or penalized':
         If penalized is set to  TRUE then method must be 'conquer'")
  }

  # If seed is provided, set seed for computations
  if(!base::is.null(seed)){base::set.seed(seed)}

  # If scores are computed in parallel, initialize num.cores
  if(parallelized.scores && is.null(num.cores)){num.cores <- parallel::detectCores() - 1}

  # Step 1: Initialize values
  n.time <- base::dim(Y)[2]
  Y.axis <- base::seq(0, 1, length.out = n.time)
  Y.mask <- !base::is.na(Y)
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
  n.basis = base::ncol(spline.basis)

  # If penalized, compute different basis with an identity second derivative matrix with 0s in first and second positions
  # Based on Garcia de la Garza et al. (2024) (Adaptive Functional Principal Component Analysis)
  if(penalized)
  {
    eigen.decomp = base::eigen(base::crossprod(base::diff(spline.basis, differences = 2)))
    first.part = eigen.decomp$vectors[, (n.basis-1):n.basis]
    second.part = eigen.decomp$vectors[, 1:(n.basis-2)] %*% (base::diag(1 / base::sqrt(eigen.decomp$values[1:(n.basis-2)])))
    U = cbind(first.part, second.part)
    spline.basis = spline.basis %*% U
  }

  loop.start.time <- base::Sys.time()

  # Randomly generate splines coefficients
  spline.coefficients <- base::matrix(stats::rnorm(n.basis * (npc+1), mean = 0, sd = 1), nrow=n.basis, ncol=npc+1)
  loadings <- spline.basis %*% spline.coefficients

  scores <- try(compute_scores(Y = Y, Y.mask = Y.mask, loadings = loadings, quantile.value = quantile.value, parallelized.scores = parallelized.scores, num.cores = num.cores), silent = FALSE)
  if(!is.matrix(scores)){stop('Iteration 1. Failed computation of scores')}

  # Rotate loadings and scores
  rotation.result <- try(rotate_scores_and_loadings(loadings = loadings, scores = scores), silent = FALSE)
  error.chr <- class(rotation.result)
  if(error.chr=="try-error")
  {
    warning('Iteration 1. Failed rotation process. Skipping rotation on this iteration.')
    function.warnings$rotation <- TRUE
    rotation.result <- list(loadings = loadings, scores = scores, rotation.matrix = diag(npc))
  }
  loadings <- rotation.result$loadings
  scores <- rotation.result$scores
  rotation.matrix = rotation.result$rotation.matrix

  # Compute objective value function
  objective.function <- compute_objective_value(quantile.value = quantile.value, Y = Y, scores = scores, loadings = loadings)
  objective.function.array <- objective.function

  last.iteration <- list(loadings = loadings, scores = scores, spline.coefficients = spline.coefficients, objective.function=objective.function, rotation.matrix = rotation.matrix, iteration = 1)

  loop.execution.time <- difftime(base::Sys.time(), loop.start.time, units = 'secs')

  if(verbose)
  {
    message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '. Iteration ', 1, ' completed in ', base::round(loop.execution.time, 3), ' seconds',
            '\n', 'Objective function   i   value: ', round(objective.function.array[1], 4),
            '\n', '____________________________________________________________')
  }

  # Iterate until convergence or until n.iters is reached
  convergence = F
  for(i in 2:n.iters)
  {
    loop.start.time <- base::Sys.time()

    # OBTAIN SPLINE COEFFICIENTS
    spline.coefficients <- try(compute_spline_coefficients(penalized = penalized, Y = Y, Y.mask = Y.mask, scores = scores, spline.basis = spline.basis, quantile.value = quantile.value, method = method, lambda.ridge = lambda.ridge), silent = FALSE)
    if(!is.matrix(spline.coefficients))
    {
      warning('Iteration ', i, '. Failed computation of spline coefficients. Providing results from previous iteration.')
      function.warnings$splines <- TRUE
      break
    }

    loadings <- compute_loadings(spline.basis=spline.basis, spline.coefficients=spline.coefficients, method=method)
    if(is_rank_deficient(loadings)){warning('Loadings matrix is singular.')}

    # OBTAIN SCORES
    scores <- try(compute_scores(Y = Y, Y.mask = Y.mask, loadings = loadings, quantile.value = quantile.value, parallelized.scores = parallelized.scores, num.cores = num.cores), silent = FALSE)
    if(!is.matrix(scores))
    {
      warning('Iteration ', i, '. Failed computation of scores. Providing results from previous iteration.')
      function.warnings$scores <- TRUE
      break
    }

    # Rotate loadings and scores
    rotation.result <- try(rotate_scores_and_loadings(loadings = loadings, scores = scores), silent = FALSE)
    error.chr <- class(rotation.result)
    if(error.chr=="try-error")
    {
      warning('Iteration ', i, '. Failed rotation process. Skipping rotation on this iteration.')
      function.warnings$rotation <- TRUE
      rotation.result <- list(loadings = loadings, scores = scores, rotation.matrix = diag(npc))
    }
    loadings <- rotation.result$loadings
    scores <- rotation.result$scores
    rotation.matrix = rotation.result$rotation.matrix

    # Compute objective function
    objective.function <- compute_objective_value(quantile.value = quantile.value, Y = Y, scores = scores, loadings = loadings)
    objective.function.array <- c(objective.function.array, objective.function)

    last.iteration <- list(loadings = loadings, scores = scores, spline.coefficients = spline.coefficients, objective.function=objective.function, rotation.matrix = rotation.matrix, iteration = i)

    convergence.criteria <- base::abs(objective.function.array[i] - objective.function.array[i-1])
    if(convergence.criteria < tol){convergence = T}

    # Measure computation time
    loop.execution.time <- difftime(base::Sys.time(), loop.start.time, units = 'secs')

    if(verbose & convergence)
    {
      message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '. Iteration ', i, ' completed in ', base::round(loop.execution.time, 3), ' seconds',
              '\n', 'Objective function (i-1) value: ', base::round(objective.function.array[i-1], 4),
              '\n', 'Objective function   i   value: ', round(objective.function.array[i], 4),
              '\n', 'Convergence criteria value    : ', base::round(convergence.criteria, 4),
              '\n', 'Algorithm converged succesfully.',
              '\n', '____________________________________________________________',
              '\n\n')
    }
    if(verbose & !convergence)
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
  if(i == n.iters & !convergence){warning('Algorithm reached maximum number of iterations without convergence. Consider increasing the value of n.iters parameter')}

  # EXPLAINED VARIABILITY -----------------------------------------------------

  pve <- compute_explained_variability(last.iteration$scores)

  global.execution.time <- difftime(base::Sys.time(), global.start.time, units = 'secs')

  results <- fqpca_structure(
    loadings = last.iteration$loadings,
    scores = last.iteration$scores,
    pve = pve,
    objective.function.value = last.iteration$objective.function,
    objective.function.array = objective.function.array,
    function.warnings = function.warnings,
    execution.time = global.execution.time,
    npc = npc,
    colname = colname,
    quantile.value = quantile.value,
    method = method,
    penalized = penalized,
    lambda.ridge = lambda.ridge,
    periodic = periodic,
    splines.df = splines.df,
    tol = tol,
    n.iters = last.iteration$iteration,
    verbose = verbose,
    seed = seed,
    parallelized.scores = parallelized.scores,
    num.cores = num.cores,
    rotation.matrix = last.iteration$rotation.matrix,
    spline.coefficients = last.iteration$spline.coefficients,
    spline.basis = spline.basis
  )
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
#' # Generate fake dataset with 150 observations and 144 time points
#'
#' Y <- matrix(rep(sin(seq(0, 2*pi, length.out = 144)), 150), byrow = TRUE, nrow = 150)
#' Y <- Y + matrix(rnorm(150*144, 0, 0.4), nrow = 150)
#'
#' # Add missing observations
#' Y[sample(150*144, as.integer(0.2*150*144))] <- NA
#'
#' results <- fqpca(data = Y[1:100,], npc = 2, quantile.value = 0.5)
#'
#' predictions <- predict(object = results, newdata = Y[101:150,])
predict.fqpca_object <- function(object, newdata, ...)
{
  if (!base::inherits(object, "fqpca_object")){stop('The object must be of class fqpca_object')}

  if(is.null(object$colname))
  {
    Y <- base::unname(base::as.matrix(newdata))
  }else{
    Y <- base::unname(base::as.matrix(dplyr::pull(newdata[object$colname])))
  }
  if(base::ncol(newdata) == 1){ newdata <- t(newdata)}

  Y.mask <- !base::is.na(newdata)
  n.obs <- base::nrow(newdata)
  n.time <- base::ncol(newdata)
  if(sum(!Y.mask) == n.obs * n.time){stop('newdata contains no information. Check that there are values different from NA')}

  # Estimation of scores
  scores <- try(compute_scores(Y = newdata, Y.mask = Y.mask, loadings = object$loadings, quantile.value = object$quantile.value, parallelized.scores = FALSE, num.cores = 1), silent = FALSE)
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
#' # Generate fake dataset with 150 observations and 144 time points
#'
#' Y <- matrix(rep(sin(seq(0, 2*pi, length.out = 144)), 150), byrow = TRUE, nrow = 150)
#' Y <- Y + matrix(rnorm(150*144, 0, 0.4), nrow = 150)
#'
#' # Add missing observations
#' Y[sample(150*144, as.integer(0.2*150*144))] <- NA
#'
#' results <- fqpca(data = Y[1:100,], npc = 2, quantile.value = 0.5)
#'
#' Yhat <- fitted(object = results, pve=0.99)
fitted.fqpca_object <- function(object, pve=0.95, ...)
{
  if (!base::inherits(object, "fqpca_object")){stop('The object must be of class fqpca_object')}
  if(pve < 1)
  {
    score.variance = cumsum(object$pve)
    n.components = min(which(score.variance > pve))
  }else{
    if(!pve == floor(pve)){stop('pve must be either a floating point number smaller than 1 or an integer number. Value provided: ', pve)}
    n.components = pve
  }
  Y.reconstruction = object$scores[, 1:(n.components+1)] %*% t(object$loadings[, 1:(n.components+1)])
  return(Y.reconstruction)
}

# BASIC PLOT ------------------------------------------------------------------

#' @title Plot fqpca loading functions
#' @description S3 method for class 'fqpca_object'. Given a fqpca object, plot the loading functions
#' @param x An object output of the fqpca function.
#' @param pve Percentage of explained variability plotted.
#' @param ... further arguments passed to or from other methods.
#' @return The plot of loadings
#' @export
#' @examples
#' # Generate fake dataset with 150 observations and 144 time points
#'
#' Y <- matrix(rep(sin(seq(0, 2*pi, length.out = 144)), 150), byrow = TRUE, nrow = 150)
#' Y <- Y + matrix(rnorm(150*144, 0, 0.4), nrow = 150)
#'
#' # Add missing observations
#' Y[sample(150*144, as.integer(0.2*150*144))] <- NA
#'
#' results <- fqpca(data = Y, npc = 2, quantile.value = 0.5)
#'
#' plot(results)
plot.fqpca_object <- function(x, pve=0.99, ...)
{
  if (!base::inherits(x, "fqpca_object")){stop('The x must be of class fqpca_object')}

  score.variance = cumsum(x$pve)
  n.components = min(which(score.variance > pve))

  loadings <- x$loadings
  graphics::par(mfrow = c(1, 2))
  graphics::plot(loadings[,1], type = 'l', lty = 1, lwd = 2, xlab = '', ylab = '')
  graphics::title('Intercept curve')

  no_intercept_loadings <- loadings[,2:(n.components+1)]
  if(base::ncol(loadings) == 2)
  {
    no_intercept_loadings <- base::matrix(no_intercept_loadings, ncol = 1)
  }
  graphics::matplot(no_intercept_loadings, type =  "l", lty = 1, lwd = 2, xlab = '', ylab = '')
  graphics::title('FQPCs ')
}
