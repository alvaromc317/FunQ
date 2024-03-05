# FUNCTIONAL QUANTILE PCA ALGORITHM -------------------------------------------

# INNER LOOP FUNCTIONS --------------------------------------------------------

#' Compute objective value function
#'
#' Inner function to compute the objective value
#'
#' @param Y \eqn{(N \times T)} matrix of observed time instants.
#' @param quantile.value The quantile considered.
#' @param scores The matrix of estimated scores.
#' @param loadings The matrix of estimated loadings.
#'
#' @return The objective value function
compute_objective_value <- function(Y, quantile.value, scores, loadings)
{
  Y.pred <- scores %*% t(loadings)
  objective.value <- FunQuantPCA::quantile_error(Y=Y, Y.pred=Y.pred, quantile.value=quantile.value)
  return(objective.value)
}

#' Compute splines coefficients
#'
#' Inner function to compute the splines coefficients of the fqpca methodology.
#'
#' @param Y \eqn{(N \times T)} matrix of observed time instants.
#' @param Y.mask Mask matrix of the same dimensions as Y indicating wich observations in Y are known.
#' @param scores Initial matrix of scores.
#' @param spline.basis The spline basis matrix.
#' @param quantile.value The quantile considered.
#' @param method Method used in the resolution of the quantile regression model.
#'                It currently accepts the methods \code{c('br', 'fn', 'pfn', 'sfn')}
#'                from \code{quantreg} package along with any available solver in \code{CVXR} package.
#' @param alpha.ridge The hyper-parameter controlling the penalization on the splines.
#'                    This parameter has no effect if method is not a cvxr-valid
#'                    method.
#' @param R.block Block diagonal matrix made of quadratic matrices used to apply
#'                a ridge based penalty.  This only has effect if method is a
#'                cvxr-valid method.
#' @param valid.in.quantreg Array of valid \code{quantreg} methods.
#'
#' @return The matrix of spline coefficients splines coefficients.
compute_spline_coefficients <- function(Y, Y.mask, scores, spline.basis, quantile.value, method, alpha.ridge, R.block, valid.in.quantreg)
{
  # At each iteration, call this function to obtain the estimation of F and Lambda
  n.obs <- base::nrow(Y)
  npc <- base::ncol(scores)
  n.basis <- base::ncol(spline.basis)
  # vectorize Y [Y11, ..., Y1T, Y21, ..., Y2T,...YNT]
  Y.vector <- c()
  for(i in base::seq(n.obs))
  {
    Y.vector <- c(Y.vector, Y[i, Y.mask[i,]])
  }
  # Compute tensor product of Lambda and spline basis
  tensor.matrix <- base::matrix(0, nrow=base::length(Y.vector), ncol=npc*n.basis)
  row.idx <- 1
  for(i in base::seq(n.obs))
  {
    lambda.i <- base::matrix(scores[i, ], nrow=1)
    tmp.splines <- spline.basis[Y.mask[i, ], ]
    n.time.i <- base::nrow(tmp.splines)
    tensor.matrix[row.idx:(row.idx + n.time.i - 1), ] <- base::kronecker(lambda.i, tmp.splines)
    row.idx <- row.idx+n.time.i
  }

  if(method %in% valid.in.quantreg)
  {
    B.vector <- quantreg::rq(Y.vector ~ -1+tensor.matrix, tau=quantile.value, method=method)$coefficients
  } else{
    B.vector <- quantile_regression_ridge(x=tensor.matrix, y=Y.vector,  R=R.block, quantile.value=quantile.value, lambda=alpha.ridge, solver=method)
  }
  spline.coefficients <- base::matrix(B.vector, byrow=FALSE, ncol=npc)
  return(spline.coefficients)
}

#' Compute loadings (aka principal components)
#'
#' Inner function to compute the loadings of the fqpca methodology.
#'
#' @param spline.basis The spline basis matrix.
#' @param spline.coefficients the matrix of spline coefficients.
#'
#' @return The matrix of loadings
compute_loadings <- function(spline.basis, spline.coefficients)
{
  loadings <- spline.basis %*% spline.coefficients
  return(loadings)
}

#' Compute scores
#'
#' Inner function to compute the scores of the fqpca methodology.
#'
#' @param Y The \eqn{(N \times T)} matrix of observed time instants.
#' @param Y.mask Mask matrix of the same dimensions as Y indicating which observations in Y are known.
#' @param loadings Matrix of loading coefficients.
#' @param quantile.value The quantile considered.
#' @param parallelized.scores Should the scores be computed in parallel?
#' @param num.cores Number of cores used in parallel executions.
#'
#' @return The matrix of scores.
#' @importFrom foreach %dopar%
#'
compute_scores <- function(Y, Y.mask, loadings, quantile.value, parallelized.scores, num.cores)
{
  n.obs <- base::nrow(Y)
  npc <- base::ncol(loadings)
  # Initialize the matrix of scores
  scores <-  base::matrix(0, n.obs, npc-1)

  if(parallelized.scores)
  {
    cl = parallel::makeCluster(num.cores)
    doParallel::registerDoParallel(cl)
    scores <-
      foreach::foreach(i = base::seq(n.obs), .combine='rbind') %dopar%
      {
        loadings.obs <- loadings[Y.mask[i, ], ]
        Yi <- Y[i, Y.mask[i, ]]
        Yi <- Yi - loadings.obs[,1]
        quantreg::rq(Yi ~ -1 + loadings.obs[,2:npc], tau=quantile.value)$coefficients
      }
    parallel::stopCluster(cl)
    foreach::registerDoSEQ()
    rownames(scores) = colnames(scores) = NULL
  } else{
    for(i in base::seq(n.obs))
    {
      # Obtain the loadings associated to observation i with a specific time grid
      loadings.obs <- loadings[Y.mask[i, ], ]
      # Obtain the observation removing missing data
      Yi <- Y[i, Y.mask[i, ]]
      # Substract intercept
      Yi <- Yi - loadings.obs[,1]
      scores[i,] <-  quantreg::rq(Yi ~ -1 + loadings.obs[,2:npc], tau=quantile.value)$coefficients
    }
  }
  # Add intercept column
  scores <- base::cbind(1, scores)
  return(scores)
}

# OUTER LOOP FUNCTIONS --------------------------------------------------------

#' FQPCA rotation
#'
#' Performs the rotation of matrices of loadings and scores in order to ensure the solution is unique.
#'
#' @param loadings Matrix of loadings.
#' @param scores Matrix of scores.
#'
#' @return The rotated matrices of loadings and scores and the rotation matrix.
rotate_scores_and_loadings <- function(loadings, scores)
{
  npc <- base::ncol(loadings)-1

  # Remove intercept from rotation
  loadings.no.intercept <- base::matrix(loadings[, 2:ncol(loadings)], ncol=npc)
  scores.no.intercept <- base::matrix(scores[, 2:ncol(scores)], ncol=npc)

  # Ensure scores are mean-centered and move the effect of the mean scores to the intercept
  means_scores <- base::matrix(colMeans(scores.no.intercept), ncol=1)
  scores.no.intercept <- base::scale(scores.no.intercept, center=TRUE, scale=FALSE)
  loadings[, 1] <- loadings[, 1] + loadings.no.intercept %*% means_scores

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

#' FQPCA explained variability
#'
#' Computes the percentage of explained variability based on the variance of the scores matrix
#'
#' @param scores Matrix of scores.
#'
#' @return The percentage of variability each component is explaining.
compute_explained_variability <- function(scores)
{
  if(base::ncol(scores) == 2)
  {
    scores <- base::matrix(scores[,-1], ncol=1)
  } else{
    scores <- scores[,-1]
  }
  variability <- base::diag(stats::var(scores))
  pve <- variability / base::sum(variability)
  return(pve)
}

# MAIN ------------------------------------------------------------------------

new_fqpca <- function(loadings, scores, pve, objective.function.value,
                      list.objective.function.values, execution.time, function.warnings,
                      colname, npc, quantile.value, method, alpha.ridge,
                      periodic, splines.df, tol, n.iters, verbose, seed,
                      parallelized.scores, num.cores,
                      rotation.matrix, spline.coefficients, spline.basis)
{
  structure(list(
    loadings = loadings,
    scores = scores,
    pve = pve,
    objective.function.value = objective.function.value,
    list.objective.function.values = list.objective.function.values,
    execution.time = execution.time,
    function.warnings = function.warnings,
    colname = colname,
    npc = npc,
    quantile.value = quantile.value,
    method = method,
    alpha.ridge = alpha.ridge,
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


#' FQPCA
#'
#' Solves the functional quantile principal component analysis methodology
#'
#' @param Y An \eqn{(N \times T)} matrix or a tf object from the tidyfun package.
#' @param data data.frame containing the functional data as a tf column.
#' @param colname The name of the column containing the functional data. If used, data argument must not be NULL.
#' @param npc The number of estimated components.
#' @param quantile.value The quantile considered.
#' @param periodic Boolean indicating if the data is expected to
#'                  be periodic (start coincides with end) or not.
#' @param splines.df Degrees of freedom for the splines.
#' @param method Method used in the resolution of the quantile
#'                regression model. It currently accepts the methods
#'                \code{c('br', 'fn', 'pfn', 'sfn')} from \code{quantreg} package along with
#'                any available solver in \code{CVXR} package.
#' @param alpha.ridge  Hyper parameter controlling the penalization
#'                      on the second derivative of the splines. It has effect
#'                      only with \code{CVXR} methods. Experimantal component.
#' @param tol Tolerance on the convergence of the algorithm.
#' @param n.iters Maximum number of iterations.
#' @param verbose Boolean indicating the verbosity.
#' @param seed Seed for the random generator number.
#' @param parallelized.scores Should the scores be computed in parallel? Experimantal component.
#' @param num.cores Number of cores to use in parallelized executions.
#'
#' @return fqpca_object
#' @export
#'
#' @examples
#' # Generate fake dataset with 150 observations and 144 time points
#'
#' Y <- matrix(rep(sin(seq(0, 2*pi, length.out = 144)), 150), byrow = TRUE, nrow = 150)
#' Y <- Y + matrix(rnorm(150*144, 0, 0.4), nrow = 150)
#'
#' # Add missing observations
#' Y[sample(150*144, as.integer(0.2*150*144))] <- NA
#'
#' results <- fqpca(Y = Y, npc = 2, quantile.value = 0.5)
#'
#' loadings <- results$loadings
#' scores <- results$scores
fqpca <- function(Y = NULL, data = NULL, colname = NULL, npc = 2,  quantile.value = 0.5,  periodic = TRUE, splines.df = 10, method = 'fn', alpha.ridge = 0, tol = 1e-3, n.iters = 20, verbose = FALSE, seed = NULL, parallelized.scores=FALSE, num.cores=NULL)
{
  global_start_time <- base::Sys.time()

  # BASIC COMPROBATIONS -------------------------------------------------------

  if(!npc == base::floor(npc)){stop('npc must be an integer number. Value provided: ', npc)}
  if(quantile.value<0 || quantile.value>1){stop('quantile.value must be a value between 0 and 1. Value provided: ', quantile.value)}
  if(!splines.df == base::floor(splines.df)){stop('splines.df must be an integer number. Value provided: ', splines.df)}
  if(!n.iters == base::floor(n.iters)){stop('n.iters must be an integer number. Value provided: ', n.iters)}
  if(alpha.ridge<0){stop('alpha.ridge must be greater or equal than 0')}
  if(tol<0){stop('tol must be greater or equal than 0')}
  if(!n.iters == base::floor(n.iters)){stop('n.iters must be an integer number. Value provided: ', n.iters)}
  if(npc > splines.df){stop('The number of componets cannot be larger than the degrees of freedom.')}

  # DETERMINE IF THE DATA PROVIDED IS A MATRIX OR FUNCTIONAL VECTOR Y, OR IF IT IS A DATAFRAME-COLUMN PAIR
  if(!is.null(Y))
  {
    if(!is.matrix(Y))
    {
      # This structure allows tf to be a suggested package rather than mandatory
      if(!tf::is_tf(Y))
      {
        stop('Y is not of a valid type. Object provided: ', class(Y))
      }
      Y <- base::unname(base::as.matrix(Y))
      colname <- NULL
    }
  }else{
    if(!is.null(data) && !is.null(colname))
    {
      if(!is.data.frame(data))
      {
        stop('data is not of a valid type. Object provided: ', class(data))
      }
      if(!is.character(colname))
      {
        stop('colname is not of a valid type. Object provided: ', class(colname))
      }
      Y <- dplyr::pull(data[colname])
      if(!tf::is_tf(Y))
      {
        stop('The colname indicated by parameter colname in the dataframe data must be a tidyfun vector')
      }
      Y <- base::unname(base::as.matrix(Y))
    }else{
      stop('Either parameter Y or parameters data and colname must be informed.')
    }
  }

  # DETERMINE IF THE ALGORITHM RUNS USING PENALIZED OR UNPENALIZED SPLINES ----

  valid.in.quantreg <- c('br', 'fn', 'sfn', 'pfn')
  if(method %in% valid.in.quantreg)
  {
    penalized <- FALSE
  } else{
    # Check available solvers in CVXR
    valid.in.cvxr <- CVXR::installed_solvers()
    if(method %in% valid.in.cvxr)
    {
      penalized <- TRUE
    } else{
      stop('Invalid method provided. Valid methods include the quantreg methods
         c("br", "fn", "sfn", "pfn") and any installed solver for CVXR.
         Value provided: ', method)
    }
  }

  if(!base::is.null(seed)){base::set.seed(seed)}

  if(parallelized.scores)
  {
    if(is.null(num.cores))
    {
      num.cores <- parallel::detectCores() - 1
    }
  }

  # Step 1: Initialize values
  Y <- base::unname(base::as.matrix(Y))
  n.time <- base::dim(Y)[2]
  Y.axis <- base::seq(0, 1, length.out = n.time)
  Y.mask <- !base::is.na(Y)
  function.warnings <- list(splines = FALSE, scores = FALSE, diverging.loop = FALSE, rotation = FALSE)

  # INITIALIZATION OF SPLINE BASIS --------------------------------------------

  # Initialize splines and second derivative ---
  if(periodic)
  {
    knots <- base::seq(0, 1, length.out = 2 + splines.df - 1)  # 2 boundary knots + df - intercept
  } else{
    knots <- base::seq(0, 1, length.out = 2 + splines.df - 3 - 1) # 2 boundary knots + df - degree - intercept
  }
  knots <- knots[2:(base::length(knots)-1)]
  spline.basis <-  splines2::mSpline(Y.axis, degree = 3, intercept = TRUE, knots = knots, periodic = periodic, Boundary.knots = c(min(Y.axis), max(Y.axis)))

  if(penalized)
  {
    deriv.basis <- stats::deriv(spline.basis, derivs = 2)
    R <- t(deriv.basis) %*% deriv.basis
    R.block <- simex::diag.block(R, npc + 1)
    # Ensure that R.block is semi definite positive: (all eigenvalues are greater or equal to 0)
    for(i in 1:5000)
    {
      R.block <- (1 - 1e-3) * R.block + 1e-3 * base::diag(base::nrow(R.block))
      if(all(base::eigen(R.block)$values > 0)){break}
    }
    if(!all(base::eigen(R.block)$values > 0)) {stop('Splines penalization matrix is not positive definite')}
  } else{R.block <- NULL}

  # INITIALIZATION OF SCORES AND LOADINGS -------------------------------------

  # Generate initial splines basis coefficients
  spline_coefficients_0 <- base::matrix(stats::rnorm(base::dim(spline.basis)[2] * (npc+1), mean = 0, sd = 1),
                                        nrow = base::dim(spline.basis)[2], ncol = npc+1)

  # Initialize scores and loadings with these spline coefficients
  loadings_0 <- compute_loadings(spline.basis = spline.basis, spline.coefficients = spline_coefficients_0)
  scores_0 <- try(compute_scores(Y = Y, Y.mask = Y.mask, loadings = loadings_0, quantile.value = quantile.value, parallelized.scores = parallelized.scores, num.cores = num.cores), silent = FALSE)
  if(!is.matrix(scores_0))
  {
    stop('Computation of scores failed at iteration ', 0)
  }

  # Rotate loadings and scores
  rotation_result <- try(rotate_scores_and_loadings(loadings = loadings_0, scores = scores_0), silent = TRUE)
  error_chr <- class(rotation_result)
  if(error_chr=="try-error")
  {
    warning('Rotation step failed at iteration: ', 0, '. Skipping rotation on this iteration.')
    function.warnings$rotation <- TRUE
    rotation_result <- list(loadings = loadings_0, scores = scores_0, rotation.matrix = diag(npc))
  }
  loadings_0 <- rotation_result$loadings
  scores_0 <- rotation_result$scores
  rotation.matrix_0 = rotation_result$rotation.matrix

  # Compute objective value function
  objective_function_0 <- compute_objective_value(quantile.value = quantile.value, Y = Y, scores = scores_0, loadings = loadings_0)

  # Initialize storage of final solution
  objective_function_array <- objective_function_0
  best_results <- list(loadings = loadings_0, scores = scores_0, objective_function = objective_function_0, spline.coefficients = spline_coefficients_0, rotation.matrix = rotation.matrix_0, iteration = 0)

  # Iterate until convergence or until max_iters is reached
  convergence_criteria <- c()

  for(i in 1:n.iters)
  {
    loop_start_time <- base::Sys.time()

    # COMPUTATION OF LOADINGS AND SCORES ---------------------------

    # Obtain splines coefficients
    spline_coefficients_1 <- try(compute_spline_coefficients(Y = Y, Y.mask = Y.mask, scores = scores_0, spline.basis = spline.basis, quantile.value = quantile.value, method = method, alpha.ridge = alpha.ridge, R.block = R.block, valid.in.quantreg = valid.in.quantreg), silent = TRUE)
    if(!is.matrix(spline_coefficients_1))
    {
      warning('Computation of spline coefficients failed at iteration ', i, '. Providing results from previous iteration.')
      function.warnings$splines <- TRUE
      break
    }

    # Compute loadings
    loadings_1 <- compute_loadings(spline.basis = spline.basis, spline.coefficients = spline_coefficients_1)

    # Compute scores
    scores_1 <- try(compute_scores(Y = Y, Y.mask = Y.mask, loadings = loadings_1, quantile.value = quantile.value, parallelized.scores = parallelized.scores, num.cores = num.cores), silent = FALSE)
    if(!is.matrix(scores_1))
    {
      warning('Computation of scores failed at iteration ', i, '. Providing results from previous iteration.')
      function.warnings$scores <- TRUE
      break
    }

    # Rotate loadings and scores
    rotation_result <- try(rotate_scores_and_loadings(loadings = loadings_1, scores = scores_1), silent = TRUE)
    error_chr <- class(rotation_result)
    if(error_chr=="try-error")
    {
      warning('Rotation step failed at iteration: ', i, '. Skipping rotation on this iteration.')
      function.warnings$rotation <- TRUE
      rotation_result <- list(loadings = loadings_1, scores = scores_1, rotation.matrix = diag(npc))
    }
    loadings_1 <- rotation_result$loadings
    scores_1 <- rotation_result$scores
    rotation.matrix_1 = rotation_result$rotation.matrix

    # Compute objective function
    objective_function_1 <- compute_objective_value(quantile.value = quantile.value, Y = Y, scores = scores_1, loadings = loadings_1)
    objective_function_array <- c(objective_function_array, objective_function_1)

    # COMPROBATION OF CONVERGENCE ---------------------------------------------

    if(objective_function_1 < best_results$objective_function)
    {
      best_results <- list(loadings = loadings_1, scores = scores_1, objective_function = objective_function_1, spline.coefficients = spline_coefficients_1, rotation.matrix = rotation.matrix_1, iteration = i)
    }

    convergence_criteria[i] <- base::abs(objective_function_1 - objective_function_0)

    # Measure computation time
    loop_end_time <- base::Sys.time()
    loop_execution_time <- difftime(loop_end_time, loop_start_time, units = 'secs')

    if(verbose)
    {
      message('Iteration ', i, ' completed in ', base::round(loop_execution_time, 3), ' seconds',
              '\n', 'Convergence criteria value: ', base::round(convergence_criteria[i], 4),
              '\n', 'Objective function (i-1) value: ', base::round(objective_function_0, 4), ' Objective function i value: ', round(objective_function_1, 4))
      message('___________________________________________________')
    }

    if(convergence_criteria[i] < tol)
    {
      if(verbose){message('Algorithm converged with value: ', base::round(convergence_criteria[i], 4))}
      break
    }

    # If the objective value function (that should be decreasing) is
    # 100 times larger than the original objective value function, break
    if( objective_function_array[i+1] > 100 * objective_function_array[2])
    {
      function.warnings$diverging.loop <- TRUE
      if(verbose)
      {
        message('Breaking diverging loop at iteration ', i,
                '\n', 'Objective function value: ', base::round(objective_function_array[i+1], 4),
                '\n', 'Initial objective value function: ', base::round(objective_function_array[2], 4))
        message('Results of the iteration with the smallest objective value: ', base::round(best_results$objective_function, 4), ' will be provided')
      }
      break
    }

    # UPDATE VALUES -----------------------------------------------------------

    objective_function_0 <- objective_function_1
    scores_0 <- scores_1
  }
  if(i == n.iters){warning('Algorithm reached maximum number of iterations without convergence: ', n.iters, ' iterations')}

  # EXPLAINED VARIABILITY -----------------------------------------------------

  pve <- compute_explained_variability(rotation_result$scores)

  global_end_time <- base::Sys.time()
  global_execution_time <- difftime(global_end_time, global_start_time, units = 'secs')

  results <- new_fqpca(
    loadings = best_results$loadings,
    scores = best_results$scores,
    pve = pve,
    objective.function.value = best_results$objective_function,
    list.objective.function.values = objective_function_array,
    execution.time = global_execution_time,
    function.warnings = function.warnings,
    npc = npc,
    colname = colname,
    quantile.value = quantile.value,
    method = method,
    alpha.ridge = alpha.ridge,
    periodic = periodic,
    splines.df = splines.df,
    tol = tol,
    n.iters = best_results$iteration,
    verbose = verbose,
    seed = seed,
    parallelized.scores = parallelized.scores,
    num.cores = num.cores,
    rotation.matrix = best_results$rotation.matrix,
    spline.coefficients = best_results$spline.coefficients,
    spline.basis = spline.basis
  )
  return(results)
}

# PREDICTIONS -----------------------------------------------------------------

#' Predict FQPCA
#'
#' ## S3 method for class 'fqpca_object'
#' Given a new matrix Y, predicts the value of the scores associated to the given matrix
#'
#' @param object An object output of the fqpca function.
#' @param newdata The N by T matrix of observed time instants to be tested.
#' @param ... further arguments passed to or from other methods.
#'
#' @return The normalized matrix of scores.
#' @export
#'
#' @examples
#' # Generate fake dataset with 150 observations and 144 time points
#'
#' Y <- matrix(rep(sin(seq(0, 2*pi, length.out = 144)), 150), byrow = TRUE, nrow = 150)
#' Y <- Y + matrix(rnorm(150*144, 0, 0.4), nrow = 150)
#'
#' # Add missing observations
#' Y[sample(150*144, as.integer(0.2*150*144))] <- NA
#'
#' results <- fqpca(Y = Y[1:100,], npc = 2, quantile.value = 0.5)
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


#' Fit Yhat
#'
#' ## S3 method for class 'fqpca_object'
#' Given an fqpca_object model, estimates Yhat for different pve values.
#'
#' @param object An object output of the fqpca function.
#' @param pve Percentage of explained variability used in Yhat estimation.
#' @param ... further arguments passed to or from other methods.
#'
#' @return The normalized matrix of scores.
#' @export
#'
#' @examples
#' # Generate fake dataset with 150 observations and 144 time points
#'
#' Y <- matrix(rep(sin(seq(0, 2*pi, length.out = 144)), 150), byrow = TRUE, nrow = 150)
#' Y <- Y + matrix(rnorm(150*144, 0, 0.4), nrow = 150)
#'
#' # Add missing observations
#' Y[sample(150*144, as.integer(0.2*150*144))] <- NA
#'
#' results <- fqpca(Y = Y[1:100,], npc = 2, quantile.value = 0.5)
#'
#' Yhat <- fitted(object = results, pve=0.99)
#'
fitted.fqpca_object <- function(object, pve=0.95, ...)
{
  score.variance = cumsum(object$pve)
  n.components = min(which(score.variance > pve))
  Y.reconstruction = object$scores[, 1:(n.components+1)] %*% t(object$loadings[, 1:(n.components+1)])
  return(Y.reconstruction)
}

# BASIC PLOT ------------------------------------------------------------------

#' Plot FQPCA loading functions
#'
#' ## S3 method for class 'fqpca'
#' Given a fqpca object, plot the loading functions
#'
#' @param x An object output of the fqpca function.
#' @param pve Percentage of explained variability plotted.
#' @param ... further arguments passed to or from other methods.
#'
#' @return The plot of loadings
#' @export
#'
#' @examples
#' # Generate fake dataset with 150 observations and 144 time points
#'
#' Y <- matrix(rep(sin(seq(0, 2*pi, length.out = 144)), 150), byrow = TRUE, nrow = 150)
#' Y <- Y + matrix(rnorm(150*144, 0, 0.4), nrow = 150)
#'
#' # Add missing observations
#' Y[sample(150*144, as.integer(0.2*150*144))] <- NA
#'
#' results <- fqpca(Y = Y, npc = 2, quantile.value = 0.5)
#'
#' plot(results)
plot.fqpca_object <- function(x, pve=0.99, ...)
{
  if (!base::inherits(x, "fqpca_object"))
  {
    stop('The object must be of class fqpca_object')
  }

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
  graphics::title('Principal components functions')
}
