
# SCORES ----------------------------------------------------------------------

#' @title Compute scores between
#' @description Inner function to compute the between level scores of the mfqpca methodology.
#' @param Y The \eqn{(N \times T)} matrix of observed time instants.
#' @param Y.mask Mask matrix of the same dimensions as Y indicating which observations in Y are known.
#' @param group An N dimensional array indicating the hierarchical structure of the data. Elements in the array with the same value indicate they are repeated measures of the same individual.
#' @param intercept A T dimensional vector of intercept values.
#' @param loadings Matrix of loading coefficients.
#' @param quantile.value The quantile considered.
#' @return The matrix of between level scores.
mfqpca_compute_scores_between <- function(Y, Y.mask, group, intercept, loadings, quantile.value)
{
  unique.subjects <- unique(group)
  n.subjects <- length(unique.subjects)
  npc <- ncol(loadings)
  scores <- matrix(0, n.subjects, npc)

  # Pre-compute subject indices for efficiency
  subject_indices <- split(seq_along(group), group)

  for(idx.g in seq_len(n.subjects))
  {
    subj_id <- unique.subjects[idx.g]
    subj_rows <- subject_indices[[as.character(subj_id)]]

    Y.subject <- Y[subj_rows, , drop = FALSE]
    Y.subject.mask <- Y.mask[subj_rows, , drop = FALSE]

    # Pre-calculate total observations for pre-allocation
    total_obs <- sum(Y.subject.mask)

    if (total_obs > 0) {
      # Pre-allocate vectors to avoid growing
      Y.vector <- numeric(total_obs)
      loadings.obs <- matrix(0, total_obs, npc)
      intercept.obs <- numeric(total_obs)

      # Fill pre-allocated vectors
      obs_idx <- 1
      for (row_idx in seq_len(nrow(Y.subject))) {
        mask_row <- Y.subject.mask[row_idx, ]
        obs_cols <- which(mask_row)
        n_obs <- length(obs_cols)

        if (n_obs > 0) {
          end_idx <- obs_idx + n_obs - 1
          Y.vector[obs_idx:end_idx] <- Y.subject[row_idx, obs_cols]
          loadings.obs[obs_idx:end_idx, ] <- loadings[obs_cols, , drop = FALSE]
          intercept.obs[obs_idx:end_idx] <- intercept[obs_cols]
          obs_idx <- end_idx + 1
        }
      }

      Y.vector <- Y.vector - intercept.obs
      scores[idx.g, ] <- quantreg::rq.fit.br(y = Y.vector, x = loadings.obs, tau = quantile.value)$coefficients
    }
  }
  return(scores)
}

#' @title Compute scores within
#' @description Inner function to compute the within level scores of the mfqpca methodology.
#' @param Y The \eqn{(N \times T)} matrix of observed time instants.
#' @param Y.mask Mask matrix of the same dimensions as Y indicating which observations in Y are known.
#' @param loadings Matrix of loading coefficients.
#' @param quantile.value The quantile considered.
#' @return The matrix of within level scores.
mfqpca_compute_scores_within <- function(Y, Y.mask, loadings, quantile.value)
{
  # Initialize the matrix of scores
  n.obs <- base::nrow(Y)
  npc <- base::ncol(loadings)
  scores <-  base::matrix(0, n.obs, npc)
  for(i in base::seq(n.obs))
  {
    Yi <- Y[i, Y.mask[i, ]]
    loadings.obs <- loadings[Y.mask[i, ], ]
    scores[i,] <-  quantreg::rq.fit.br(y=Yi, x=loadings.obs, tau=quantile.value)$coefficients
  }
  return(scores)
}

# SPLINES ---------------------------------------------------------------------

#' @title Compute spline coefficients between
#' @description Inner function to compute the spline coefficients of the mfqpca methodology using quantreg
#' @param Y.vector vectorised version of Y, the \eqn{(N \times T)} matrix of observed time instants.
#' @param Y.mask Mask matrix of the same dimensions as Y indicating which observations in Y are known.
#' @param scores Initial matrix of scores.
#' @param spline.basis The spline basis matrix.
#' @param quantile.value The quantile considered.
#' @param method Method used in the resolution of the quantile regression model. It currently accepts the methods \code{c('conquer', 'quantreg')}.
#' @return The matrix of spline coefficients splines coefficients.
mfqpca_compute_splines_between <- function(Y.vector, Y.mask, scores, spline.basis, quantile.value, method)
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

#' @title Compute spline coefficients within
#' @description Inner function to compute the spline coefficients of the mfqpca methodology using conquer2
#' @param Y \eqn{(N \times T)} matrix of observed time instants.
#' @param Y.mask Mask matrix of the same dimensions as Y indicating which observations in Y are known.
#' @param scores Initial matrix of scores.
#' @param spline.basis The spline basis matrix.
#' @param quantile.value The quantile considered.
#' @param method Method used in the resolution of the quantile regression model. It currently accepts the methods \code{c('conquer', 'quantreg')}.
#' @return The matrix of spline coefficients splines coefficients.
mfqpca_compute_splines_within <- function(Y, Y.mask, scores, spline.basis, quantile.value, method)
{
  n.obs <- base::nrow(Y)
  npc <- base::ncol(scores)
  n.basis <- base::ncol(spline.basis)
  Y.list <- lapply(base::seq(n.obs), function(j) Y[j, Y.mask[j, ]])
  Y.vector <- unlist(Y.list, use.names = FALSE)
  tensor.matrix <- base::matrix(0, nrow=base::length(Y.vector), ncol=npc*n.basis)
  row.idx <- 1
  for(i in base::seq(n.obs))
  {
    scores.i <- scores[i, , drop=FALSE]
    tmp.splines <- spline.basis[Y.mask[i, ], , drop=FALSE]
    n.time.i <- base::nrow(tmp.splines)
    tensor.matrix[row.idx:(row.idx + n.time.i - 1), ] <- base::kronecker(scores.i, tmp.splines)
    row.idx <- row.idx+n.time.i
  }
  if(method == 'conquer'){
    B.vector <- conquer2::conquer(X=tensor.matrix, Y=Y.vector, tau=quantile.value)$coeff
  }else if(method == 'quantreg')
  {
    B.vector <- quantreg::rq.fit.fnb(y=Y.vector, x=tensor.matrix, tau=quantile.value)$coefficients
  }
  spline.coefficients <- base::matrix(B.vector, byrow=FALSE, ncol=npc)
  return(spline.coefficients)
}

# LOADINGS --------------------------------------------------------------------

#' @title Compute loadings (aka principal components)
#' @description Inner function to compute the loadings of the mfqpca methodology when using conquer or other intercept alternatives
#' @param spline.basis The spline basis matrix.
#' @param spline.coefficients the matrix of spline coefficients.
#' @return The matrix of loadings
mfqpca_compute_loadings_between <- function(spline.basis, spline.coefficients)
{
  intercept.spline.basis <-  spline.basis[, -1]
  intercept.part <- cbind(1, intercept.spline.basis) %*% matrix(spline.coefficients[,1], ncol=1)
  fqpc.part <- spline.basis %*% spline.coefficients[, -1, drop=FALSE]
  return(list(intercept=intercept.part, loadings=fqpc.part))
}

# ROTATION --------------------------------------------------------------------

#' @title Rotation of fqpca loadings and scores
#' @description Performs the rotation of loadings and scores in order to ensure the solution is unique.
#' @param scores Matrix of scores.
#' @param intercept population intercept
#' @param loadings Matrix of loadings.
#' @return The rotated matrices of loadings and scores and the rotation matrix.
mfqpca_rotate_scores_and_loadings_between <- function(scores, intercept, loadings)
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

#' @title Rotation of mfqpca within-level loadings and scores
#' @description Performs the rotation of loadings and scores in order to ensure the solution is unique.
#' @param loadings Matrix of loadings.
#' @param scores Matrix of scores.
#' @return The rotated matrices of loadings and scores and the rotation matrix.
mfqpca_rotate_scores_and_loadings_within <- function(scores, loadings)
{
  npc <- base::ncol(loadings)

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

  results <- list(scores=scores, loadings=loadings, rotation.matrix=rotation.matrix)
  return(results)
}

# OTHER FUNCTIONS -------------------------------------------------------------

#' @title Check input parameters
#' @description Check input parameters
#' @param npc.between The number of estimated between components.
#' @param npc.within The number of estimated within components.
#' @param quantile.value The quantile considered.
#' @param periodic Boolean indicating if the data is expected to be periodic (start coincides with end) or not.
#' @param splines.df Degrees of freedom for the splines.
#' @param splines.method Method used in the resolution of the splines quantile regression model. It currently accepts the methods \code{c('conquer', 'quantreg')}.
#' @param tol Tolerance on the convergence of the algorithm.
#' @param max.iters Maximum number of iterations.
#' @param verbose Boolean indicating the verbosity.
#' @param seed Seed for the random generator number.
#' @return No return
mfqpca_check_params <- function(npc.between, npc.within, quantile.value, periodic, splines.df, splines.method, tol, max.iters, verbose, seed)
{
  # Check 'npc.between': integer number, positive
  if (!is.numeric(npc.between) || length(npc.between) != 1 || npc.between %% 1 != 0 || npc.between <= 0) {
    stop("Invalid input for 'npc.between': ", npc.between,
         ". Expected a positive integer number.")
  }

  if (!is.numeric(npc.within) || length(npc.within) != 1 || npc.within %% 1 != 0 || npc.within <= 0) {
    stop("Invalid input for 'npc.within': ", npc.within,
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
      splines.df %% 1 != 0 || splines.df < npc.between) {
    stop("Invalid input for 'splines.df': ", splines.df,
         ". Expected an integer number larger or equal than 'npc.between' (", npc.between, ").")
  }

  # Check 'method': must be either 'conquer' or 'quantreg'
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

#' @title Check input data
#' @description Check input data
#' @param data An \eqn{(N \times T)} matrix, a tf object from the tidyfun package or a data.frame containing the functional data as a tf column.
#' @param colname The name of the column containing the functional data. Use only if data is a dataframe and colname is a column in the dataframe.
#' @param group The grouping structure of the hierarchical data. If data is a dataframe this can be a name pointing to a column in the dataframe. If data is any other object it must be an array with length equal to the number of rows of data.
#' @return A list containing two elements:
#'         1. An unnamed matrix of the functional data.
#'         2. The processed group vector.
mfqpca_check_input_data <- function(data, group, colname)
{

  # Helper function to check that the length of group matches a given number of rows
  check_group_length <- function(grp, n) {
    if (length(grp) != n) {
      stop("Length of 'group' (", length(grp),
           ") does not match the number of rows (", n, ").")}}

  # Case 1: data is a matrix
  if (is.matrix(data)) {
    # For matrix, a single string is not allowed for group.
    if (is.character(group) && length(group) == 1) {
      stop("For matrix input, 'group' cannot be a single string; it must be an array matching the number of rows.")
    }
    check_group_length(group, nrow(data))
    return(list(unname(data), group))
  }

  # Case 2: data is a 'tf' object
  if (tf::is_tf(data)) {
    # For tf objects, group must be an array matching the rows, not a single string.
    if (is.character(group) && length(group) == 1) {
      stop("For 'tf' input, 'group' cannot be a single string; it must be an array matching the number of rows.")
    }
    mdata <- unname(as.matrix(data))
    check_group_length(group, nrow(mdata))
    return(list(mdata, group))
  }

  # Case 3: data is a data.frame
  if (is.data.frame(data)) {
    # Verify that colname is provided as a single string and it exists in data.
    if (missing(colname) || !is.character(colname) || length(colname) != 1) {
      stop("Data is a data.frame, but 'colname' is not provided as a single string corresponding to a column name.")
    }
    if (!(colname %in% colnames(data))) {
      stop("Data is a data.frame, but 'colname' ('", colname, "') does not match any column in the data.")
    }

    # Extract the functional data column and check its type.
    func_data_column <- data[[colname]]
    if (!tf::is_tf(func_data_column)) {
      stop("The column '", colname, "' in the data.frame is not a 'tf' object.")
    }

    processed_func_data <- unname(as.matrix(func_data_column))

    # Process the 'group' input for data.frame:
    actual_group_vector <- NULL # Initialize

    if (is.character(group) && length(group) == 1) {
      # Option 1: If group is a string, it must correspond to a column name in data.
      if (group %in% colnames(data)) {
        actual_group_vector <- data[[group]]
      } else {
        stop("For data.frame input, when 'group' is provided as a string ('", group, "'), it must match a column name in the data.")
      }
    } else {
      # Option 2: group is treated as an array/vector.
      actual_group_vector <- group
    }

    # Check length of the processed group vector against the number of rows in the data.frame
    check_group_length(actual_group_vector, nrow(data))

    return(list(processed_func_data, actual_group_vector))
  }

  # If none of the supported types, throw an error.
  stop("Data is not a matrix, 'tf' object, or data.frame. Cannot process the input data.")
}

# MAIN ------------------------------------------------------------------------

#' @title MFQPCA (Multilevel Functional Quantile Principal Component Analysis)
#' @description Solves the multilevel functional quantile principal component analysis methodology
#' @param data An \eqn{(N \times T)} matrix, a tf object from the tidyfun package or a data.frame containing the functional data as a tf column.
#' @param group Either a string or an array. If it is a string, it must point to the grouping variable in data only if data is a dataframe. If an array, it must be the N dimensional array indicating the hierarchical structure of the data. Elements in the array with the same value indicate they are repeated measures of the same individual.
#' @param colname The name of the column containing the functional data. Use only if data is a dataframe and colname is a column in the dataframe.
#' @param npc.between The number of estimated between level components.
#' @param npc.within The number of estimated within level components.
#' @param quantile.value The quantile considered.
#' @param periodic Boolean indicating if the data is expected to be periodic (start coincides with end) or not.
#' @param splines.df Degrees of freedom for the splines.
#' @param splines.method Method used in the resolution of the splines quantile regression model. It currently accepts the methods \code{c('conquer', 'quantreg')}.
#' @param tol Tolerance on the convergence of the algorithm.
#' @param max.iters Maximum number of iterations.
#' @param verbose Boolean indicating the verbosity.
#' @param seed Seed for the random generator number.
#' @return mfqpca_object
#' @export
#' @examples
#'
#' n.individuals <- 20
#' n.repeated <- 10
#' n.time = 144
#' N <- n.repeated * n.individuals
#'
#' group <- rep(1:n.individuals, each=n.repeated)
#'
#' # Define score values using a normal distribution
#' c1.vals <- rnorm(n.individuals)
#' c1.vals <- c1.vals[match(group, unique(group))]
#' c2.vals <- rnorm(N)
#'
#' # Define principal components
#' pcb <- sin(seq(0, 2*pi, length.out = n.time))
#' pcw <- cos(seq(0, 2*pi, length.out = n.time))
#'
#' # Generate a data matrix and add missing observations
#' Y <- c1.vals * matrix(pcb, nrow = N, ncol=n.time, byrow = TRUE) +
#' c2.vals * matrix(pcw, nrow = N, ncol=n.time, byrow = TRUE)
#'
#' # Add noise
#' Y <- Y + matrix(rnorm(N*n.time, 0, 0.4), nrow = N)
#' Y[sample(N*n.time, as.integer(0.2*N))] <- NA
#'
#' results <- mfqpca(data = Y, group = group, npc.between = 1, npc.within=1, quantile.value = 0.5)
mfqpca <- function(
    data, group, colname=NULL, npc.between = 2,  npc.within = 2, quantile.value = 0.5,
    periodic = TRUE, splines.df = 10, splines.method = 'conquer',
    tol = 1e-2, max.iters = 10, verbose = FALSE, seed = NULL)
{

  global.start.time <- base::Sys.time()

  # Get the input parameters
  formal_names <- base::setdiff(names(formals(sys.function())), "...")
  inputs <- base::mget(formal_names, envir = environment())
  inputs$data = NULL
  inputs$group = NULL
  inputs$colname = NULL

  # Check the input parameters except Y
  mfqpca_check_params(npc.between=npc.between, npc.within=npc.within, quantile.value=quantile.value, periodic=periodic,
                     splines.df=splines.df, splines.method=splines.method, tol=tol, max.iters=max.iters,
                     verbose=verbose, seed=seed)

  # Check Y and colname and return an unnamed matrix
  r <- mfqpca_check_input_data(data=data, colname=colname, group=group)
  Y <- r[[1]]
  group <- r[[2]]

  # If seed is provided, set seed for computations
  if(!base::is.null(seed)){base::set.seed(seed)}

  # Step 1: Initialize values
  n.obs <- base::nrow(Y)
  n.time <- base::dim(Y)[2]
  Y.axis <- base::seq(0, 1, length.out = n.time)
  Y.mask <- !base::is.na(Y)
  Y.list <- lapply(base::seq(n.obs), function(j) Y[j, Y.mask[j, ]])
  Y.vector <- unlist(Y.list, use.names = FALSE)

  # Initialize splines knots
  if(periodic)
  {
    knots <- base::seq(0, 1, length.out = 2 + splines.df - 1)  # 2 boundary knots + df - intercept
  } else{
    knots <- base::seq(0, 1, length.out = 2 + splines.df - 3 - 1) # 2 boundary knots + df - degree - intercept
  }
  knots <- knots[2:(base::length(knots)-1)]

  # Initialize spline basis
  spline.basis <- pbs::pbs(Y.axis, degree = 3, knots = knots, intercept = TRUE, periodic=periodic, Boundary.knots = c(min(Y.axis), max(Y.axis)))
  n.basis = base::ncol(spline.basis)

  loop.start.time <- base::Sys.time()

  partial.execution.time <- list(
    scores.between = NULL,
    scores.within = NULL,
    loadings.between = NULL,
    loadings.within = NULL)

  # Between level
  spline.coefficients.between <- base::matrix(stats::rnorm(n.basis * (npc.between+1), mean = 0, sd = 1), nrow=n.basis, ncol=npc.between+1)
  loadings_list <- spline.basis %*% spline.coefficients.between
  intercept <- loadings_list[, 1]
  loadings.between <- loadings_list[, -1, drop=FALSE]
  scores.between <- mfqpca_compute_scores_between(Y=Y, Y.mask=Y.mask, group=group, intercept=intercept, loadings=loadings.between, quantile.value=quantile.value)
  scores.between.full <- scores.between[match(group, unique(group)), , drop=FALSE]

  # Residuals
  Y.between.fitted <- sweep(scores.between.full %*% t(loadings.between), MARGIN = 2, STATS = intercept, FUN = "+")
  Y.within <- Y - Y.between.fitted

  # Within level
  spline.coefficients.within <- base::matrix(stats::rnorm(n.basis * (npc.within), mean = 0, sd = 1), nrow=n.basis, ncol=npc.within)
  loadings.within <- spline.basis %*% spline.coefficients.within
  scores.within <- mfqpca_compute_scores_within(Y=Y.within, Y.mask=Y.mask, loadings=loadings.within, quantile.value=quantile.value)

  # Rotation
  rotation.result <- mfqpca_rotate_scores_and_loadings_between(scores=scores.between, intercept=intercept, loadings=loadings.between)
  intercept = rotation.result$intercept
  loadings.between <- rotation.result$loadings
  scores.between <- rotation.result$scores
  scores.between.full <- scores.between[match(group, unique(group)), , drop=FALSE]
  rotation.matrix.between <- rotation.result$rotation.matrix

  rotation.result <- mfqpca_rotate_scores_and_loadings_within(scores=scores.within, loadings=loadings.within)
  loadings.within <- rotation.result$loadings
  scores.within <- rotation.result$scores
  rotation.matrix.within <- rotation.result$rotation.matrix

  # Compute objective value function
  objective.function.array <- compute_objective_value(quantile.value=quantile.value, Y=Y, scores=cbind(scores.between.full, scores.within), intercept=intercept, loadings=cbind(loadings.between, loadings.within))

  loop.execution.time <- difftime(base::Sys.time(), loop.start.time, units = 'secs')

  if(verbose)
  {
    message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '. Iteration: ', 1, ' completed in ', base::round(loop.execution.time, 3), ' seconds',
            '\n', 'Objective function   i   value: ', round(objective.function.array[1], 4),
            '\n', '____________________________________________________________')
  }

  convergence = F
  for(i in 2:max.iters)
  {
    loop.start.time <- base::Sys.time()

    # Between level
    tmp <- Sys.time()
    spline.coefficients.between <- mfqpca_compute_splines_between(Y.vector=Y.vector, Y.mask=Y.mask, scores=scores.between.full, spline.basis=spline.basis, quantile.value=quantile.value, method=splines.method)
    partial.execution.time$loadings.between <- c(partial.execution.time$loadings.between, as.numeric(difftime(base::Sys.time(), tmp, units = 'secs')))

    loadings_list <- mfqpca_compute_loadings_between(spline.basis = spline.basis, spline.coefficients = spline.coefficients.between)
    intercept <- loadings_list$intercept
    loadings.between <- loadings_list$loadings

    tmp <- Sys.time()
    scores.between <- mfqpca_compute_scores_between(Y=Y, Y.mask=Y.mask, group=group, intercept=intercept, loadings=loadings.between, quantile.value=quantile.value)
    partial.execution.time$scores.between <- c(partial.execution.time$scores.between, as.numeric(difftime(base::Sys.time(), tmp, units = 'secs')))

    scores.between.full <- scores.between[match(group, unique(group)), , drop=FALSE]

    # Compute residuals
    Y.between.fitted <- sweep(scores.between.full %*% t(loadings.between), MARGIN = 2, STATS = intercept, FUN = "+")
    Y.within <- Y - Y.between.fitted

    # Within level
    tmp <- Sys.time()
    spline.coefficients.within <- mfqpca_compute_splines_within(Y=Y.within, Y.mask=Y.mask, scores=scores.within, spline.basis=spline.basis, quantile.value=quantile.value, method=splines.method)
    partial.execution.time$loadings.within <- c(partial.execution.time$loadings.within, as.numeric(difftime(base::Sys.time(), tmp, units = 'secs')))

    loadings.within <- spline.basis %*% spline.coefficients.within

    tmp <- Sys.time()
    scores.within <- mfqpca_compute_scores_within(Y=Y.within, Y.mask=Y.mask, loadings=loadings.within, quantile.value=quantile.value)
    partial.execution.time$scores.within <- c(partial.execution.time$scores.within, as.numeric(difftime(base::Sys.time(), tmp, units = 'secs')))

    # Rotation
    rotation.result <- mfqpca_rotate_scores_and_loadings_between(scores=scores.between, intercept=intercept, loadings=loadings.between)
    intercept = rotation.result$intercept
    loadings.between <- rotation.result$loadings
    scores.between <- rotation.result$scores
    scores.between.full <- scores.between[match(group, unique(group)), , drop=FALSE]
    rotation.matrix.between <- rotation.result$rotation.matrix

    rotation.result <- mfqpca_rotate_scores_and_loadings_within(scores=scores.within, loadings=loadings.within)
    loadings.within <- rotation.result$loadings
    scores.within <- rotation.result$scores
    rotation.matrix.within <- rotation.result$rotation.matrix

    # Compute objective value function
    objective.function <- compute_objective_value(quantile.value=quantile.value, Y=Y, scores=cbind(scores.between.full, scores.within), intercept=intercept, loadings=cbind(loadings.between, loadings.within))
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

    # Avoid computational issues and perform early stopping if reaching diverging loops.
    if(i>=3 & (objective.function.array[i] > 10 * objective.function.array[2]))
    {
      message('Iteration ', i, '. Breaking diverging loop')
      break
    }
  }

  if(verbose & convergence){message("\u2705 Algorithm converged succesfully")}

  if(i == max.iters & !convergence){warning('\u274C Algorithm reached maximum number of iterations without convergence. Consider increasing the value of max.iters parameter')}

  # Compute explained variability
  pve.between <- compute_explained_variability(scores.between)
  pve.within <- compute_explained_variability(scores.within)

  global.execution.time <- difftime(base::Sys.time(), global.start.time, units = 'secs')

  results <- structure(list(
    intercept = intercept,
    loadings.between = loadings.between,
    scores.between = scores.between,
    scores.between.full = scores.between.full,
    pve.between = pve.between,
    loadings.within = loadings.within,
    scores.within = scores.within,
    pve.within = pve.within,
    group = group,
    objective.function.value = objective.function,
    objective.function.array = objective.function.array,
    execution.time = global.execution.time,
    rotation.matrix.between = rotation.matrix.between,
    rotation.matrix.within = rotation.matrix.within,
    spline.coefficients.between = spline.coefficients.between,
    spline.coefficients.within = spline.coefficients.within,
    spline.basis = spline.basis,
    n.iters = i,
    partial.execution.time = partial.execution.time,
    inputs = inputs),
    class = "mfqpca_object")

  return(results)
}

# ADDITIONAL FUNCTIONS --------------------------------------------------------

#' @title Helper function to compute fitted values
#' @description Helper function to compute fitted values
#' @param scores matrix of scores
#' @param loadings matrix of loadings
#' @param n.comp number of components used. If 0 then a matrix of 0s is returned.
mfqpca_reconstruct <- function(scores, loadings, n.comp)
{
  if (n.comp > 0)
  {
    r = scores[, seq_len(n.comp), drop = FALSE] %*% t(loadings[, seq_len(n.comp), drop = FALSE])
  }else{
    # Create a zero matrix
    r = matrix(0, nrow = nrow(scores), ncol = nrow(loadings))
  }
  return(r)
}

#' @title Predict mfqpca scores
#' @description S3 method for class 'fqpca_object' Given a new matrix Y, predicts the value of the scores associated to the given matrix.
#' @param object An object output of the fqpca function.
#' @param newdata The N by T matrix of observed time instants to be tested
#' @param newdata.group An N dimensional array indicating the hierarchical structure of the data. Elements in the array with the same value indicate they are repeated measures of the same individual.
#' @param ... further arguments passed to or from other methods.
#' @return The normalized matrix of scores.
#' @export
#' @examples
#'
#' n.individuals <- 20
#' n.repeated <- 10
#' n.time = 144
#' N <- n.repeated * n.individuals
#'
#' group <- rep(1:n.individuals, each=n.repeated)
#'
#' # Define score values using a normal distribution
#' c1.vals <- rnorm(n.individuals)
#' c1.vals <- c1.vals[match(group, unique(group))]
#' c2.vals <- rnorm(N)
#'
#' # Define principal components
#' pcb <- sin(seq(0, 2*pi, length.out = n.time))
#' pcw <- cos(seq(0, 2*pi, length.out = n.time))
#'
#' # Generate a data matrix and add missing observations
#' Y <- c1.vals * matrix(pcb, nrow = N, ncol=n.time, byrow = TRUE) +
#' c2.vals * matrix(pcw, nrow = N, ncol=n.time, byrow = TRUE)
#' Y <- Y + matrix(rnorm(N*n.time, 0, 0.4), nrow = N)
#' Y[sample(N*n.time, as.integer(0.2*N))] <- NA
#'
#' results <- mfqpca(data = Y, group=group, npc.between = 1, npc.within=1, quantile.value = 0.5)
#' predictions <- predict(object = results, newdata = Y[101:150,], newdata.group = group[101:150])
#'
predict.mfqpca_object <- function(object, newdata, newdata.group, ...)
{
  if(!base::inherits(object, "mfqpca_object")){stop('The object must be of class mfqpca_object')}

  r <- mfqpca_check_input_data(data=newdata, colname=object$inputs$colname, group=newdata.group)
  Y <- r[[1]]
  group <- r[[2]]

  Y.mask <- !base::is.na(Y)
  n.obs <- base::nrow(Y)
  n.time <- base::ncol(Y)
  if(sum(!Y.mask) == n.obs * n.time){stop('newdata contains no information. Check that there are values different from NA')}

  scores.between <- mfqpca_compute_scores_between(Y=Y, Y.mask=Y.mask, group=group, intercept=object$intercept, loadings=object$loadings.between, quantile.value=object$inputs$quantile.value)
  scores.between.full <- scores.between[match(group, unique(group)), , drop=FALSE]

  Y.between.fitted <- sweep(scores.between.full %*% t(object$loadings.between), MARGIN = 2, STATS = object$intercept, FUN = "+")
  Y.within <- Y - Y.between.fitted

  scores.within <- mfqpca_compute_scores_within(Y=Y.within, Y.mask=Y.mask, loadings=object$loadings.within, quantile.value=object$inputs$quantile.value)

  return(list(scores.between=scores.between, scores.between.full=scores.between.full, scores.within=scores.within))
}

#' @title Fit Yhat
#' @description S3 method for class 'mfqpca_object'. Given an mfqpca_object model, estimates Yhat for different pve values.
#' @param object An object output of the fqpca function.
#' @param pve.between If smaller than 1, taken as percentage of explained variability used in Yhat estimation. If greater than 1, taken as number of components  used in Yhat estimation.
#' @param pve.within If smaller than 1, taken as percentage of explained variability used in Yhat estimation. If greater than 1, taken as number of components  used in Yhat estimation.
#' @param ... further arguments passed to or from other methods.
#' @return The normalized matrix of scores.
#' @export
#' @examples
#'
#' n.individuals <- 20
#' n.repeated <- 10
#' n.time = 144
#' N <- n.repeated * n.individuals
#'
#' group <- rep(1:n.individuals, each=n.repeated)
#'
#' # Define score values using a normal distribution
#' c1.vals <- rnorm(n.individuals)
#' c1.vals <- c1.vals[match(group, unique(group))]
#' c2.vals <- rnorm(N)
#'
#' # Define principal components
#' pcb <- sin(seq(0, 2*pi, length.out = n.time))
#' pcw <- cos(seq(0, 2*pi, length.out = n.time))
#'
#' # Generate a data matrix and add missing observations
#' Y <- c1.vals * matrix(pcb, nrow = N, ncol=n.time, byrow = TRUE) +
#' c2.vals * matrix(pcw, nrow = N, ncol=n.time, byrow = TRUE)
#' Y <- Y + matrix(rnorm(N*n.time, 0, 0.4), nrow = N)
#' Y[sample(N*n.time, as.integer(0.2*N))] <- NA
#'
#' results <- mfqpca(data = Y, group=group, npc.between = 1, npc.within=1, quantile.value = 0.5)
#' Yhat <- fitted(object = results, pve.between=0.95, pve.within=0)
#'
fitted.mfqpca_object <- function(object, pve.between=0.95, pve.within=0.95, ...)
{
  if (!base::inherits(object, "mfqpca_object")){stop('The object must be of class mfqpca_object')}

  n.components.between <- obtain_npc(scores=object$scores.between, pve=pve.between)
  n.components.within <- obtain_npc(scores=object$scores.within, pve=pve.within)

  Y.fitted =
    mfqpca_reconstruct(object$scores.between.full, object$loadings.between, n.components.between) +
    mfqpca_reconstruct(object$scores.within, object$loadings.within, n.components.within)

  Y.fitted = sweep(Y.fitted, MARGIN = 2, STATS = object$intercept, FUN = "+")

  return(Y.fitted)
}

#' @title Plot fqpca loading functions with ggplot2
#' @description S3 method for class 'fqpca_object'. Given a fqpca object, plot
#'              the loading functions using ggplot2 facets.
#' @param x An object output of the fqpca function.
#' @param pve.between Percentage of explained variability to determine the number of components.
#' @param pve.within Percentage of explained variability to determine the number of components.
#' @param ... Further arguments passed to or from other methods.
#' @return A ggplot object plotting the intercept and FQPC curves.
#' @importFrom rlang .data
#' @export
#' @examples
#' n.individuals <- 20
#' n.repeated <- 10
#' n.time = 144
#' N <- n.repeated * n.individuals
#'
#' group <- rep(1:n.individuals, each=n.repeated)
#'
#' # Define score values using a normal distribution
#' c1.vals <- rnorm(n.individuals)
#' c1.vals <- c1.vals[match(group, unique(group))]
#' c2.vals <- rnorm(N)
#'
#' # Define principal components
#' pcb <- sin(seq(0, 2*pi, length.out = n.time))
#' pcw <- cos(seq(0, 2*pi, length.out = n.time))
#'
#' # Generate a data matrix and add missing observations
#' Y <- c1.vals * matrix(pcb, nrow = N, ncol=n.time, byrow = TRUE) +
#' c2.vals * matrix(pcw, nrow = N, ncol=n.time, byrow = TRUE) +
#' matrix(seq(1, 10, length.out=n.time), nrow = N, ncol=n.time, byrow = TRUE)
#' Y <- Y + matrix(rnorm(N*n.time, 0, 0.4), nrow = N)
#' Y[sample(N*n.time, as.integer(0.2*N))] <- NA
#'
#' results <- mfqpca(
#'    data = Y, group=group, npc.between = 1, npc.within=1, quantile.value = 0.5,
#'    periodic=FALSE, seed=1)
#' plot(results)
plot.mfqpca_object <- function(x, pve.between=0.95, pve.within=0.95, ...)
{
  if (!inherits(x, "mfqpca_object")) {stop("The x must be of class mfqpca_object")}

  n.components.between <- obtain_npc(scores=x$scores.between, pve=pve.between)
  n.components.within <- obtain_npc(scores=x$scores.within, pve=pve.within)

  intercept <- x$intercept

  # Create a data frame for the intercept curve
  intercept_df <- data.frame(
    Time = seq_len(length(intercept)),
    Loading = intercept,
    Component = "Intercept")

  # Extract Between components
  loadings.between <- x$loadings.between
  loadings.between.filter <- loadings.between[, 1:n.components.between, drop = FALSE]

  # Reshape loadings-between into a long format data frame
  fqpc_between_df <- data.frame(
    Time = rep(seq_len(nrow(loadings.between.filter)), times = ncol(loadings.between.filter)),
    Loading = as.vector(loadings.between.filter),
    Component = rep(paste0("MFQPC Between ", seq_len(ncol(loadings.between.filter))), each = nrow(loadings.between.filter)))

  loadings.within <- x$loadings.within

  # Extract Within components
  loadings.within.filter <- loadings.within[, 1:n.components.within, drop = FALSE]

  fqpc_within_df <- data.frame(
    Time = rep(seq_len(nrow(loadings.within.filter)), times = ncol(loadings.within.filter)),
    Loading = as.vector(loadings.within.filter),
    Component = rep(paste0("MFQPC Within ", seq_len(ncol(loadings.within.filter))), each = nrow(loadings.within.filter)))

  # Combine the intercept and FQPC data frames
  plot_data <- rbind(intercept_df, fqpc_between_df, fqpc_within_df)
  plot_data$Component = factor(plot_data$Component, levels = c('Intercept',
                                                               paste0('MFQPC Between ', 1:ncol(loadings.between.filter)),
                                                               paste0('MFQPC Within ', 1:ncol(loadings.within.filter))))

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
