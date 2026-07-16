
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
mfqpca_compute_scores_between <- function(
    Y,
    Y.mask,
    group,
    intercept,
    loadings,
    quantile.value)
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
      if(length(Y.vector) < ncol(loadings.obs)){stop('All subject must have a larger number of timepoints than npcs.')}
      scores[idx.g, ] <- quantreg::rq.fit.br(y = Y.vector, x = loadings.obs, tau = quantile.value)$coefficients
    }
  }
  return(scores)
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
    data,
    group,
    colname=NULL,
    npc.between = 2,
    npc.within = 2,
    quantile.value = 0.5,
    periodic = TRUE,
    splines.df = 10,
    splines.method = 'conquer',
    tol = 1e-2,
    max.iters = 10,
    verbose = FALSE,
    seed = NULL)
{

  global.start.time <- Sys.time()

  # Get the input parameters
  formal_names <- setdiff(names(formals(sys.function())), "...")
  inputs <- mget(formal_names, envir = environment())
  inputs$data = NULL
  inputs$group = NULL

  # Check the input parameters except Y
  mfqpca_check_params(
    npc.between=npc.between,
    npc.within=npc.within,
    quantile.value=quantile.value,
    periodic=periodic,
    splines.df=splines.df,
    splines.method=splines.method,
    tol=tol,
    max.iters=max.iters,
    verbose=verbose,
    seed=seed)

  # Check Y and colname and return an unnamed matrix
  if(is.null(group)){stop("'group' must be provided for the mfqpca methodology.")}
  r <- coerce_functional_input(data=data, colname=colname, group=group)
  Y <- r$Y
  group <- r$group

  # If seed is provided, set seed for computations
  if(!is.null(seed)){set.seed(seed)}

  # Step 1: Initialize values
  n.obs <- nrow(Y)
  n.time <- dim(Y)[2]
  Y.axis <- seq(0, 1, length.out = n.time)
  Y.mask <- !is.na(Y)
  Y.list <- lapply(seq_len(n.obs), function(j) Y[j, Y.mask[j, ]])
  Y.vector <- unlist(Y.list, use.names = FALSE)
  function.warnings <- list(
    splines.between = FALSE,
    scores.between = FALSE,
    splines.within = FALSE,
    scores.within = FALSE,
    rotation = FALSE,
    diverging.loop = FALSE)

  spline.basis <- pbs::pbs(
    Y.axis,
    degree = 3,
    df = splines.df,
    intercept = TRUE,
    periodic=periodic,
    Boundary.knots = c(min(Y.axis), max(Y.axis)))
  n.basis <- ncol(spline.basis)

  i <- 1
  loop.start.time <- Sys.time()

  partial.execution.time <- list(
    scores.between = NULL,
    scores.within = NULL,
    splines.between = NULL,
    splines.within = NULL)

  # Between level
  spline.coefficients.between <- matrix(stats::rnorm(n.basis * (npc.between+1), mean = 0, sd = 1), nrow=n.basis, ncol=npc.between+1)
  loadings_list <- spline.basis %*% spline.coefficients.between
  intercept <- loadings_list[, 1]
  loadings.between <- loadings_list[, -1, drop=FALSE]
  scores.between <- mfqpca_compute_scores_between(
    Y=Y,
    Y.mask=Y.mask,
    group=group,
    intercept=intercept,
    loadings=loadings.between,
    quantile.value=quantile.value)
  scores.between.full <- scores.between[match(group, unique(group)), , drop=FALSE]

  # Residuals
  Y.between.fitted <- sweep(
    scores.between.full %*% t(loadings.between),
    MARGIN = 2,
    STATS = intercept,
    FUN = "+")
  Y.within <- Y - Y.between.fitted

  # Within level
  spline.coefficients.within <- matrix(stats::rnorm(n.basis * (npc.within), mean = 0, sd = 1), nrow=n.basis, ncol=npc.within)
  loadings.within <- spline.basis %*% spline.coefficients.within
  scores.within <- compute_scores(
    Y=Y.within,
    Y.mask=Y.mask,
    loadings=loadings.within,
    quantile.value=quantile.value)

  # Rotation
  rotation.result <- try(rotate_factors(scores=scores.between, loadings=loadings.between, intercept=intercept), silent = FALSE)
  if(inherits(rotation.result, "try-error"))
  {
    warning('Iteration 1. Failed between level rotation process. Skipping rotation on this iteration.')
    function.warnings$rotation <- TRUE
    rotation.result <- list(
      scores=scores.between,
      intercept=intercept,
      loadings=loadings.between,
      rotation.matrix=diag(npc.between))
  }
  intercept = rotation.result$intercept
  loadings.between <- rotation.result$loadings
  scores.between <- rotation.result$scores
  scores.between.full <- scores.between[match(group, unique(group)), , drop=FALSE]
  rotation.matrix.between <- rotation.result$rotation.matrix

  rotation.result <- try(rotate_factors(scores=scores.within, loadings=loadings.within), silent = FALSE)
  if(inherits(rotation.result, "try-error"))
  {
    warning('Iteration 1. Failed within level rotation process. Skipping rotation on this iteration.')
    function.warnings$rotation <- TRUE
    rotation.result <- list(scores=scores.within, loadings=loadings.within, rotation.matrix=diag(npc.within))
  }
  loadings.within <- rotation.result$loadings
  scores.within <- rotation.result$scores
  rotation.matrix.within <- rotation.result$rotation.matrix

  # Compute objective value function
  objective.function <- compute_objective_value(
    quantile.value=quantile.value,
    Y=Y,
    scores=cbind(scores.between.full, scores.within),
    intercept=intercept,
    loadings=cbind(loadings.between, loadings.within))
  objective.function.array <- objective.function

  loop.execution.time <- difftime(Sys.time(), loop.start.time, units = 'secs')

  if(verbose)
  {
    message(
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      '. Iteration: ', 1, ' completed in ', round(loop.execution.time, 3), ' seconds', '\n',
      'Objective function   i   value: ', round(objective.function.array[1], 4), '\n',
      '____________________________________________________________')
  }

  convergence = FALSE
  for(i in 2:max.iters)
  {
    loop.start.time <- Sys.time()

    # Snapshot the previous iteration state so it can be fully restored if any step fails
    spline.coefficients.between.prev <- spline.coefficients.between
    intercept.prev <- intercept
    loadings.between.prev <- loadings.between
    scores.between.prev <- scores.between
    scores.between.full.prev <- scores.between.full
    spline.coefficients.within.prev <- spline.coefficients.within
    loadings.within.prev <- loadings.within
    scores.within.prev <- scores.within

    # Between level
    tmp <- Sys.time()
    spline.coefficients.between <- try(fit_tensor_quantile_regression(
      Y.vector=Y.vector,
      Y.mask=Y.mask,
      scores=scores.between.full,
      spline.basis=spline.basis,
      quantile.value=quantile.value,
      method=splines.method), silent = FALSE)
    partial.execution.time$splines.between <- c(partial.execution.time$splines.between, as.numeric(difftime(Sys.time(), tmp, units = 'secs')))
    if(!is.matrix(spline.coefficients.between))
    {
      warning('Iteration ', i, '. Failed computation of between level spline coefficients. Providing results from previous iteration.')
      function.warnings$splines.between <- TRUE
      spline.coefficients.between <- spline.coefficients.between.prev
      break
    }

    loadings_list <- compute_loadings(spline.basis = spline.basis, spline.coefficients = spline.coefficients.between)
    intercept <- loadings_list$intercept
    loadings.between <- loadings_list$loadings

    tmp <- Sys.time()
    scores.between <- try(mfqpca_compute_scores_between(
      Y=Y,
      Y.mask=Y.mask,
      group=group,
      intercept=intercept,
      loadings=loadings.between,
      quantile.value=quantile.value), silent = FALSE)
    partial.execution.time$scores.between <- c(partial.execution.time$scores.between, as.numeric(difftime(Sys.time(), tmp, units = 'secs')))
    if(!is.matrix(scores.between))
    {
      warning('Iteration ', i, '. Failed computation of between level scores. Providing results from previous iteration.')
      function.warnings$scores.between <- TRUE
      spline.coefficients.between <- spline.coefficients.between.prev
      intercept <- intercept.prev
      loadings.between <- loadings.between.prev
      scores.between <- scores.between.prev
      break
    }

    scores.between.full <- scores.between[match(group, unique(group)), , drop=FALSE]

    # Compute residuals
    Y.between.fitted <- sweep(
      scores.between.full %*% t(loadings.between),
      MARGIN = 2,
      STATS = intercept,
      FUN = "+")
    Y.within <- Y - Y.between.fitted

    # Within level
    tmp <- Sys.time()
    Y.within.vector <- unlist(lapply(seq_len(n.obs), function(j) Y.within[j, Y.mask[j, ]]), use.names = FALSE)
    spline.coefficients.within <- try(fit_tensor_quantile_regression(
      Y.vector=Y.within.vector,
      Y.mask=Y.mask,
      scores=scores.within,
      spline.basis=spline.basis,
      quantile.value=quantile.value,
      method=splines.method,
      intercept=FALSE), silent = FALSE)
    partial.execution.time$splines.within <- c(partial.execution.time$splines.within, as.numeric(difftime(Sys.time(), tmp, units = 'secs')))
    if(!is.matrix(spline.coefficients.within))
    {
      warning('Iteration ', i, '. Failed computation of within level spline coefficients. Providing results from previous iteration.')
      function.warnings$splines.within <- TRUE
      spline.coefficients.between <- spline.coefficients.between.prev
      intercept <- intercept.prev
      loadings.between <- loadings.between.prev
      scores.between <- scores.between.prev
      scores.between.full <- scores.between.full.prev
      spline.coefficients.within <- spline.coefficients.within.prev
      break
    }

    loadings.within <- spline.basis %*% spline.coefficients.within

    tmp <- Sys.time()
    scores.within <- try(compute_scores(
      Y=Y.within,
      Y.mask=Y.mask,
      loadings=loadings.within,
      quantile.value=quantile.value), silent = FALSE)
    partial.execution.time$scores.within <- c(partial.execution.time$scores.within, as.numeric(difftime(Sys.time(), tmp, units = 'secs')))
    if(!is.matrix(scores.within))
    {
      warning('Iteration ', i, '. Failed computation of within level scores. Providing results from previous iteration.')
      function.warnings$scores.within <- TRUE
      spline.coefficients.between <- spline.coefficients.between.prev
      intercept <- intercept.prev
      loadings.between <- loadings.between.prev
      scores.between <- scores.between.prev
      scores.between.full <- scores.between.full.prev
      spline.coefficients.within <- spline.coefficients.within.prev
      loadings.within <- loadings.within.prev
      scores.within <- scores.within.prev
      break
    }

    # Rotation
    rotation.result <- try(rotate_factors(scores=scores.between, loadings=loadings.between, intercept=intercept), silent = FALSE)
    if(inherits(rotation.result, "try-error"))
    {
      warning('Iteration ', i, '. Failed between level rotation process. Skipping rotation on this iteration.')
      function.warnings$rotation <- TRUE
      rotation.result <- list(
        scores=scores.between,
        intercept=intercept,
        loadings=loadings.between,
        rotation.matrix=diag(npc.between))
    }
    intercept = rotation.result$intercept
    loadings.between <- rotation.result$loadings
    scores.between <- rotation.result$scores
    scores.between.full <- scores.between[match(group, unique(group)), , drop=FALSE]
    rotation.matrix.between <- rotation.result$rotation.matrix

    rotation.result <- try(rotate_factors(scores=scores.within, loadings=loadings.within), silent = FALSE)
    if(inherits(rotation.result, "try-error"))
    {
      warning('Iteration ', i, '. Failed within level rotation process. Skipping rotation on this iteration.')
      function.warnings$rotation <- TRUE
      rotation.result <- list(scores=scores.within, loadings=loadings.within, rotation.matrix=diag(npc.within))
    }
    loadings.within <- rotation.result$loadings
    scores.within <- rotation.result$scores
    rotation.matrix.within <- rotation.result$rotation.matrix

    # Compute objective value function
    objective.function <- compute_objective_value(
      quantile.value=quantile.value,
      Y=Y,
      scores=cbind(scores.between.full, scores.within),
      intercept=intercept,
      loadings=cbind(loadings.between, loadings.within))
    objective.function.array <- c(objective.function.array, objective.function)

    convergence.criteria <- abs(objective.function.array[i] - objective.function.array[i-1])
    if(convergence.criteria < tol){convergence = TRUE}

    # Measure computation time
    loop.execution.time <- difftime(Sys.time(), loop.start.time, units = 'secs')

    if(verbose)
    {
      message(
        format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        '. Iteration ', i, ' completed in ', round(loop.execution.time, 3), ' seconds', '\n',
        'Objective function (i-1) value: ', round(objective.function.array[i-1], 4), '\n',
        'Objective function   i   value: ', round(objective.function.array[i], 4), '\n',
        'Convergence criteria value    : ', round(convergence.criteria, 4), '\n',
        '____________________________________________________________')
    }

    if(convergence){break}

    # Avoid computational issues and perform early stopping if reaching diverging loops.
    if(i>=3 && (objective.function.array[i] > 10 * objective.function.array[2]))
    {
      function.warnings$diverging.loop <- TRUE
      message('Iteration ', i, '. Breaking diverging loop')
      break
    }
  }

  if(verbose && convergence){message("\u2705 Algorithm converged succesfully")}

  if(i == max.iters && !convergence && !function.warnings$diverging.loop){warning('\u274C Algorithm reached maximum number of iterations without convergence. Consider increasing the value of max.iters parameter')}

  # Compute explained variability
  pve.between <- explained_variance_ratio(scores.between)
  pve.within <- explained_variance_ratio(scores.within)

  global.execution.time <- difftime(Sys.time(), global.start.time, units = 'secs')

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
      function.warnings = function.warnings,
      execution.time = global.execution.time,
      rotation.matrix.between = rotation.matrix.between,
      rotation.matrix.within = rotation.matrix.within,
      spline.coefficients.between = spline.coefficients.between,
      spline.coefficients.within = spline.coefficients.within,
      spline.basis = spline.basis,
      n.iters = i,
      partial.execution.time = partial.execution.time,
      inputs = inputs), class = "mfqpca_object")

  return(results)
}

# ADDITIONAL FUNCTIONS --------------------------------------------------------

#' @title Predict mfqpca scores
#' @description S3 method for class 'mfqpca_object' Given a new matrix Y, predicts the value of the scores associated to the given matrix.
#' @param object An object output of the fqpca function.
#' @param newdata The N by T matrix of observed time instants to be tested
#' @param newdata.group An N dimensional array indicating the hierarchical structure of the data. Elements in the array with the same value indicate they are repeated measures of the same individual.
#' @param ... further arguments passed to or from other methods.
#' @return List of matrices of predicted scores.
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
  if(!inherits(object, "mfqpca_object")){stop('The object must be of class mfqpca_object')}

  r <- coerce_functional_input(data=newdata, colname=object$inputs$colname, group=newdata.group)
  Y <- r$Y
  group <- r$group

  Y.mask <- !is.na(Y)
  n.obs <- nrow(Y)
  n.time <- ncol(Y)
  if(sum(!Y.mask) == n.obs * n.time){stop('newdata contains no information. Check that there are values different from NA')}

  scores.between <- mfqpca_compute_scores_between(
    Y=Y,
    Y.mask=Y.mask,
    group=group,
    intercept=object$intercept,
    loadings=object$loadings.between,
    quantile.value=object$inputs$quantile.value)
  scores.between.full <- scores.between[match(group, unique(group)), , drop=FALSE]

  Y.between.fitted <- sweep(
    scores.between.full %*% t(object$loadings.between),
    MARGIN = 2,
    STATS = object$intercept,
    FUN = "+")
  Y.within <- Y - Y.between.fitted

  scores.within <- compute_scores(
    Y=Y.within,
    Y.mask=Y.mask,
    loadings=object$loadings.within,
    quantile.value=object$inputs$quantile.value)

  return(list(scores.between=scores.between, scores.between.full=scores.between.full, scores.within=scores.within))
}

#' @title Fit Yhat
#' @description S3 method for class 'mfqpca_object'. Given an mfqpca_object model, estimates Yhat for different pve values.
#' @param object An object output of the fqpca function.
#' @param pve.between Percentage of explained variability (between 0 and 1) used to select the number of between level components in Yhat estimation. Set to NULL to use all components.
#' @param pve.within Percentage of explained variability (between 0 and 1) used to select the number of within level components in Yhat estimation. Set to NULL to use all components.
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
  if (!inherits(object, "mfqpca_object")){stop('The object must be of class mfqpca_object')}

  n.components.between <- select_npc(scores=object$scores.between, pve=pve.between)
  n.components.within <- select_npc(scores=object$scores.within, pve=pve.within)

  Y.fitted =
    reconstruct_scores_loadings(object$scores.between.full, object$loadings.between, n.components.between) +
    reconstruct_scores_loadings(object$scores.within, object$loadings.within, n.components.within)

  Y.fitted = sweep(
    Y.fitted,
    MARGIN = 2,
    STATS = object$intercept,
    FUN = "+")

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
plot.mfqpca_object <- function(x, pve.between = 0.95, pve.within = 0.95, ...)
{
  if (!inherits(x, "mfqpca_object")) stop("The x must be of class mfqpca_object")

  # Ensure zero components when pve == 0
  n.components.between <- select_npc(scores = x$scores.between, pve = pve.between)
  n.components.within  <- select_npc(scores = x$scores.within, pve = pve.within)

  intercept <- x$intercept

  # Intercept data frame
  intercept_df <- data.frame(
    Time = seq_along(intercept),
    Loading = intercept,
    Component = "Intercept",
    stringsAsFactors = FALSE)

  # Between components
  loadings.between <- x$loadings.between
  if (is.null(loadings.between) || n.components.between <= 0) {
    fqpc_between_df <- data.frame(
      Time = integer(0),
      Loading = numeric(0),
      Component = character(0),
      stringsAsFactors = FALSE)
    n.components.between <- 0
  } else {
    loadings.between.filter <- loadings.between[, seq_len(n.components.between), drop = FALSE]
    fqpc_between_df <- data.frame(
      Time = rep(seq_len(nrow(loadings.between.filter)), times = ncol(loadings.between.filter)),
      Loading = as.vector(loadings.between.filter),
      Component = rep(paste0("MFQPC Between ", seq_len(ncol(loadings.between.filter))), each = nrow(loadings.between.filter)),
      stringsAsFactors = FALSE)
  }

  # Within components
  loadings.within <- x$loadings.within
  if (is.null(loadings.within) || n.components.within <= 0) {
    fqpc_within_df <- data.frame(
      Time = integer(0),
      Loading = numeric(0),
      Component = character(0),
      stringsAsFactors = FALSE)
    n.components.within <- 0
  } else {
    loadings.within.filter <- loadings.within[, seq_len(n.components.within), drop = FALSE]
    fqpc_within_df <- data.frame(
      Time = rep(seq_len(nrow(loadings.within.filter)), times = ncol(loadings.within.filter)),
      Loading = as.vector(loadings.within.filter),
      Component = rep(paste0("MFQPC Within ", seq_len(ncol(loadings.within.filter))), each = nrow(loadings.within.filter)),
      stringsAsFactors = FALSE)
  }

  # Combine data
  plot_data <- rbind(intercept_df, fqpc_between_df, fqpc_within_df)

  # Safe factor levels (avoid 1:0 issue by using seq_len)
  between_levels <- if (n.components.between > 0) paste0("MFQPC Between ", seq_len(n.components.between)) else character(0)
  within_levels  <- if (n.components.within  > 0) paste0("MFQPC Within ", seq_len(n.components.within))  else character(0)
  plot_data$Component <- factor(plot_data$Component, levels = c("Intercept", between_levels, within_levels))

  # Plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$Time, y = .data$Loading)) +
    ggplot2::geom_line(color = "steelblue") +
    ggplot2::facet_wrap(~ Component, scales = "free_y", ncol = 4) +
    ggplot2::labs(title = "Loading Functions", x = "Time", y = "Loading") +
    ggplot2::theme_bw()

  return(p)
}

# CROSS VALIDATION ------------------------------------------------------------

#' @title Cross validation of degrees of freedom
#' @description Performs cross validation on degrees of freedom parameter of mfqpca
#' @param data An \eqn{(N \times T)} matrix, a tf object from the tidyfun package or a data.frame containing the functional data as a tf column.
#' @param group Either a string or an array. If it is a string, it must point to the grouping variable in data only if data is a dataframe. If an array, it must be the N dimensional array indicating the hierarchical structure of the data. Elements in the array with the same value indicate they are repeated measures of the same individual.
#' @param colname The name of the column containing the functional data. Use only if data is a dataframe and colname is a column in the dataframe.
#' @param npc.between The number of estimated between level components.
#' @param npc.within The number of estimated within level components.
#' @param pve.between Float between 0 and 1. Percentage of variability explained by between level components. This affects the number of components used in the curve reconstruction and error estimation. Set to NULL to avoid this behavior.
#' @param pve.within  Float between 0 and 1. Percentage of variability explained by within level components. This affects the number of components used in the curve reconstruction and error estimation. Set to NULL to avoid this behavior.
#' @param quantile.value The quantile considered.
#' @param periodic Boolean indicating if the data is expected to be periodic (start coincides with end) or not.
#' @param splines.df.grid Grid of possible values for the degrees of freedom.
#' @param splines.method Method used in the resolution of the splines quantile regression model. It currently accepts the methods \code{c('conquer', 'quantreg')}.
#' @param n.folds Number of folds to be used on cross validation.
#' @param return.models Should the list of all the models built be returned?
#' @param criteria Criteria used to divide the data. Valid values are \code{'rows'}, which considers the division based on full rows, or \code{'points'}, which considers the division based on points within the matrix.
#' @param tol Tolerance on the convergence of the algorithm.
#' @param max.iters Maximum number of iterations.
#' @param verbose.mfqpca Boolean indicating verbosity of the fqpca function.
#' @param verbose.cv Boolean indicating verbosity of the cross-validation process.
#' @param seed Seed for the random generator number.
#' @return A list containing the matrices of scores, the matrices of loadings, and a secondary list with extra information.
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
#' c2.vals * matrix(pcw, nrow = N, ncol=n.time, byrow = TRUE)
#' Y <- Y + matrix(rnorm(N*n.time, 0, 0.4), nrow = N)
#' Y[sample(N*n.time, as.integer(0.2*N))] <- NA
#'
#' cv_results <- mfqpca_cv_df(data = Y, group = group,  splines.df.grid = c(5, 10, 15), n.folds = 2)
#'
mfqpca_cv_df <- function(
    data, group, colname=NULL, npc.between = 2,  npc.within = 2, pve.between = NULL, pve.within = NULL,
    quantile.value = 0.5, n.folds = 3, return.models = TRUE, criteria = 'points', periodic = TRUE,
    splines.df.grid = c(
      5,
      10,
      15,
      20), tol = 1e-3, max.iters = 20,
    splines.method = 'conquer', verbose.mfqpca = FALSE, verbose.cv = TRUE, seed = NULL)
{
  if(!is.null(seed)){set.seed(seed)}

  if(n.folds != floor(n.folds)){stop('n.folds must be an integer number. Value provided: ', n.folds)}
  if(!(criteria %in% c('rows', 'points'))){stop('Invalid criteria. Valid criterias are c("rows", "points". Value provided: ', criteria)}
  if(!all(splines.df.grid == floor(splines.df.grid))){stop('splines.df.grid must be a positive integer array.')}

  # Check the input parameters except Y
  mfqpca_check_params(
    npc.between=npc.between,
    npc.within=npc.within,
    quantile.value=quantile.value,
    periodic=periodic,
    splines.df=max(npc.within, npc.between, splines.df.grid),
    splines.method=splines.method,
    tol=tol,
    max.iters=max.iters,
    verbose=verbose.mfqpca,
    seed=seed)

  # Check Y and colname and return an unnamed matrix
  r <- coerce_functional_input(data=data, colname=colname, group=group)
  Y <- r$Y
  group <- r$group

  fit_fun <- function(Y.train, splines.df, train.idx)
  {
    if(npc.between > splines.df)warning('The npc.between is larger than splines.df\nCurrent splines.df value: ', splines.df, '\nTaking npc.between=splines.df for this iteration of the Cross-Validation process.')
    if(npc.within > splines.df)warning('The npc.within is larger than splines.df\nCurrent splines.df value: ', splines.df, '\nTaking npc.within=splines.df for this iteration of the Cross-Validation process.')

    npc.between.iteration <- min(splines.df, npc.between)
    npc.within.iteration <- min(splines.df, npc.within)

    group.train <- if(criteria == 'rows'){group[train.idx]}else{group}

    mfqpca(
      data = Y.train,
      group = group.train,
      colname = colname,
      npc.between = npc.between.iteration,
      npc.within = npc.within.iteration,
      quantile.value = quantile.value,
      periodic = periodic,
      splines.df = splines.df,
      splines.method = splines.method,
      tol = tol,
      max.iters = max.iters,
      verbose = verbose.mfqpca,
      seed = seed)
  }

  reconstruct_fun <- function(model, Y.test, train.idx)
  {
    # Obtain optimal number of PC
    n.components.between <- select_npc(scores=model$scores.between, pve=pve.between)
    n.components.within <- select_npc(scores=model$scores.within, pve=pve.within)

    if(criteria == 'points')
    {
      scores.between.full <- model$scores.between.full
      scores.within <- model$scores.within
    }else{
      group.test <- group[-train.idx]
      test.scores <- predict.mfqpca_object(model, newdata.group=group.test, newdata=Y.test)
      scores.between.full <- test.scores$scores.between.full
      scores.within <- test.scores$scores.within
    }

    Y.pred =
      reconstruct_scores_loadings(scores.between.full, model$loadings.between, n.components.between) +
      reconstruct_scores_loadings(scores.within, model$loadings.within, n.components.within)

    sweep(
      Y.pred,
      MARGIN = 2,
      STATS = model$intercept,
      FUN = "+")
  }

  r <- cv_framework(
    Y = Y,
    grid = splines.df.grid,
    n.folds = n.folds,
    criteria = criteria,
    quantile.value = quantile.value,
    fit_fun = fit_fun,
    reconstruct_fun = reconstruct_fun,
    return.models = return.models,
    verbose.cv = verbose.cv,
    seed = seed,
    grid.label = 'Degrees of freedom: ',
    model.prefix = 'df_idx')

  return(list(
    error.matrix = r$error.matrix,
    execution.time = r$execution.time,
    splines.df.grid = splines.df.grid,
    criteria = criteria,
    list.models = r$list.models))
}
