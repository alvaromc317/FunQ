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
  global.start.time <- Sys.time()

  # Get the input parameters
  formal_names <- setdiff(names(formals(sys.function())), "...")
  inputs <- mget(formal_names, envir = environment())
  # regressors are stored at the top level of the returned object; avoid storing them twice
  inputs$regressors = NULL

  # Check the input parameters except Y and regressors
  check_positive_integer(npc, 'npc')
  check_quantile_value(quantile.value)
  check_flag(periodic, 'periodic')
  if (!is.numeric(splines.df) || length(splines.df) != 1 ||
      splines.df %% 1 != 0 || splines.df < npc) {
    stop(
      "Invalid input for 'splines.df': ",
      splines.df,
      ". Expected an integer number larger or equal than 'npc' (",
      npc,
      ").")
  }
  check_splines_method(splines.method)
  check_positive_float(tol, 'tol')
  check_iters(max.iters, 'max.iters')
  check_flag(verbose, 'verbose')
  check_seed(seed)

  # Check Y and regressors
  Y <- coerce_functional_input(data=Y)$Y
  regressors_checked <- check_regressors(regressors)
  regressors <- regressors_checked$regressors

  if(nrow(Y) != nrow(regressors)){stop('nrow(Y) and nrow(regressors) must be the same.')}

  # If seed is provided, set seed for computations
  if(!is.null(seed)){set.seed(seed)}

  n.obs <- dim(Y)[1]
  n.time <- dim(Y)[2]
  n.regressors <- ncol(regressors)
  Y.axis <- seq(0, 1, length.out = n.time)

  Y.mask <- !is.na(Y)
  Y.list <- lapply(seq_len(n.obs), function(j) Y[j, Y.mask[j, ]])
  Y.vector <- unlist(Y.list, use.names = FALSE)

  function.warnings <- list(
    splines = FALSE,
    scores = FALSE,
    diverging.loop = FALSE,
    rotation = FALSE)

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

  # fosqr LEVEL
  spline.coef.result <- try(fit_tensor_quantile_regression(
    Y.vector = Y.vector,
    Y.mask = Y.mask,
    scores = regressors,
    spline.basis = spline.basis,
    quantile.value = quantile.value,
    method=splines.method,
    return.model = TRUE), silent = FALSE)
  if(inherits(spline.coef.result, "try-error")){stop('Iteration 1. Failed computation of spline coefficients.')}
  fosqr.spline.coef <- spline.coef.result$spline.coefficients
  fosqr.loadings.list <- compute_loadings(spline.basis = spline.basis, spline.coefficients = fosqr.spline.coef)
  fosqr.intercept <- fosqr.loadings.list$intercept
  fosqr.loadings <- fosqr.loadings.list$loadings
  Y.fosqr.hat <- sweep(
    regressors %*% t(fosqr.loadings),
    MARGIN = 2,
    STATS = fosqr.intercept,
    FUN = "+")

  # FQPCA level
  fqpca.spline.coef <- matrix(stats::rnorm(dim(spline.basis)[2] * npc, mean = 0, sd = 1), nrow = dim(spline.basis)[2], ncol = npc)
  fqpca.loadings <- spline.basis %*% fqpca.spline.coef
  fqpca.scores <- try(compute_scores(
    Y = Y,
    Y.mask = Y.mask,
    loadings = fqpca.loadings,
    quantile.value = quantile.value,
    offset = Y.fosqr.hat), silent = FALSE)
  if(!is.matrix(fqpca.scores)){stop('Iteration 1. Failed computation of scores')}

  # ROTATION PROCESS
  rotation.result <- try(rotate_factors(scores = fqpca.scores, loadings = fqpca.loadings, intercept = fosqr.intercept), silent = FALSE)
  if(inherits(rotation.result, "try-error"))
  {
    warning('Iteration 1. Failed rotation process. Skipping rotation on this iteration.')
    function.warnings$rotation <- TRUE
    rotation.result <- list(
      scores = fqpca.scores,
      intercept = fosqr.intercept,
      loadings = fqpca.loadings,
      rotation.matrix = diag(npc))
  }
  fqpca.loadings <- rotation.result$loadings
  fqpca.scores <- rotation.result$scores
  fosqr.intercept <- rotation.result$intercept
  rotation.matrix = rotation.result$rotation.matrix

  # Compute objective value function
  objective.function <- compute_objective_value(
    quantile.value=quantile.value,
    Y=Y,
    scores=cbind(regressors, fqpca.scores),
    intercept=fosqr.intercept,
    loadings=cbind(fosqr.loadings, fqpca.loadings))
  objective.function.array <- objective.function

  loop.execution.time <- difftime(Sys.time(), loop.start.time, units = 'secs')

  if(verbose)
  {
    message(
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      '. Iteration ', 1, ' completed in ', round(loop.execution.time, 3), ' seconds', '\n',
      'Objective function   i   value: ', round(objective.function.array[1], 4), '\n',
      '____________________________________________________________')
  }

  convergence <- FALSE
  for(i in 2:max.iters)
  {
    loop.start.time <- Sys.time()

    # Snapshot the previous iteration state so it can be fully restored if any step fails
    fqpca.loadings.prev <- fqpca.loadings
    fqpca.scores.prev <- fqpca.scores
    fosqr.intercept.prev <- fosqr.intercept
    fosqr.loadings.prev <- fosqr.loadings
    Y.fosqr.hat.prev <- Y.fosqr.hat

    # Obtain splines coefficients
    spline.coef.result <- try(fit_tensor_quantile_regression(
      Y.vector = Y.vector,
      Y.mask = Y.mask,
      scores = cbind(regressors, fqpca.scores),
      spline.basis = spline.basis,
      quantile.value = quantile.value,
      method=splines.method,
      return.model = TRUE), silent = FALSE)
    if(inherits(spline.coef.result, "try-error"))
    {
      warning('Iteration ', i, '. Failed computation of spline coefficients. Providing results from previous iteration.')
      function.warnings$splines <- TRUE
      break
    }
    spline.coefficients <- spline.coef.result$spline.coefficients
    fqpca.spline.coef <- spline.coefficients[, (n.regressors+2):ncol(spline.coefficients), drop=FALSE]
    fosqr.spline.coef <- spline.coefficients[, seq_len(n.regressors+1), drop=FALSE]

    # Compute loadings
    fqpca.loadings <- spline.basis %*% fqpca.spline.coef
    fosqr.loadings.list <- compute_loadings(spline.basis = spline.basis, spline.coefficients = fosqr.spline.coef)
    fosqr.intercept <- fosqr.loadings.list$intercept
    fosqr.loadings <- fosqr.loadings.list$loadings
    Y.fosqr.hat <- sweep(
      regressors %*% t(fosqr.loadings),
      MARGIN = 2,
      STATS = fosqr.intercept,
      FUN = "+")

    # Compute scores
    fqpca.scores <- try(compute_scores(
      Y = Y,
      Y.mask = Y.mask,
      loadings = fqpca.loadings,
      quantile.value = quantile.value,
      offset = Y.fosqr.hat), silent = FALSE)
    if(!is.matrix(fqpca.scores))
    {
      warning('Iteration ', i, '. Failed computation of scores. Providing results from previous iteration.')
      function.warnings$scores <- TRUE
      fqpca.loadings <- fqpca.loadings.prev
      fqpca.scores <- fqpca.scores.prev
      fosqr.intercept <- fosqr.intercept.prev
      fosqr.loadings <- fosqr.loadings.prev
      Y.fosqr.hat <- Y.fosqr.hat.prev
      break
    }

    # Rotate loadings and scores
    rotation.result <- try(rotate_factors(scores = fqpca.scores, loadings = fqpca.loadings, intercept = fosqr.intercept), silent = FALSE)
    if(inherits(rotation.result, "try-error"))
    {
      warning('Iteration ', i, '. Failed rotation process. Skipping rotation on this iteration.')
      function.warnings$rotation <- TRUE
      rotation.result <- list(
        scores = fqpca.scores,
        intercept = fosqr.intercept,
        loadings = fqpca.loadings,
        rotation.matrix = diag(npc))
    }
    fqpca.loadings <- rotation.result$loadings
    fqpca.scores <- rotation.result$scores
    fosqr.intercept <- rotation.result$intercept
    rotation.matrix = rotation.result$rotation.matrix

    # Compute objective value function
    objective.function <- compute_objective_value(
      quantile.value=quantile.value,
      Y=Y,
      scores=cbind(regressors, fqpca.scores),
      intercept=fosqr.intercept,
      loadings=cbind(fosqr.loadings, fqpca.loadings))
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

    # Avoid computational issues and do early stops if reaching diverging loops.
    if(i>=3 && (objective.function.array[i] > 10 * objective.function.array[2]))
    {
      function.warnings$diverging.loop <- TRUE
      message('Iteration ', i, '. Breaking diverging loop')
      break
    }
  }

  if(verbose && convergence){message("\u2705 Algorithm converged successfully")}
  if(i == max.iters && !convergence && !function.warnings$diverging.loop){warning('\u274C Algorithm reached maximum number of iterations without convergence. Consider increasing the value of max.iters')}

  # EXPLAINED VARIABILITY
  pve <- explained_variance_ratio(fqpca.scores)

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
      inputs = inputs), class = "fosqr_fqpca_object")
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
  if (!inherits(object, "fosqr_fqpca_object")){stop('The object must be of class fosqr_fqpca_object')}

  Y <- coerce_functional_input(data=newdata)$Y
  regressors_checked <- check_regressors(newregressors)
  regressors <- regressors_checked$regressors
  Y.fosqr.hat <- sweep(
    regressors %*% t(object$fosqr.loadings),
    MARGIN = 2,
    STATS = object$fosqr.intercept,
    FUN = "+")

  n.time.model <- nrow(object$fqpca.loadings)
  if(ncol(Y) != n.time.model) {
    stop(
      "Number of time points in newdata (",
      ncol(Y),
      ") does not match the model dimension (",
      n.time.model,
      ").")
  }

  Y.mask <- !is.na(Y)
  n.obs <- nrow(Y)
  n.time <- ncol(Y)
  if(sum(!Y.mask) == n.obs * n.time){stop('newdata contains no information. Check that there are values different from NA')}

  # Estimation of scores
  scores <- try(compute_scores(
    Y = Y,
    Y.mask = Y.mask,
    loadings = object$fqpca.loadings,
    quantile.value = object$inputs$quantile.value,
    offset = Y.fosqr.hat), silent = FALSE)
  if(!is.matrix(scores))
  {
    stop('Computation of scores failed.')
  }
  return(scores)
}

#' @title Fit Yhat
#' @description S3 method for class 'fosqr_fqpca_object'. Given a fosqr_fqpca_object model, estimates Yhat for different pve values.
#' @param object An object output of the fosqr_fqpca function.
#' @param pve Percentage of explained variability (between 0 and 1) used to select the number of components in Yhat estimation. Set to NULL to use all components.
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
  if (!inherits(object, "fosqr_fqpca_object")){stop('The object must be of class fosqr_fqpca_object')}
  n.components <- select_npc(scores=object$fqpca.scores, pve=pve)

  # Build fosqr fitted value
  fosqr.fitted <- sweep(
    object$regressors %*% t(object$fosqr.loadings),
    MARGIN = 2,
    STATS = object$fosqr.intercept,
    FUN = "+")

  # Build FQPCA fitted value
  fqpca.fitted <- reconstruct_scores_loadings(object$fqpca.scores, object$fqpca.loadings, n.components)
  Y.pred <- fosqr.fitted + fqpca.fitted
  return(Y.pred)
}

# BASIC PLOT ------------------------------------------------------------------

#' @title Plot fosqr and fqpca loading functions with ggplot2
#' @description S3 method for class 'fosqr_fqpca_object'. Given a fosqr_fqpca_object model, plots the principal component and fosqr loadings.
#' @param x An object output of the fosqr_fqpca function.
#' @param pve Percentage of explained variability (between 0 and 1) used to select the number of components in Yhat estimation. Set to NULL to use all components.
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
  if (!inherits(x, "fosqr_fqpca_object")){stop('The object x must be of class fosqr_fqpca_object')}

  # Ensure zero components when pve == 0
  n.components <- select_npc(scores = x$fqpca.scores, pve = pve)

  intercept <- x$fosqr.intercept

  # Intercept data frame
  intercept_df <- data.frame(
    Time = seq_along(intercept),
    Loading = intercept,
    Component = "Intercept",
    stringsAsFactors = FALSE)

  # Between components
  fosqr.loadings <- x$fosqr.loadings
  names.regressors = x$names.regressors
  if(is.null(names.regressors)){names.regressors = paste0("B ", seq_len(ncol(fosqr.loadings)))}
  fosqr_df <- data.frame(
    Time = rep(seq_len(nrow(fosqr.loadings)), times = ncol(fosqr.loadings)),
    Loading = as.vector(fosqr.loadings),
    Component = rep(names.regressors, each = nrow(fosqr.loadings)),
    stringsAsFactors = FALSE)

  # FQPCA components
  fqpca.loadings <- x$fqpca.loadings
  if (is.null(fqpca.loadings) || n.components <= 0) {
    fqpca_df <- data.frame(
      Time = integer(0),
      Loading = numeric(0),
      Component = character(0),
      stringsAsFactors = FALSE)
    n.components <- 0
  } else {
    fqpca.loadings.filter <- fqpca.loadings[, seq_len(n.components), drop = FALSE]
    fqpca_df <- data.frame(
      Time = rep(seq_len(nrow(fqpca.loadings.filter)), times = ncol(fqpca.loadings.filter)),
      Loading = as.vector(fqpca.loadings.filter),
      Component = rep(paste0("FQPC ", seq_len(ncol(fqpca.loadings.filter))), each = nrow(fqpca.loadings.filter)),
      stringsAsFactors = FALSE)
  }

  # Combine data
  plot_data <- rbind(intercept_df, fosqr_df, fqpca_df)

  # Safe factor levels (avoid 1:0 issue by using seq_len)
  fosqr_levels <- names.regressors
  fqpca_levels  <- if (n.components  > 0) paste0("FQPC ", seq_len(n.components))  else character(0)
  plot_data$Component <- factor(plot_data$Component, levels = c("Intercept", fosqr_levels, fqpca_levels))

  # Plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$Time, y = .data$Loading)) +
    ggplot2::geom_line(color = "steelblue") +
    ggplot2::facet_wrap(~ Component, scales = "free_y", ncol = 4) +
    ggplot2::labs(title = "Loading Functions", x = "Time", y = "Loading") +
    ggplot2::theme_bw()
  return(p)
}

# INFERENCE -------------------------------------------------------------------

#' @title Compute variance associated to fosqr loadings
#' @description S3 method for class 'fosqr_fqpca_object'. Given a fosqr_fqpca_object model, estimates variances associated to fosqr loadings.
#' @param object An object output of the fosqr_fqpca function.
#' @param pve Percentage of explained variability (between 0 and 1) used to select the number of components in Yhat estimation. Set to NULL to use all components.
#' @param n.bootstrap An integer number used in the bootstrapping part of the variance estimation.
#' @param ... further arguments passed to or from other methods.
#' @return A list containing the two parts of the variance computation and the final variances.
#' @export
#' @method compute_variance fosqr_fqpca_object
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
#' variances <- compute_variance(object=results, pve=0.95)
compute_variance.fosqr_fqpca_object <- function(object, pve=0.95, n.bootstrap=1000, ...)
{
  if (!inherits(object, "fosqr_fqpca_object")){stop('The object must be of class fosqr_fqpca_object')}
  Y <- coerce_functional_input(data=object$inputs$Y)$Y
  if(nrow(Y) != nrow(object$regressors)){stop('nrow(Y) does not match nrow(regressors).')}
  if(ncol(Y) != nrow(object$fosqr.loadings)){stop('ncol(Y) does not match nrow(fosqr.loadings)')}

  n.obs <- dim(Y)[1]
  n.regressors <- ncol(object$regressors)
  Y.mask <- !is.na(Y)
  Y.list <- lapply(seq_len(n.obs), function(j) Y[j, Y.mask[j, ]])
  Y.vector <- unlist(Y.list, use.names = FALSE)
  n.components <- select_npc(scores=object$fqpca.scores, pve=pve)
  spline.coefficients.results <- try(fit_tensor_quantile_regression(
    Y.vector = Y.vector,
    Y.mask = Y.mask,
    scores = cbind(object$regressors, object$fqpca.scores[, seq_len(n.components), drop = FALSE]),
    spline.basis = object$spline.basis,
    quantile.value = object$inputs$quantile.value,
    method = object$inputs$splines.method,
    return.model = TRUE), silent = FALSE)

  spline.coefficients <- spline.coefficients.results$spline.coefficients
  B.vector <- spline.coefficients.results$B.vector
  tensor.matrix <- spline.coefficients.results$tensor.matrix

  cov.model <- get_kernel_cov(
    x = cbind(1, tensor.matrix),
    y = Y.vector,
    coefficients = B.vector,
    tau = object$inputs$quantile.value)

  var.analytical <- compute_var_sandwich(cov.model=cov.model, n.regressors=n.regressors, spline.basis=object$spline.basis)

  if (n.components == 0) {
    var.correction <- matrix(0, nrow = nrow(object$fosqr.loadings), ncol = n.regressors + 1)
  } else {
    fqpca.spline.coef <- spline.coefficients[, (n.regressors+2):ncol(spline.coefficients), drop=FALSE]
    fqpca.loadings <- (object$spline.basis %*% fqpca.spline.coef)[, seq_len(n.components), drop = FALSE]
    var.correction <- compute_var_correction(
      n.bootstrap = n.bootstrap,
      fqpca.scores = object$fqpca.scores[, seq_len(n.components), drop = FALSE],
      regressors = object$regressors,
      fqpca.loadings = fqpca.loadings)
  }
  results <- list(variance = var.analytical+var.correction, var.analytical=var.analytical, var.correction=var.correction)
  return(results)
}

# CROSS VALIDATION ------------------------------------------------------------

#' @title Cross validation of degrees of freedom
#' @description Performs cross validation on degrees of freedom parameter of fosqr_fqpca
#' @param Y An \eqn{(N \times T)} matrix or a tf object from the tidyfun package.
#' @param regressors An \eqn{(N \times P)} matrix of covariates.
#' @param npc The number of estimated components on the FQPCA part (used for modeling residual correlation)
#' @param pve Float between 0 and 1. Percentage of variability explained by components. This affects the number of components used in the curve reconstruction and error estimation. Set to NULL to avoid this behavior.
#' @param quantile.value The quantile considered.
#' @param periodic Boolean indicating if the data is expected to be periodic (start coincides with end) or not.
#' @param splines.df.grid Grid of possible values for the degrees of freedom.
#' @param splines.method Method used in the resolution of the splines quantile regression model. It currently accepts the methods \code{c('conquer', 'quantreg')}.
#' @param n.folds Number of folds to be used on cross validation.
#' @param return.models Should the list of all the models built be returned?
#' @param criteria Criteria used to divide the data. Valid values are \code{'rows'}, which considers the division based on full rows, or \code{'points'}, which considers the division based on points within the matrix.
#' @param tol Tolerance on the convergence of the algorithm.
#' @param max.iters Maximum number of iterations.
#' @param verbose.fosqr_fqpca Boolean indicating verbosity of the fosqr-fqpca function.
#' @param verbose.cv Boolean indicating verbosity of the cross-validation process.
#' @param seed Seed for the random generator number.
#' @return A list containing the matrices of scores, the matrices of loadings, and a secondary list with extra information.
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
#' cv_results <- fosqr_fqpca_cv_df(
#'     Y = Y, regressors = regressors,
#'     splines.df.grid = c(5, 10, 15), n.folds = 2)
#'
fosqr_fqpca_cv_df <- function(
    Y = NULL,
    regressors = NULL,
    npc = 2,
    pve = 0.95,
    quantile.value = 0.5,
    periodic = TRUE,
    splines.df.grid = c(5, 10, 15),
    splines.method = 'conquer',
    n.folds = 3,
    return.models = TRUE,
    criteria = 'points',
    tol = 1e-3,
    max.iters = 20,
    verbose.fosqr_fqpca = FALSE,
    verbose.cv = TRUE,
    seed = NULL)
{
  if(!is.null(seed)){set.seed(seed)}

  if(n.folds != floor(n.folds)){stop('n.folds must be an integer number. Value provided: ', n.folds)}
  if(!(criteria %in% c('rows', 'points'))){stop('Invalid criteria. Valid criterias are c("rows", "points". Value provided: ', criteria)}
  if(!all(splines.df.grid == floor(splines.df.grid))){stop('splines.df.grid must be a positive integer array.')}

  # Check the input parameters except Y and regressors
  check_positive_integer(npc, 'npc')
  check_quantile_value(quantile.value)
  check_flag(periodic, 'periodic')
  check_splines_method(splines.method)
  check_positive_float(tol, 'tol')
  check_iters(max.iters, 'max.iters')
  check_flag(verbose.fosqr_fqpca, 'verbose')
  check_seed(seed)

  # Check Y and regressors
  Y <- coerce_functional_input(data=Y)$Y
  regressors_checked <- check_regressors(regressors)
  regressors <- regressors_checked$regressors

  fit_fun <- function(Y.train, splines.df, train.idx)
  {
    if(npc > splines.df)warning('The npc is larger than splines.df\nCurrent splines.df value: ', splines.df, '\nTaking npc=splines.df for this iteration of the Cross-Validation process.')
    npc.iteration <- min(splines.df, npc)

    regressors.train <- if(criteria == 'rows'){regressors[train.idx, , drop = FALSE]}else{regressors}

    fosqr_fqpca(
      Y = Y.train,
      regressors = regressors.train,
      npc = npc.iteration,
      quantile.value = quantile.value,
      periodic = periodic,
      splines.df = splines.df,
      splines.method = splines.method,
      tol = tol,
      max.iters = max.iters,
      verbose = verbose.fosqr_fqpca,
      seed = seed)
  }

  reconstruct_fun <- function(model, Y.test, train.idx)
  {
    # Obtain optimal number of PC
    npc.reconstruction <- select_npc(scores=model$fqpca.scores, pve=pve)

    regressors.test <- if(criteria == 'rows'){regressors[-train.idx, , drop = FALSE]}else{regressors}

    if(criteria == 'points')
    {
      scores <- model$fqpca.scores
    }else{
      scores <- predict.fosqr_fqpca_object(model, newdata=Y.test, newregressors=regressors.test)
    }

    # Build fosqr fitted value
    fosqr.fitted <- sweep(
      regressors.test %*% t(model$fosqr.loadings),
      MARGIN = 2,
      STATS = model$fosqr.intercept,
      FUN = "+")

    # Build FQPCA fitted value
    fqpca.fitted <- reconstruct_scores_loadings(scores, model$fqpca.loadings, npc.reconstruction)
    fosqr.fitted + fqpca.fitted
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
