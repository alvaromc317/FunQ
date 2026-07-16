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
    data,
    colname = NULL,
    npc = 2,
    quantile.value = 0.5,
    periodic = TRUE,
    splines.df = 10,
    splines.method = 'conquer',
    penalized = FALSE,
    lambda.ridge = 0,
    tol = 1e-3,
    max.iters = 20,
    verbose = FALSE,
    seed = NULL)
{
  global.start.time <- Sys.time()

  # Get the input parameters
  formal_names <- setdiff(names(formals(sys.function())), "...")
  inputs <- mget(formal_names, envir = environment())
  inputs$data = NULL

  # Check the input parameters except Y and colname
  check_fqpca_params(
    npc=npc,
    quantile.value=quantile.value,
    periodic=periodic,
    splines.df=splines.df,
    splines.method=splines.method,
    penalized=penalized,
    lambda.ridge=lambda.ridge,
    tol=tol,
    max.iters=max.iters,
    verbose=verbose,
    seed=seed)

  # Check Y and colname and return an unnamed matrix
  Y <- coerce_functional_input(data=data, colname=colname)$Y

  # Check if splines.method is valid
  if(penalized && splines.method != 'conquer')
  {
    stop("Invalid input for 'splines.method and/or penalized':
         If penalized is set to  TRUE then splines.method must be 'conquer'")
  }

  # If seed is provided, set seed for computations
  if(!is.null(seed)){set.seed(seed)}

  # Step 1: Initialize values
  n.obs <- nrow(Y)
  n.time <- dim(Y)[2]
  Y.axis <- seq(0, 1, length.out = n.time)
  Y.mask <- !is.na(Y)

  # vectorize Y [Y11, ..., Y1T, Y21, ..., Y2T,...YNT]
  Y.list <- lapply(seq_len(n.obs), function(j) Y[j, Y.mask[j, ]])
  Y.vector <- unlist(Y.list, use.names = FALSE)
  function.warnings <- list(
    splines = FALSE,
    scores = FALSE,
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

  # If penalized, compute different basis with an identity second derivative matrix with 0s in first and second positions
  # Based on Garcia de la Garza et al. (2024) (Adaptive Functional Principal Component Analysis)
  if(penalized)
  {
    eigen.decomp <- eigen(crossprod(diff(spline.basis, differences = 2)))
    first.part <- eigen.decomp$vectors[, (n.basis-1):n.basis]
    safe_eigen_values <- pmax(eigen.decomp$values[seq_len(n.basis-2)], 1e-8) # or a suitable small constant
    second.part <- eigen.decomp$vectors[, seq_len(n.basis-2)] %*% diag(1 / sqrt(safe_eigen_values))
    U <- cbind(first.part, second.part)
    spline.basis <- spline.basis %*% U
  }

  i <- 1
  loop.start.time <- Sys.time()

  # Randomly generate splines coefficients
  spline.coefficients <- matrix(stats::rnorm(n.basis * (npc+1), mean = 0, sd = 1), nrow=n.basis, ncol=npc+1)

  loadings_list <- spline.basis %*% spline.coefficients
  intercept <- loadings_list[, 1]
  loadings <- loadings_list[, -1, drop=FALSE]

  scores <- try(compute_scores(
    Y = Y,
    Y.mask = Y.mask,
    loadings = loadings,
    quantile.value = quantile.value,
    offset = intercept), silent = FALSE)
  if(!is.matrix(scores)){stop('Iteration 1. Failed computation of scores')}

  # Rotate loadings and scores
  rotation.result <- try(rotate_factors(scores = scores, loadings = loadings, intercept = intercept), silent = FALSE)
  if(inherits(rotation.result, "try-error"))
  {
    warning('Iteration 1. Failed rotation process. Skipping rotation on this iteration.')
    function.warnings$rotation <- TRUE
    rotation.result <- list(
      scores = scores,
      intercept = intercept,
      loadings = loadings,
      rotation.matrix = diag(npc))
  }
  scores <- rotation.result$scores
  intercept <- rotation.result$intercept
  loadings <- rotation.result$loadings
  rotation.matrix = rotation.result$rotation.matrix

  # Compute objective value function
  objective.function <- compute_objective_value(
    quantile.value = quantile.value,
    Y = Y,
    scores = scores,
    intercept = intercept,
    loadings = loadings)
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

  # Iterate until convergence or until max.iters is reached
  convergence <- FALSE
  for(i in 2:max.iters)
  {
    loop.start.time <- Sys.time()

    # Snapshot the previous iteration state so it can be fully restored if any step fails
    spline.coefficients.prev <- spline.coefficients
    intercept.prev <- intercept
    loadings.prev <- loadings
    scores.prev <- scores

    # OBTAIN SPLINE COEFFICIENTS
    spline.coefficients <- try(fit_tensor_quantile_regression(
      Y.vector = Y.vector,
      Y.mask = Y.mask,
      scores = scores,
      spline.basis = spline.basis,
      quantile.value = quantile.value,
      method = splines.method,
      penalized = penalized,
      lambda.ridge = lambda.ridge), silent = FALSE)
    if(!is.matrix(spline.coefficients))
    {
      warning('Iteration ', i, '. Failed computation of spline coefficients. Providing results from previous iteration.')
      function.warnings$splines <- TRUE
      spline.coefficients <- spline.coefficients.prev
      break
    }

    loadings_list <- compute_loadings(spline.basis=spline.basis, spline.coefficients=spline.coefficients)
    intercept <- loadings_list$intercept
    loadings <- loadings_list$loadings

    # OBTAIN SCORES
    scores <- try(compute_scores(
      Y = Y,
      Y.mask = Y.mask,
      loadings = loadings,
      quantile.value = quantile.value,
      offset = intercept), silent = FALSE)
    if(!is.matrix(scores))
    {
      warning('Iteration ', i, '. Failed computation of scores. Providing results from previous iteration.')
      function.warnings$scores <- TRUE
      spline.coefficients <- spline.coefficients.prev
      intercept <- intercept.prev
      loadings <- loadings.prev
      scores <- scores.prev
      break
    }

    # Rotate loadings and scores
    rotation.result <- try(rotate_factors(scores = scores, loadings = loadings, intercept = intercept), silent = FALSE)
    if(inherits(rotation.result, "try-error"))
    {
      warning('Iteration ', i, '. Failed rotation process. Skipping rotation on this iteration.')
      function.warnings$rotation <- TRUE
      rotation.result <- list(
        scores = scores,
        intercept = intercept,
        loadings = loadings,
        rotation.matrix = diag(npc))
    }
    scores <- rotation.result$scores
    intercept <- rotation.result$intercept
    loadings <- rotation.result$loadings
    rotation.matrix = rotation.result$rotation.matrix

    # Compute objective function
    objective.function <- compute_objective_value(
      quantile.value = quantile.value,
      Y = Y,
      scores = scores,
      intercept = intercept,
      loadings = loadings)
    objective.function.array <- c(objective.function.array, objective.function)

    convergence.criteria <- abs(objective.function.array[i] - objective.function.array[i-1])
    if(convergence.criteria < tol){convergence = TRUE}

    # Measure computation time
    loop.execution.time <- difftime(Sys.time(), loop.start.time, units = 'secs')

    if(verbose)
    {
      message(
        format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        '. Iteration ', i, ' completed in ',  round(loop.execution.time, 3), ' seconds', '\n',
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

  if(verbose && convergence){message("\u2705 Algorithm converged succesfully")}

  if(i == max.iters && !convergence && !function.warnings$diverging.loop){warning('\u274C Algorithm reached maximum number of iterations without convergence. Consider increasing the value of max.iters parameter')}

  # EXPLAINED VARIABILITY -----------------------------------------------------

  pve <- explained_variance_ratio(scores)

  global.execution.time <- difftime(Sys.time(), global.start.time, units = 'secs')

  results <- structure(list(
      scores = scores,
      intercept = intercept,
      loadings = loadings,
      pve = pve,
      objective.function.value = objective.function,
      objective.function.array = objective.function.array,
      function.warnings = function.warnings,
      execution.time = global.execution.time,
      rotation.matrix = rotation.matrix,
      spline.coefficients = spline.coefficients,
      spline.basis = spline.basis,
      n.iters = i,
      inputs = inputs), class = "fqpca_object")

  return(results)
}

# PREDICTIONS -----------------------------------------------------------------

#' @title Predict fqpca scores
#' @description S3 method for class 'fqpca_object' Given a new matrix Y, predicts the value of the scores associated to the given matrix.
#' @param object An object output of the fqpca function.
#' @param newdata The N by T matrix of observed time instants to be tested, or the dataframe storing the tf functional vector using the same colname as the one used in the fqpca function.
#' @param ... further arguments passed to or from other methods.
#' @return The predicted matrix of scores.
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
  if (!inherits(object, "fqpca_object")){stop('The object must be of class fqpca_object')}

  Y <- coerce_functional_input(data=newdata, colname=object$inputs$colname)$Y

  Y.mask <- !is.na(Y)
  n.obs <- nrow(Y)
  n.time <- ncol(Y)
  if(sum(!Y.mask) == n.obs * n.time){stop('newdata contains no information. Check that there are values different from NA')}

  # Estimation of scores
  scores <- try(compute_scores(
    Y = Y,
    Y.mask = Y.mask,
    loadings = object$loadings,
    quantile.value = object$inputs$quantile.value,
    offset = object$intercept), silent = FALSE)
  if(!is.matrix(scores))
  {
    stop('Computation of scores failed.')
  }
  return(scores)
}

#' @title Fit Yhat
#' @description S3 method for class 'fqpca_object'. Given an fqpca_object model, estimates Yhat for different pve values.
#' @param object An object output of the fqpca function.
#' @param pve Percentage of explained variability (between 0 and 1) used to select the number of components in Yhat estimation. Set to NULL to use all components.
#' @param ... further arguments passed to or from other methods.
#' @return The matrix of fitted values.
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
  if (!inherits(object, "fqpca_object")){stop('The object must be of class fqpca_object')}
  n.components <- select_npc(scores=object$scores, pve=pve)
  Y.pred <- sweep(
    reconstruct_scores_loadings(object$scores, object$loadings, n.components),
    MARGIN = 2,
    STATS = object$intercept,
    FUN = "+")
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
  n.components <- select_npc(scores=x$scores, pve=pve)

  intercept <- x$intercept
  loadings <- x$loadings

  # Create a data frame for the intercept curve
  intercept_df <- data.frame(Time = seq_len(length(intercept)), Loading = intercept, Component = "Intercept")

  # Handle the case when n.components = 0 (only plot intercept)
  if (n.components == 0) {
    plot_data <- intercept_df
    plot_data$Component <- factor(plot_data$Component, levels = "Intercept")
  } else {
    non_int <- loadings[, seq_len(n.components), drop = FALSE]

    # Reshape non-intercept loadings into a long format data frame
    fqpc_df <- data.frame(Time = rep(seq_len(nrow(non_int)), times = ncol(non_int)), Loading = as.vector(non_int), Component = rep(paste0("FQPC ", seq_len(ncol(non_int))), each = nrow(non_int)))

    # Combine the intercept and FQPC data frames
    plot_data <- rbind(intercept_df, fqpc_df)
    plot_data$Component <- factor(plot_data$Component, levels = c('Intercept', paste0('FQPC ', seq_len(ncol(non_int)))))
  }

  # Create the GGPlot object with facets displaying each component separately.
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$Time, y = .data$Loading)) +
    ggplot2::geom_line(color = "steelblue") +
    ggplot2::facet_wrap(~ Component, scales = "free_y", ncol = 3) +
    ggplot2::labs(title = "Loading Functions", x = "Time", y = "Loading") +
    ggplot2::theme_bw()
  return(p)
}

# CROSS VALIDATION ------------------------------------------------------------

#' @title Cross validation of lambda.ridge
#' @description Performs cross validation on lambda parameter of fqpca. Only valid if method is set to conquer
#' @param data An \eqn{(N \times T)} matrix, a tf object from the tidyfun package or a data.frame containing the functional data as a tf column.
#' @param colname The name of the column containing the functional data. Use only if data is a dataframe and colname is a column in the dataframe.
#' @param npc The number of estimated components.
#' @param pve Float between 0 and 1. Percentage of variability explained by components. This affects the number of components used in the curve reconstruction and error estimation. Set to NULL to avoid this behavior.
#' @param quantile.value The quantile considered.
#' @param periodic Boolean indicating if the data is expected to be periodic (start coincides with end) or not.
#' @param splines.df Degrees of freedom for the splines.
#' @param splines.method Method used in the resolution of the splines quantile regression model. It currently accepts the methods \code{c('conquer', 'quantreg')}.
#' @param penalized Boolean indicating if the smoothness should be controlled using a second derivative penalty. This functionality is experimental.
#' @param lambda.grid  Grid of hyper parameter values controlling the penalization on the second derivative of the splines. It has effect only with \code{penalized=TRUE} and \code{method='conquer'}.
#' @param n.folds Number of folds to be used on cross validation.
#' @param return.models Should the list of all the models built be returned?
#' @param criteria Criteria used to divide the data. Valid values are \code{'rows'}, which considers the division based on full rows, or \code{'points'}, which considers the division based on points within the matrix.
#' @param tol Tolerance on the convergence of the algorithm.
#' @param max.iters Maximum number of iterations.
#' @param verbose.fqpca Boolean indicating verbosity of the fqpca function.
#' @param verbose.cv Boolean indicating verbosity of the cross-validation process.
#' @param seed Seed for the random generator number.
#' @return A list containing the matrix of scores, the matrix of loadings, a list with all the trained models (if the return_models param is TRUE) and a secondary list with extra information.
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
#' pc1 = sin(seq(0, 2*pi, length.out = 144))
#' pc2 = cos(seq(0, 2*pi, length.out = 144))
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
#' cv_result <- fqpca_cv_lambda(data=Y, lambda.grid = c(0, 1e-6), n.folds = 2)
#'
fqpca_cv_lambda <- function(
    data, colname=NULL, npc = 2,  pve=NULL, quantile.value = 0.5,
    lambda.grid  =  c(0, 1e-10, 1e-5), n.folds = 3, return.models=TRUE,
    criteria = 'points', periodic = TRUE, splines.df = 10, tol = 1e-3,
    max.iters = 20, splines.method = 'conquer', penalized=TRUE, verbose.fqpca = FALSE,
    verbose.cv = TRUE, seed = NULL)
{
  if(!is.null(seed)){set.seed(seed)}

  if(splines.method != 'conquer'){stop('For a penalized second derivative approach, the splines.method must be conquer')}
  if(!penalized){stop('For a penalized second derivative approach, the penalized parameter must be set to TRUE.')}
  if(n.folds != floor(n.folds)){stop('n.folds must be an integer number. Value provided: ', n.folds)}
  if(!(criteria %in% c('rows', 'points'))){stop('Invalid criteria. Valid criterias are c("rows", "points". Value provided: ', criteria)}
  if(!is.numeric(lambda.grid) || length(lambda.grid) < 1 || any(lambda.grid < 0)){stop('lambda.grid must be a numeric array of positive numbers or zero.')}

  # Check the input parameters except Y and colname
  check_fqpca_params(
    npc=npc,
    quantile.value=quantile.value,
    periodic=periodic,
    splines.df=splines.df,
    splines.method=splines.method,
    penalized=penalized,
    lambda.ridge=max(lambda.grid),
    tol=tol,
    max.iters=max.iters,
    verbose=verbose.fqpca,
    seed=seed)

  # Check Y and colname and return an unnamed matrix
  Y <- coerce_functional_input(data=data, colname=colname)$Y

  fit_fun <- function(Y.train, lambda.ridge, train.idx)
  {
    fqpca(
      data = Y.train,
      npc = npc,
      quantile.value = quantile.value,
      periodic = periodic,
      splines.df = splines.df,
      splines.method = splines.method,
      penalized = TRUE,
      lambda.ridge = lambda.ridge,
      tol = tol,
      max.iters = max.iters,
      verbose = verbose.fqpca,
      seed = seed)
  }

  reconstruct_fun <- function(model, Y.test, train.idx)
  {
    npc.reconstruction <- select_npc(scores=model$scores, pve=pve)
    if(npc.reconstruction == 0){stop('pve cannot be 0.')}
    if(criteria == 'points')
    {
      scores <- model$scores
    }else{
      scores <- predict.fqpca_object(model, Y.test)
    }
    sweep(
      reconstruct_scores_loadings(scores, model$loadings, npc.reconstruction),
      MARGIN = 2,
      STATS = model$intercept,
      FUN = "+")
  }

  r <- cv_framework(
    Y = Y,
    grid = lambda.grid,
    n.folds = n.folds,
    criteria = criteria,
    quantile.value = quantile.value,
    fit_fun = fit_fun,
    reconstruct_fun = reconstruct_fun,
    return.models = return.models,
    verbose.cv = verbose.cv,
    seed = seed,
    grid.label = 'lambda=',
    model.prefix = 'lambda_idx')

  return(list(
    error.matrix = r$error.matrix,
    execution.time = r$execution.time,
    lambda.grid = lambda.grid,
    criteria = criteria,
    list.models = r$list.models))
}

#' @title Cross validation of degrees of freedom
#' @description Performs cross validation on degrees of freedom parameter of fqpca
#' @param data An \eqn{(N \times T)} matrix, a tf object from the tidyfun package or a data.frame containing the functional data as a tf column.
#' @param colname The name of the column containing the functional data. Use only if data is a dataframe and colname is a column in the dataframe.
#' @param npc The number of estimated components.
#' @param pve Float between 0 and 1. Percentage of variability explained by components. This affects the number of components used in the curve reconstruction and error estimation. Set to NULL to avoid this behavior.
#' @param quantile.value The quantile considered.
#' @param periodic Boolean indicating if the data is expected to be periodic (start coincides with end) or not.
#' @param splines.df.grid Grid of possible values for the degrees of freedom.
#' @param splines.method Method used in the resolution of the splines quantile regression model. It currently accepts the methods \code{c('conquer', 'quantreg')}.
#' @param penalized Boolean indicating if the smoothness should be controlled using a second derivative penalty. This functionality is experimental.
#' @param lambda.ridge  Hyper parameter controlling the penalization on the second derivative of the splines. It has effect only with \code{penalized=TRUE} and \code{method='conquer'}.
#' @param n.folds Number of folds to be used on cross validation.
#' @param return.models Should the list of all the models built be returned?
#' @param criteria Criteria used to divide the data. Valid values are \code{'rows'}, which considers the division based on full rows, or \code{'points'}, which considers the division based on points within the matrix.
#' @param tol Tolerance on the convergence of the algorithm.
#' @param max.iters Maximum number of iterations.
#' @param verbose.fqpca Boolean indicating verbosity of the fqpca function.
#' @param verbose.cv Boolean indicating verbosity of the cross-validation process.
#' @param seed Seed for the random generator number.
#' @return A list containing the matrix of scores, the matrix of loadings, and a secondary list with extra information.
#' @export
#' @examples
#' n.obs = 150
#' n.time = 144
#'
#' # Generate scores
#' c1.vals = rnorm(n.obs)
#' c2.vals = rnorm(n.obs)
#'
#' # Generate pc's
#' pc1 = sin(seq(0, 2*pi, length.out = 144))
#' pc2 = cos(seq(0, 2*pi, length.out = 144))
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
#' cv_result <- fqpca_cv_df(data=Y, splines.df.grid = c(5, 10, 15), n.folds = 2)
fqpca_cv_df <- function(
    data, colname=NULL, npc = 2,  pve=NULL, quantile.value = 0.5,
    lambda.ridge = 0, n.folds = 3, return.models = TRUE, criteria = 'points',
    periodic = TRUE, splines.df.grid = c(
      5,
      10,
      15,
      20), tol = 1e-3,
    max.iters = 20, splines.method = 'conquer', penalized=FALSE, verbose.fqpca = FALSE,
    verbose.cv = TRUE, seed = NULL)
{
  if(!is.null(seed)){set.seed(seed)}

  if(n.folds != floor(n.folds)){stop('n.folds must be an integer number. Value provided: ', n.folds)}
  if(!(criteria %in% c('rows', 'points'))){stop('Invalid criteria. Valid criterias are c("rows", "points". Value provided: ', criteria)}
  if(!all(splines.df.grid == floor(splines.df.grid))){stop('splines.df.grid must be a positive integer array.')}

  # Check the input parameters except Y and colname
  check_fqpca_params(
    npc=npc,
    quantile.value=quantile.value,
    periodic=periodic,
    splines.df=max(npc, splines.df.grid),
    splines.method=splines.method,
    penalized=penalized,
    lambda.ridge=lambda.ridge,
    tol=tol,
    max.iters=max.iters,
    verbose=verbose.fqpca,
    seed=seed)

  # Check Y and colname and return an unnamed matrix
  Y <- coerce_functional_input(data=data, colname=colname)$Y

  fit_fun <- function(Y.train, splines.df, train.idx)
  {
    if(npc > splines.df)warning('The npc is larger than splines.df\nCurrent splines.df value: ', splines.df, '\nTaking npc=splines.df for this iteration of the Cross-Validation process.')
    npc.iteration <- min(splines.df, npc)
    fqpca(
      data = Y.train,
      npc = npc.iteration,
      quantile.value = quantile.value,
      periodic = periodic,
      splines.df = splines.df,
      splines.method = splines.method,
      penalized=penalized,
      lambda.ridge = lambda.ridge,
      tol = tol,
      max.iters = max.iters,
      verbose = verbose.fqpca,
      seed = seed)
  }

  reconstruct_fun <- function(model, Y.test, train.idx)
  {
    npc.reconstruction <- select_npc(scores=model$scores, pve=pve)
    if(npc.reconstruction == 0){stop('pve cannot be 0.')}
    if(criteria == 'points')
    {
      scores <- model$scores
    }else{
      scores <- predict.fqpca_object(model, Y.test)
    }
    sweep(
      reconstruct_scores_loadings(scores, model$loadings, npc.reconstruction),
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
