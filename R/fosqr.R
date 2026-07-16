# MAIN ------------------------------------------------------------------------

#' @title FOSQR main function
#' @description Solves the Function on Scalar Quantile Regression (FOSQR) model.
#' @param Y An \eqn{(N \times T)} matrix or a tf object from the tidyfun package.
#' @param regressors An \eqn{(N \times P)} matrix of covariates.
#' @param quantile.value The quantile considered.
#' @param periodic Boolean indicating if the data is expected to be periodic (start coincides with end) or not.
#' @param splines.df Degrees of freedom for the splines.
#' @param splines.method Method used in the resolution of the splines quantile regression model. It currently accepts the methods \code{c('conquer', 'quantreg')}.
#' @param verbose Boolean indicating the verbosity.
#' @param seed Seed for the random generator number.
#' @return fosqr_object
#' @export
fosqr <- function(
    Y = NULL,
    regressors = NULL,
    quantile.value = 0.5,
    periodic = TRUE,
    splines.df = 10,
    splines.method = 'conquer',
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
  check_quantile_value(quantile.value)
  check_flag(periodic, 'periodic')
  check_positive_integer(splines.df, 'splines.df')
  check_splines_method(splines.method)
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

  function.warnings <- list(splines = FALSE)

  spline.basis <- pbs::pbs(
    Y.axis,
    degree = 3,
    df = splines.df,
    intercept = TRUE,
    periodic=periodic,
    Boundary.knots = c(min(Y.axis), max(Y.axis)))
  n.basis <- ncol(spline.basis)

  # FOSQR LEVEL
  spline.coef.result <- try(fit_tensor_quantile_regression(
    Y.vector = Y.vector,
    Y.mask = Y.mask,
    scores = regressors,
    spline.basis = spline.basis,
    quantile.value = quantile.value,
    method = splines.method,
    return.model = TRUE), silent = FALSE)

  if(inherits(spline.coef.result, "try-error")) {
    function.warnings$splines <- TRUE
    stop("Failed computation of spline coefficients.")
  }

  fosqr.loadings.list <- compute_loadings(spline.basis = spline.basis, spline.coefficients = spline.coef.result$spline.coefficients)

  fosqr.intercept <- fosqr.loadings.list$intercept
  fosqr.loadings <- fosqr.loadings.list$loadings

  if(verbose) { message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " FOSQR computation completed.") }

  # FINAL RESULTS
  global.execution.time <- difftime(Sys.time(), global.start.time, units = 'secs')

  results <- structure(list(
      fosqr.intercept = fosqr.intercept,
      fosqr.loadings = fosqr.loadings,
      regressors = regressors,
      execution.time = global.execution.time,
      function.warnings = function.warnings,
      spline.basis = spline.basis,
      splines.model = spline.coef.result$model,
      tensor.matrix = spline.coef.result$tensor.matrix,
      names.regressors = regressors_checked$names.regressors,
      inputs = inputs), class = "fosqr_object")

  return(results)
}

# PREDICTIONS -----------------------------------------------------------------

#' @title Predict FOSQR values
#' @description S3 method for class 'fosqr_object'. Predicts the fitted curves associated to new regressors.
#' @param object An object output of the fosqr function.
#' @param newregressors The N by p matrix of regressors to be tested.
#' @param ... further arguments passed to or from other methods.
#' @return The predicted matrix of fitted curves.
#' @export
predict.fosqr_object <- function(object, newregressors, ...)
{
  if (!inherits(object, "fosqr_object")){stop('The object must be of class fosqr_object')}

  regressors_checked <- check_regressors(newregressors)
  regressors <- regressors_checked$regressors

  if (ncol(regressors) != ncol(object$regressors)) {
    stop(
      "Number of columns in newregressors (",
      ncol(regressors),
      ") does not match the model regressors (",
      ncol(object$regressors),
      ").")
  }

  Y.fosqr.hat <- sweep(
    regressors %*% t(object$fosqr.loadings),
    MARGIN = 2,
    STATS = object$fosqr.intercept,
    FUN = "+")
  return(Y.fosqr.hat)
}

#' @title Fit Yhat
#' @description S3 method for class 'fosqr_object'. Given a fosqr_object model, extracts the fitted curves.
#' @param object An object output of the fosqr function.
#' @param ... further arguments passed to or from other methods.
#' @return The matrix of fitted values.
#' @export
fitted.fosqr_object <- function(object, ...)
{
  if (!inherits(object, "fosqr_object")){stop('The object must be of class fosqr_object')}

  # Build fosqr fitted value
  Y.pred <- sweep(
    object$regressors %*% t(object$fosqr.loadings),
    MARGIN = 2,
    STATS = object$fosqr.intercept,
    FUN = "+")
  return(Y.pred)
}

#' @title Compute variance
#' @description S3 method for class 'fosqr_object'. Given a fosqr_object model, estimates the variance associated to the FOSQR loadings.
#' @param object An object output of the fosqr function.
#' @param ... further arguments passed to or from other methods.
#' @return A list containing the estimated analytical variance and the kernel covariance model.
#' @export
#' @method compute_variance fosqr_object
compute_variance.fosqr_object <- function(object, ...)
{
  if (!inherits(object, "fosqr_object")){stop('The object must be of class fosqr_object')}

  coefficients <- extract_coefficients(object$splines.model)

  Y <- object$inputs$Y
  n.obs <- nrow(Y)
  Y.mask <- !is.na(Y)
  Y.list <- lapply(seq_len(n.obs), function(j) Y[j, Y.mask[j, ]])
  Y.vector <- unlist(Y.list, use.names = FALSE)

  cov.model <- get_kernel_cov(
    x = cbind(1, object$tensor.matrix),
    y = Y.vector,
    coefficients = coefficients,
    tau = object$inputs$quantile.value)

  var.analytical <- compute_var_sandwich(n.regressors=ncol(object$regressors), cov.model=cov.model, spline.basis=object$spline.basis)

  fosqr.variances <- list(variance = var.analytical, var.analytical = var.analytical, cov.model = cov.model)
  return(fosqr.variances)
}

# BASIC PLOT ------------------------------------------------------------------

#' @title Plot fosqr loading functions with ggplot2
#' @description S3 method for class 'fosqr_object'. Given a fosqr_object model, plots the fosqr loadings.
#' @param x An object output of the fosqr function.
#' @param ... further arguments passed to or from other methods.
#' @return A ggplot object plotting the intercept and FOSQR loading curves.
#' @export
plot.fosqr_object <- function(x, ...)
{
  if (!inherits(x, "fosqr_object")){stop('The object x must be of class fosqr_object')}

  intercept <- x$fosqr.intercept

  # Intercept data frame
  intercept_df <- data.frame(
    Time = seq_along(intercept),
    Loading = intercept,
    Component = "Intercept",
    stringsAsFactors = FALSE)

  # FOSQR components
  fosqr.loadings <- x$fosqr.loadings
  names.regressors = x$names.regressors
  if(is.null(names.regressors)){names.regressors = paste0("B ", seq_len(ncol(fosqr.loadings)))}

  fosqr_df <- data.frame(
    Time = rep(seq_len(nrow(fosqr.loadings)), times = ncol(fosqr.loadings)),
    Loading = as.vector(fosqr.loadings),
    Component = rep(names.regressors, each = nrow(fosqr.loadings)),
    stringsAsFactors = FALSE)

  # Combine data
  plot_data <- rbind(intercept_df, fosqr_df)

  # Safe factor levels
  fosqr_levels <- names.regressors
  plot_data$Component <- factor(plot_data$Component, levels = c("Intercept", fosqr_levels))

  # Plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$Time, y = .data$Loading)) +
    ggplot2::geom_line(color = "steelblue") +
    ggplot2::facet_wrap(~ Component, scales = "free_y", ncol = 4) +
    ggplot2::labs(title = "FOSQR Loading Functions", x = "Time", y = "Loading") +
    ggplot2::theme_bw()
  return(p)
}
