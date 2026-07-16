
# COMPOSABLE PARAMETER VALIDATORS ----------------------------------------------

#' @title Check that a parameter is a positive integer
#' @param x The value to check.
#' @param name The parameter name used in the error message.
#' @return No return
check_positive_integer <- function(x, name)
{
  if (!is.numeric(x) || length(x) != 1 || x %% 1 != 0 || x <= 0) {
    stop(
      "Invalid input for '",
      name,
      "': ",
      x,
      ". Expected a positive integer number.")
  }
}

#' @title Check that a maximum-iterations parameter is an integer >= 2
#' @param x The value to check.
#' @param name The parameter name used in the error message.
#' @return No return
check_iters <- function(x, name)
{
  if (!is.numeric(x) || length(x) != 1 || x %% 1 != 0 || x < 2) {
    stop(
      "Invalid input for '",
      name,
      "': ",
      x,
      ". Expected an integer number >= 2.")
  }
}

#' @title Check that the quantile value lies in (0, 1)
#' @param x The value to check.
#' @return No return
check_quantile_value <- function(x)
{
  if (!is.numeric(x) || length(x) != 1 || x <= 0 || x >= 1) {
    stop("Invalid input for 'quantile.value': ", x, ". Expected a float number in (0, 1).")
  }
}

#' @title Check that a parameter is a single Boolean
#' @param x The value to check.
#' @param name The parameter name used in the error message.
#' @return No return
check_flag <- function(x, name)
{
  if (!is.logical(x) || length(x) != 1) {
    stop(
      "Invalid input for '",
      name,
      "': ",
      x,
      ". Expected a Boolean value (TRUE or FALSE).")
  }
}

#' @title Check that the splines method is supported
#' @param x The value to check.
#' @return No return
check_splines_method <- function(x)
{
  if (!is.character(x) || length(x) != 1 || !x %in% c("conquer", "quantreg")) {
    stop("Invalid input for 'splines.method': '", x, "'. Expected 'conquer' or 'quantreg'.")
  }
}

#' @title Check that a parameter is a non-negative number
#' @param x The value to check.
#' @param name The parameter name used in the error message.
#' @return No return
check_nonnegative_number <- function(x, name)
{
  if (!is.numeric(x) || length(x) != 1 || x < 0) {
    stop(
      "Invalid input for '",
      name,
      "': ",
      x,
      ". Expected a positive number or zero.")
  }
}

#' @title Check that a parameter is a positive float
#' @param x The value to check.
#' @param name The parameter name used in the error message.
#' @return No return
check_positive_float <- function(x, name)
{
  if (!is.numeric(x) || length(x) != 1 || x <= 0) {
    stop(
      "Invalid input for '",
      name,
      "': ",
      x,
      ". Expected a positive float number.")
  }
}

#' @title Check that the seed is a positive integer or NULL
#' @param seed The value to check.
#' @return No return
check_seed <- function(seed)
{
  if (!(is.null(seed) || (is.numeric(seed) && length(seed) == 1 &&
                          seed %% 1 == 0 && seed > 0))) {
    stop("Invalid input for 'seed': ", seed, ". Expected a positive integer number or NULL.")
  }
}

# PER-METHODOLOGY PARAMETER CHECKERS -------------------------------------------

#' @title Check fqpca input parameters
#' @description Check input parameters of the fqpca methodology (also used by its CV drivers).
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
check_fqpca_params <- function(
    npc,
    quantile.value,
    periodic,
    splines.df,
    splines.method,
    penalized,
    lambda.ridge,
    tol,
    max.iters,
    verbose,
    seed)
{
  check_positive_integer(npc, 'npc')
  check_quantile_value(quantile.value)
  check_flag(periodic, 'periodic')

  # Check 'splines.df': integer positive number, larger than 'npc'
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
  check_flag(penalized, 'penalized')
  check_nonnegative_number(lambda.ridge, 'lambda.ridge')
  check_positive_float(tol, 'tol')
  check_iters(max.iters, 'max.iters')
  check_flag(verbose, 'verbose')
  check_seed(seed)
}

#' @title Check mfqpca input parameters
#' @description Check input parameters of the mfqpca methodology (also used by its CV driver).
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
mfqpca_check_params <- function(
    npc.between,
    npc.within,
    quantile.value,
    periodic,
    splines.df,
    splines.method,
    tol,
    max.iters,
    verbose,
    seed)
{
  check_positive_integer(npc.between, 'npc.between')
  check_positive_integer(npc.within, 'npc.within')
  check_quantile_value(quantile.value)
  check_flag(periodic, 'periodic')

  # Check 'splines.df': integer positive number, larger than both npc values
  if (!is.numeric(splines.df) || length(splines.df) != 1 ||
      splines.df %% 1 != 0 || splines.df < npc.between || splines.df < npc.within) {
    stop(
      "Invalid input for 'splines.df': ",
      splines.df,
      ". Expected an integer number larger or equal than 'npc.between' (",
      npc.between,
      ") and 'npc.within' (",
      npc.within,
      ").")
  }

  check_splines_method(splines.method)
  check_positive_float(tol, 'tol')
  check_iters(max.iters, 'max.iters')
  check_flag(verbose, 'verbose')
  check_seed(seed)
}

# INPUT DATA COERCION ----------------------------------------------------------

#' @title Coerce functional input data
#' @description Coerces the supported functional-data containers (matrix, tf object, or data.frame with a tf column) into an unnamed matrix, optionally resolving a grouping structure.
#' @param data An \eqn{(N \times T)} matrix, a tf object from the tidyfun package or a data.frame containing the functional data as a tf column.
#' @param colname The name of the column containing the functional data. Use only if data is a dataframe and colname is a column in the dataframe.
#' @param group Optional grouping structure of hierarchical data. If data is a dataframe this can be a name pointing to a column in the dataframe. If data is any other object it must be an array with length equal to the number of rows of data. Set to NULL when no grouping is used.
#' @return A list with elements \code{Y} (unnamed matrix of functional data) and \code{group} (the processed group vector, or NULL).
coerce_functional_input <- function(
    data,
    colname = NULL,
    group = NULL)
{
  # Helper function to check that the length of group matches a given number of rows
  check_group_length <- function(grp, n) {
    if (length(grp) != n) {
      stop(
        "Length of 'group' (",
        length(grp),
        ") does not match the number of rows (",
        n,
        ").")}}

  # Case 1: data is a matrix
  if (is.matrix(data)) {
    if (!is.null(group)) {
      # For matrix, a single string is not allowed for group.
      if (is.character(group) && length(group) == 1) {
        stop("For matrix input, 'group' cannot be a single string; it must be an array matching the number of rows.")
      }
      check_group_length(group, nrow(data))
    }
    return(list(Y = unname(data), group = group))
  }

  # Case 2: data is a 'tf' object
  if (tf::is_tf(data)) {
    mdata <- unname(as.matrix(data))
    if (!is.null(group)) {
      # For tf objects, group must be an array matching the rows, not a single string.
      if (is.character(group) && length(group) == 1) {
        stop("For 'tf' input, 'group' cannot be a single string; it must be an array matching the number of rows.")
      }
      check_group_length(group, nrow(mdata))
    }
    return(list(Y = mdata, group = group))
  }

  # Case 3: data is a data.frame
  if (is.data.frame(data)) {
    # Verify that colname is provided as a single string and it exists in data.
    if (is.null(colname) || !is.character(colname) || length(colname) != 1) {
      stop("Data is a data.frame, but 'colname' is not provided as a single string corresponding to a column name.")
    }
    if (!(colname %in% colnames(data))) {
      stop("Data is a data.frame, but 'colname' ('", colname, "') does not correspond to any column in the data.")
    }

    # Extract the functional data column and check its type.
    column_data <- data[[colname]]
    if (!tf::is_tf(column_data)) {
      stop("The column '", colname, "' in the data.frame is not a 'tf' object.")
    }
    Y <- unname(as.matrix(column_data))

    if (!is.null(group)) {
      if (is.character(group) && length(group) == 1) {
        # If group is a string, it must correspond to a column name in data.
        if (group %in% colnames(data)) {
          group <- data[[group]]
        } else {
          stop("For data.frame input, when 'group' is provided as a string ('", group, "'), it must match a column name in the data.")
        }
      }
      check_group_length(group, nrow(data))
    }
    return(list(Y = Y, group = group))
  }

  # If none of the supported types, throw an error.
  stop("Data is not a matrix, 'tf' object, or data.frame. Cannot process the input data.")
}

#' @title Checks input regressors
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
    stop("regressors cannot contain constant columns.")
  }
  names.regressors <- colnames(regressors)
  regressors <- unname(regressors)
  return(list(regressors=regressors, names.regressors=names.regressors))
}
