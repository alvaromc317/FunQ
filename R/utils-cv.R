
# GENERIC CROSS-VALIDATION FRAMEWORK -------------------------------------------

#' @title Generic cross-validation framework
#' @description Inner engine shared by all the cross-validation drivers: splits the data into folds, loops over a hyper-parameter grid fitting a model per (grid value, fold) pair via \code{fit_fun}, reconstructs the test predictions via \code{reconstruct_fun}, and collects the quantile error matrix.
#' @param Y The \eqn{(N \times T)} matrix of functional data.
#' @param grid Numeric vector of hyper-parameter values to cross-validate.
#' @param n.folds Number of folds to be used on cross validation.
#' @param criteria Criteria used to divide the data. Valid values are \code{'rows'} or \code{'points'}.
#' @param quantile.value The quantile considered.
#' @param fit_fun Function \code{(Y.train, grid.value, train.idx)} returning a fitted model. \code{train.idx} contains the training row indices when \code{criteria='rows'} and is NULL otherwise.
#' @param reconstruct_fun Function \code{(model, Y.test, train.idx)} returning the predicted \eqn{(N \times T)} matrix used to score the fold.
#' @param return.models Should the list of all the models built be returned?
#' @param verbose.cv Boolean indicating verbosity of the cross-validation process.
#' @param seed Seed for the random generator number (used to build the folds).
#' @param grid.label Label prefixing the grid value in verbose messages.
#' @param model.prefix Prefix used to name the stored models (\code{'<prefix>=<i>.fold=<j>'}).
#' @return A list containing the error matrix, the execution time, the grid, the criteria and the list of models.
cv_framework <- function(
    Y,
    grid,
    n.folds,
    criteria,
    quantile.value,
    fit_fun,
    reconstruct_fun,
    return.models = TRUE,
    verbose.cv = TRUE,
    seed = NULL,
    grid.label = 'Degrees of freedom: ',
    model.prefix = 'df_idx')
{
  start_time <- Sys.time()

  # KFOLDS
  Y.folds <- create_folds(
    Y,
    criteria = criteria,
    folds = n.folds,
    seed = seed)

  # Initialize error storage
  error.matrix <- matrix(-1, nrow = length(grid), ncol = n.folds)
  colnames(error.matrix) <- paste('Fold', seq_len(n.folds))
  list.models <- list()

  # CROSS VALIDATION
  for(i in seq_along(grid))
  {
    start_loop_time <- Sys.time()
    if(verbose.cv){message(grid.label, grid[i], ' ---------------------')}

    for(j in seq_len(n.folds))
    {
      if(verbose.cv){message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '. Fold: ', j)}

      Y.train <- Y.folds$Y.train.list[[j]]
      Y.test <- Y.folds$Y.test.list[[j]]
      train.idx <- if(criteria == 'rows'){Y.folds$train.idx.list[[j]]}else{NULL}

      # Execute model
      model <- fit_fun(Y.train, grid[i], train.idx)

      if(return.models)
      {
        name.model <- paste0(
          model.prefix,
          '=',
          i,
          '.fold=',
          j)
        list.models[[name.model]] <- model
      }

      # Obtain reconstruction and error
      Y.pred <- reconstruct_fun(model, Y.test, train.idx)
      error.matrix[i, j] <- quantile_error(Y = Y.test, Y.pred = Y.pred, quantile.value = quantile.value)
    }

    end_loop_time <-  Sys.time()
    if(verbose.cv){message(
      grid.label,
      grid[i],
      '. Execution completed in: ',
      round(difftime(end_loop_time, start_loop_time, units = "secs"), 3),
      ' seconds.')}
  }

  end_time <- Sys.time()
  execution.time <- difftime(end_time, start_time, units = "secs")
  return(list(
    error.matrix = error.matrix,
    execution.time = execution.time,
    grid = grid,
    criteria = criteria,
    list.models = list.models))
}
