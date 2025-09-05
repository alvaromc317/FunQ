
# CROSS VALIDATION ------------------------------------------------------------

#' @title CROSS VALIDATION OF LAMBDA.RIDGE
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
  start_time <- Sys.time()
  if(!base::is.null(seed))
  {
    base::set.seed(seed)
  } else{seed <- sample(1e9, 1)}

  if(splines.method != 'conquer'){stop('For a penalized second derivative approach, the splines.method must be conquer')}
  if(!penalized){stop('For a penalized second derivative approach, the penalized parameter must be set to TRUE.')}
  if(!n.folds == floor(n.folds)){stop('n.folds must be an integer number. Value provided: ', n.folds)}
  if(!(criteria %in% c('rows', 'points'))){stop('Invalid criteria. Valid criterias are c("rows", "points". Value provided: ', criteria)}

  # Check the input parameters except Y and colname
  check_fqpca_params(npc=npc, quantile.value=quantile.value, periodic=periodic,
                     splines.df=splines.df, splines.method=splines.method, penalized=penalized,
                     lambda.ridge=lambda.grid[1], tol=tol, max.iters=max.iters,
                     verbose=verbose.fqpca, seed=seed)

  # Check Y and colname and return an unnamed matrix
  Y <- check_input_data(data=data, colname=colname)

  # KFOLDS
  Y.folds <- create_folds(Y, criteria = criteria, folds = n.folds, seed = seed)

  # Initialize error storage
  error.matrix <- matrix(-1, nrow = length(lambda.grid), ncol = n.folds)
  colnames(error.matrix) <- paste('Fold', 1:n.folds)
  list.models <- list()

  # CROSS VALIDATION
  for(i in seq_along(lambda.grid))
  {
    start_loop_time <- Sys.time()
    if(verbose.cv){message('lambda=', lambda.grid[i], ' ---------------------')}
    lambda.ridge <- lambda.grid[i]
    for(j in 1:n.folds)
    {
      if(verbose.cv){message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '. Fold: ', j)}

      # Access data depending on criteria
      if(criteria=='rows')
      {
        Y.train <- Y.folds$Y.train.list[[j]]
        Y.test <- Y.folds$Y.test.list[[j]]
      } else if(criteria == 'points')
      {
        Y.train <- Y.folds$Y.train.list[[j]]
        Y.test <-  Y.folds$Y.test.list[[j]]
      } else{stop('Invalid value for criteria. Valid values are observations or curves')}

      # Execute model
      fqpca_results <- fqpca(
        data = Y.train, npc = npc,  quantile.value = quantile.value,
        periodic = periodic, splines.df = splines.df, splines.method = splines.method,
        penalized = TRUE, lambda.ridge = lambda.ridge, tol = tol, max.iters = max.iters,
        verbose = verbose.fqpca, seed = seed)

      if(return.models)
      {
        name.model <- paste0('lambda_idx=', i, '.fold=', j)
        list.models[[name.model]] <- fqpca_results
      }

      npc.reconstruction <- obtain_npc(scores=fqpca_results$scores, pve=pve)

      # Obtain reconstruction
      if(criteria == 'points')
      {
        scores <- fqpca_results$scores[, 1:npc.reconstruction, drop=F]
      }else if(criteria=='rows')
      {
        test.scores <- predict.fqpca_object(fqpca_results, Y.test)
        scores <- test.scores[, 1:npc.reconstruction, drop=F]
      }
      loadings <- fqpca_results$loadings[, 1:npc.reconstruction, drop=F]
      intercept <- fqpca_results$intercept
      Y.pred <- sweep(scores %*% t(loadings), MARGIN = 2, STATS = intercept, FUN = "+")
      error.matrix[i, j] <- quantile_error(Y = Y.test, Y.pred = Y.pred, quantile.value = quantile.value)
    }
    end_loop_time <-  Sys.time()
    if(verbose.cv){message('lambda: ', lambda.ridge, '. Execution completed in: ', round(difftime(end_loop_time, start_loop_time, units = "secs"), 3), ' seconds.')}
  }
  end_time <- Sys.time()
  execution.time <- difftime(end_time, start_time, units = "secs")
  return(list(error.matrix = error.matrix, execution.time = execution.time, lambda.grid = lambda.grid, criteria = criteria, list.models = list.models))
}


#' @title CROSS VALIDATION OF DEGREES OF FREEDOM
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
    periodic = TRUE, splines.df.grid = c(5, 10, 15, 20), tol = 1e-3,
    max.iters = 20, splines.method = 'conquer', penalized=FALSE, verbose.fqpca = FALSE,
    verbose.cv = TRUE, seed = NULL)
{
  start_time <- Sys.time()
  if(!base::is.null(seed)){base::set.seed(seed)}

  if(!n.folds == floor(n.folds)){stop('n.folds must be an integer number. Value provided: ', n.folds)}
  if(!(criteria %in% c('rows', 'points'))){stop('Invalid criteria. Valid criterias are c("rows", "points". Value provided: ', criteria)}

  # Check the input parameters except Y and colname
  check_fqpca_params(npc=npc, quantile.value=quantile.value, periodic=periodic,
                     splines.df=max(npc, splines.df.grid), splines.method=splines.method, penalized=penalized,
                     lambda.ridge=lambda.ridge, tol=tol, max.iters=max.iters,
                     verbose=verbose.fqpca, seed=seed)

  # Check Y and colname and return an unnamed matrix
  Y <- check_input_data(data=data, colname=colname)

  # KFOLDS
  Y.folds <- create_folds(Y, criteria = criteria, folds = n.folds, seed = seed)

  # Initialize error storage
  error.matrix <- matrix(-1, nrow = length(splines.df.grid), ncol = n.folds)
  colnames(error.matrix) <- paste('Fold', 1:n.folds)
  list.models <- list()
  # CROSS VALIDATION
  for(i in seq_along(splines.df.grid))
  {
    start_loop_time <- Sys.time()
    if(verbose.cv){message('Degrees of freedom: ', splines.df.grid[i], ' ---------------------')}
    splines.df <- splines.df.grid[i]

    if(npc > splines.df)warning('The npc is larger than splines.df\nCurrent splines.df value: ', splines.df, '\nTaking npc=splines.df for this iteration of the Cross-Validation process.')

    npc.iteration <- min(splines.df, npc)

    for(j in 1:n.folds)
    {
      if(verbose.cv){message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '. Fold: ', j)}

      # Access data depending on criteria
      if(criteria=='rows')
      {
        Y.train <- Y.folds$Y.train.list[[j]]
        Y.test <- Y.folds$Y.test.list[[j]]
      } else if(criteria == 'points')
      {
        Y.train <- Y.folds$Y.train.list[[j]]
        Y.test <-  Y.folds$Y.test.list[[j]]
      } else{stop('Invalid value for criteria. Valid values are observations or curves.')}

      # Execute model
      fqpca_results <- fqpca(
        data = Y.train, npc = npc.iteration,  quantile.value = quantile.value,
        periodic = periodic, splines.df = splines.df, splines.method = splines.method,
        penalized=penalized, lambda.ridge = lambda.ridge, tol = tol,
        max.iters = max.iters, verbose = verbose.fqpca, seed = seed)

      if(return.models)
      {
        name.model <- paste0('df_idx=', i, '.fold=', j)
        list.models[[name.model]] <- fqpca_results
      }

      npc.reconstruction <- obtain_npc(scores=fqpca_results$scores, pve=pve)

      # Obtain reconstruction
      if(criteria == 'points')
      {
        scores <- fqpca_results$scores[, 1:npc.reconstruction, drop=F]
      }else if(criteria=='rows')
      {
        test.scores <- predict.fqpca_object(fqpca_results, Y.test)
        scores <- test.scores[, 1:npc.reconstruction, drop=F]
      }
      loadings <- fqpca_results$loadings[, 1:npc.reconstruction, drop=F]
      intercept <- fqpca_results$intercept
      Y.pred <- sweep(scores %*% t(loadings), MARGIN = 2, STATS = intercept, FUN = "+")
      error.matrix[i, j] <- quantile_error(Y = Y.test, Y.pred = Y.pred, quantile.value = quantile.value)
    }

    end_loop_time <-  Sys.time()
    if(verbose.cv){message('Degrees of freedom: ', splines.df, '. Execution completed in: ', round(difftime(end_loop_time, start_loop_time, units = "secs"), 3), ' seconds.')}
  }

  end_time <- Sys.time()
  execution.time <- difftime(end_time, start_time, units = "secs")
  return(list(error.matrix = error.matrix, execution.time = execution.time, splines.df.grid = splines.df.grid, criteria = criteria, list.models = list.models))
}
