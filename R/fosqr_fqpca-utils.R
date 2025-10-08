
#' @title CROSS VALIDATION OF DEGREES OF FREEDOM
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
#' @param verbose.fosqr_fqpca Boolean indicating verbosity of the fqpca function.
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
  start_time <- Sys.time()
  if(!base::is.null(seed)){base::set.seed(seed)}

  if(n.folds != floor(n.folds)){stop('n.folds must be an integer number. Value provided: ', n.folds)}
  if(!(criteria %in% c('rows', 'points'))){stop('Invalid criteria. Valid criterias are c("rows", "points". Value provided: ', criteria)}
  if(!all(splines.df.grid == floor(splines.df.grid))){stop('splines.df.grid must be a positive integer array.')}

  # Check the input parameters except Y
  check_fosqr_fqpca_params(npc=npc, quantile.value=quantile.value, periodic=periodic,
                           splines.df=max(npc, splines.df.grid), splines.method=splines.method,
                           tol=tol, max.iters=max.iters,verbose=verbose.cv, seed=seed)

  # Check Y and regressors
  Y <- check_Y(Y)
  regressors_checked <- check_regressors(regressors)
  regressors <- regressors_checked$regressors

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

      Y.train <- Y.folds$Y.train.list[[j]]
      Y.test <- Y.folds$Y.test.list[[j]]
      if(criteria == 'rows')
      {
        regressors.train <- regressors[Y.folds$train.idx.list[[j]], ]
        regressors.test <- regressors[-Y.folds$train.idx.list[[j]], ]
      }else{
        regressors.train <- regressors
        regressors.test <- regressors
      }

      # Execute model
      fosqr_fqpca_results <- fosqr_fqpca(
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

      if(return.models)
      {
        name.model <- paste0('df_idx=', i, '.fold=', j)
        list.models[[name.model]] <- fosqr_fqpca_results
      }

      # Obtain optimal number of PC
      npc.reconstruction <- obtain_npc(scores=fosqr_fqpca_results$fqpca.scores, pve=pve)

      # Obtain reconstruction
      if(criteria == 'points')
      {
        scores <- fosqr_fqpca_results$fqpca.scores[, 1:npc.reconstruction, drop=FALSE]
      }else if(criteria=='rows')
      {
        test.scores <- predict.fosqr_fqpca_object(fosqr_fqpca_results, newdata=Y.test, newregressors=regressors.test)
        scores <- test.scores[, 1:npc.reconstruction, drop=FALSE]
      }

      # Build fosqr fitted value
      fosqr.fitted <- sweep(regressors.test %*% t(fosqr_fqpca_results$fosqr.loadings), MARGIN = 2, STATS = fosqr_fqpca_results$fosqr.intercept, FUN = "+")

      # Build FQPCA fitted value
      if (npc.reconstruction > 0){
        fqpca.fitted <- scores[, seq_len(npc.reconstruction), drop = FALSE] %*% t(fosqr_fqpca_results$fqpca.loadings[, seq_len(npc.reconstruction), drop = FALSE])
      }else{
        fqpca.fitted <- matrix(0, nrow = nrow(scores), ncol = nrow(fosqr_fqpca_results$fqpca.loadings))
      }
      Y.pred <- fosqr.fitted + fqpca.fitted
      error.matrix[i, j] <- quantile_error(Y = Y.test, Y.pred = Y.pred, quantile.value = quantile.value)
    }

    end_loop_time <-  Sys.time()
    if(verbose.cv){message('Degrees of freedom: ', splines.df, '. Execution completed in: ', round(difftime(end_loop_time, start_loop_time, units = "secs"), 3), ' seconds.')}
  }

  end_time <- Sys.time()
  execution.time <- difftime(end_time, start_time, units = "secs")
  return(list(error.matrix = error.matrix, execution.time = execution.time, splines.df.grid = splines.df.grid, criteria = criteria, list.models = list.models))
}
