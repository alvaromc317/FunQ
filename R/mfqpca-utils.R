
#' @title CROSS VALIDATION OF DEGREES OF FREEDOM
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
    splines.df.grid = c(5, 10, 15, 20), tol = 1e-3, max.iters = 20,
    splines.method = 'conquer', verbose.mfqpca = FALSE, verbose.cv = TRUE, seed = NULL)
{
  start_time <- Sys.time()
  if(!base::is.null(seed)){base::set.seed(seed)}

  if(!(n.folds == floor(n.folds))){stop('n.folds must be an integer number. Value provided: ', n.folds)}
  if(!(criteria %in% c('rows', 'points'))){stop('Invalid criteria. Valid criterias are c("rows", "points". Value provided: ', criteria)}

  # Check the input parameters except Y
  mfqpca_check_params(npc.between=npc.between, npc.within=npc.within, quantile.value=quantile.value, periodic=periodic,
                      splines.df=max(npc.within, npc.between, splines.df.grid), splines.method=splines.method, tol=tol, max.iters=max.iters,
                      verbose=verbose.mfqpca, seed=seed)

  # Check Y and colname and return an unnamed matrix
  r <- mfqpca_check_input_data(data=data, colname=colname, group=group)
  Y <- r[[1]]
  group <- r[[2]]

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

    if(npc.between > splines.df)warning('The npc.between is larger than splines.df\nCurrent splines.df value: ', splines.df, '\nTaking npc.between=splines.df for this iteration of the Cross-Validation process.')
    if(npc.within > splines.df)warning('The npc.within is larger than splines.df\nCurrent splines.df value: ', splines.df, '\nTaking npc.within=splines.df for this iteration of the Cross-Validation process.')

    npc.between.iteration <- min(splines.df, npc.between)
    npc.within.iteration <- min(splines.df, npc.within)

    for(j in 1:n.folds)
    {
      if(verbose.cv){message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '. Fold: ', j)}

      # Access data depending on criteria
      if(criteria=='rows')
      {
        Y.train <- Y.folds$Y.train.list[[j]]
        Y.test <- Y.folds$Y.test.list[[j]]
        group.train = group[Y.folds$train.idx.list[[j]]]
        group.test = group[-Y.folds$train.idx.list[[j]]]
      } else if(criteria == 'points')
      {
        Y.train <- Y.folds$Y.train.list[[j]]
        Y.test <-  Y.folds$Y.test.list[[j]]
        group.train = group
      } else{stop('Invalid value for criteria. Valid values are observations or curves.')}

      # Execute model
      mfqpca_results <- mfqpca(
        data = Y.train, group = group.train, colname = colname, npc.between = npc.between.iteration,
        npc.within = npc.within.iteration, quantile.value = quantile.value,  periodic = periodic,
        splines.df = splines.df, splines.method = splines.method, tol = tol,
        max.iters = max.iters, verbose = verbose.mfqpca, seed = seed)

      if(return.models)
      {
        name.model <- paste0('df_idx=', i, '.fold=', j)
        list.models[[name.model]] <- mfqpca_results
      }

      # Obtain optimal number of PC
      n.components.between <- obtain_npc(scores=mfqpca_results$scores.between, pve=pve.between)
      n.components.within <- obtain_npc(scores=mfqpca_results$scores.within, pve=pve.within)

      if(criteria == 'points')
      {
        scores.between.full <- mfqpca_results$scores.between.full
        scores.within <- mfqpca_results$scores.within
      }else if(criteria=='rows')
      {
        test.scores <- predict.mfqpca_object(mfqpca_results, newdata.group=group.test, Y.test)
        scores.between.full <- test.scores$scores.between.full
        scores.within <- test.scores$scores.within
      }
      intercept <- mfqpca_results$intercept
      loadings.between <- mfqpca_results$loadings.between
      loadings.within <- mfqpca_results$loadings.within

      Y.pred =
        mfqpca_reconstruct(scores.between.full, loadings.between, n.components.between) +
        mfqpca_reconstruct(scores.within, loadings.within, n.components.within)

      Y.pred = sweep(Y.pred, MARGIN = 2, STATS = intercept, FUN = "+")

      error.matrix[i, j] <- quantile_error(Y = Y.test, Y.pred = Y.pred, quantile.value = quantile.value)
    }

    end_loop_time <-  Sys.time()
    if(verbose.cv){message('Degrees of freedom: ', splines.df, '. Execution completed in: ', round(difftime(end_loop_time, start_loop_time, units = "secs"), 3), ' seconds.')}
  }

  end_time <- Sys.time()
  execution.time <- difftime(end_time, start_time, units = "secs")
  return(list(error.matrix = error.matrix, execution.time = execution.time, splines.df.grid = splines.df.grid, criteria = criteria, list.models = list.models))
}
