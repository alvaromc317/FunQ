# MATRIX OPERATIONS -----------------------------------------------------------

#' @title Check if a matrix is singular
#' @description Given a matrix, return TRUE if singular
#' @param x Input matrix.
#' @return Boolean (true if)
is_rank_deficient <- function(x)
{
  # Compute the rank of the matrix
  qr_mat <- base::qr(x)
  # Check if the rank is less than the minimum dimension
  min_dim <- min(nrow(x), ncol(x))
  if(qr_mat$rank < min_dim){return(TRUE)}else{return(FALSE)}
}

#' @title Pivot long a matrix
#' @description Given a new matrix x of dimensions (n, m) pivot into a (n*m, 3) matrix indexing (row, column, value)
#' @param x Input matrix.
#' @return The pivoted matrix.
pivot_long_x <- function(x)
{
  rowid <- column <- NULL # Required to avoid warning when compiling package
  x.df <- data.frame(x)
  colnames(x.df) <- seq_len(ncol(x.df))
  x.df <- tibble::rowid_to_column(x.df)
  x.df <- tidyr::pivot_longer(data = x.df, cols = -rowid, names_to = 'column')
  x.df <- dplyr::mutate(.data = x.df, column = as.integer(column))
  x.df <- dplyr::rename(.data = x.df, row = rowid)
  x.df <- as.matrix(x.df)
  return(x.df)
}

#' @title Pivot short a matrix
#' @description Given a (n*m, 3) matrix indexing (row, column, value) pivots it into an (n, m) matrix.
#' @param x.long Input matrix.
#' @param dimensions array of two elements indicating the dimensions of the original data.
#' @return The pivoted matrix as a tibble
pivot_wide_x <- function(x.long, dimensions = NULL)
{
  x.long <- as.matrix(x.long)
  if(!is.null(dimensions))
  {
    dimensions <- c(max(x.long[, 1], dimensions[1]), max(x.long[, 2], dimensions[2]))
  } else {
    dimensions <- c(max(x.long[, 1]), max(x.long[, 2]))
  }
  x <- matrix(NA, nrow = dimensions[1], ncol = dimensions[2])
  x <- base::replace(x, x.long[,c(1,2)], x.long[, 3])
  return(x)
}

# QUANTILE BASED FUNCTIONS ----------------------------------------------------

#' @title Quantile regression objective value
#' @title Objective value function of a quantile regression model
#' @param quantile.value Default<-0.5. The quantile considered.
#' @param x The vector.
#' @return Objective value function for a quantile regression model.
quantile_function <- function(x, quantile.value = 0.5)
{
  # Quantile check function
  return(0.5 * base::abs(x) + (quantile.value - 0.5) * x)
}

#' @title Quantile error computation
#' @description Quantile error metric computation
#' @param Y \eqn{(N \times T)} matrix of observed time instants.
#' @param Y.pred \eqn{(N \times T)} matrix of predicted time instants.
#' @param quantile.value The quantile considered.
#' @return The quantile error between the two matrices
#' @export
quantile_error <- function(Y, Y.pred, quantile.value)
{
  mean(quantile_function(quantile.value = quantile.value, x = (Y - Y.pred)), na.rm = TRUE)
}

#' @title Proportion of points under quantile
#' @description Proportion of points under quantile.
#' @param Y \eqn{(N \times T)} matrix of observed time instants.
#' @param Y.pred \eqn{(N \times T)} matrix of predicted time instants.
#' @return The proportion of points under the estimated quantile function
#' @export
proportion_under_quantile <- function(Y, Y.pred)
{
  mean(Y < Y.pred, na.rm = T)
}

# TRAIN TEST SPLIT ------------------------------------------------------------

#' @title Split rows of a given matrix Y into train / test
#' @description Split rows of a given matrix Y into train / test based on parameters train.pct and train.size. If both are informed train.pct takes preference
#' @param Y \eqn{(N \times P)} matrix of predictive variables.
#' @param train.pct Float number indicating the % of rows used for training.
#' @param train.size Integer number indicating number of rows used for training. \code{train.size} is superseded by \code{train.pct}.
#' @param seed Seed for the random generator number.
#' @return A list containing a matrix Y.train and a matrix Y.test
train_test_split_rows <- function(Y, train.pct = NULL, train.size = NULL, seed = NULL)
{
  # train.pct takes preference over train.size
  if(is.null(train.pct) && is.null(train.size)){stop('Either train.pct or train.size must be filled.')}
  if(!is.null(seed)){set.seed(seed)}
  if(is.null(train.pct)){train.pct <- train.size / nrow(Y) }
  Y <- as.matrix(Y)
  Y.mask <- base::is.na(Y)
  curve_NAs <- base::rowSums(Y.mask) == ncol(Y)
  train_idx <- caret::createDataPartition(curve_NAs, p = train.pct, list = FALSE)
  Y.train <- Y[train_idx, ]
  Y.test <- Y[-train_idx,]
  results <- list(Y.train = Y.train, Y.test = Y.test)
  return(results)
}

#' @title Split points of a given matrix Y into train / test
#' @description Split points of a given matrix Y into train / test based on parameters train.pct and train.size. If both are informed train.pct takes preference. This function keeps the dimensions of the original Y in both train and test and substitutes the values in both splits by NAs
#' @param Y \eqn{(N \times P)} matrix of predictive variables.
#' @param train.pct Float number indicating the % of rows used for training.
#' @param train.size Integer number indicating number of rows used for training. \code{train.size} is superseded by \code{train.pct}.
#' @param seed Seed for the random generator number.
#' @return A list containing a matrix Y.train and a matrix Y.test
train_test_split_points <- function(Y, train.pct = NULL, train.size = NULL, seed = NULL)
{
  value <- NULL
  if(is.null(train.pct) && is.null(train.size)){stop('Either train.pct or train.size must be filled.')}
  if(!is.null(seed)){set.seed(seed)}
  if(is.null(train.pct)){train.pct <- train.size / length(Y) }
  Y.long <- pivot_long_x(Y)
  Y.long <- tibble::as_tibble(Y.long)
  Y.long <- dplyr::mutate(.data = Y.long, is_na = is.na(value))
  train_idx <- caret::createDataPartition(Y.long$is_na, p = train.pct, list = FALSE)
  train_idx <- as.vector(train_idx)
  Y.train <- Y.long[train_idx, c('row', 'column', 'value')]
  Y.train <- pivot_wide_x(x.long = Y.train, dimensions = c(nrow(Y), ncol(Y)))
  Y.test <- Y.long[-train_idx, c('row', 'column', 'value')]
  Y.test <- pivot_wide_x(x.long = Y.test, dimensions = c(nrow(Y), ncol(Y)))
  results <- list(Y.train = Y.train, Y.test = Y.test)
  return(results)
}

#' @title Split a given matrix Y into train / test
#' @description Splits a given matrix into a train / test split based on two possible criteria: either based on the total number of rows or the total number of data points.
#' @param Y \eqn{(N \times P)} matrix of predictive variables.
#' @param criteria Criteria used to divide the data. Valid values are \code{'rows'}, which considers the division based on full rows, or \code{'points'}, which considers the division based on points within the matrix.
#' @param train.pct Float number indicating the % of rows used for training. This takes precedence over \code{train.size}.
#' @param train.size Integer number indicating number of rows used for training. \code{train.size} is superseded by \code{train.pct}.
#' @param seed Seed for the random generator number.
#' @return A list containing a matrix Y.train and a matrix Y.test
#' @export
#' @examples
#' # Generate a small matrix
#' Y <- matrix(rnorm(50), nrow = 10)
#'
#' # Divide based on full rows
#' train_test_rows <- train_test_split(Y, criteria = 'rows', train.pct = 0.5)
#' Y.train <- train_test_rows$Y.train
#' Y.test <- train_test_rows$Y.test
#'
#' train_test_points <- train_test_split(Y, criteria = 'points', train.pct = 0.5)
#' Y.train <- train_test_points$Y.train
#' Y.test <- train_test_points$Y.test
train_test_split <- function(Y, criteria = 'points', train.pct = NULL, train.size = NULL, seed = NULL)
{
  Y <- base::unname(base::as.matrix(Y))
  if(criteria == 'points')
  {
    return(train_test_split_points(Y = Y, train.pct = train.pct, train.size = train.size, seed = seed))
  } else if(criteria == 'rows')
  {
    return(train_test_split_rows(Y = Y, train.pct = train.pct, train.size = train.size, seed = seed))
  } else{stop('Invalid criteria')}
}

# K-FOLD DIVISION -------------------------------------------------------------

#' @title Split a given matrix Y into k-folds
#' @description Splits full rows of a given matrix Y into k-folds
#' @param Y \eqn{(N \times P)} matrix of predictive variables.
#' @param folds Integer number indicating number of folds.
#' @param seed Seed for the random generator number.
#' @return A list containing two inside lists, one for training and one for testing. The length of the inside lists is equal to the number of folds
kfold_cv_rows <- function(Y, folds = 3, seed = NULL)
{
  if(!is.null(seed)){set.seed(seed)}
  Y <- as.matrix(Y)
  Y.mask <- base::is.na(Y)
  curve_NAs <- base::rowSums(Y.mask) == ncol(Y)
  test_idx <- caret::createFolds(curve_NAs, k = folds, list = FALSE)
  Y.train.list <- list()
  Y.test.list <- list()
  for(idx in 1:folds)
  {
    Y.train.list[[idx]] <- Y[test_idx != idx,]
    Y.test.list[[idx]] <- Y[test_idx == idx,]
  }
  results <- list(Y.train.list = Y.train.list, Y.test.list = Y.test.list)
  return(results)
}

#' @title Split a given matrix Y into k-folds
#' @description Splits a given matrix Y into k-folds. This function keeps the dimensions of the original Y in both train and test and substitutes the values in both split by NAs
#' @param Y \eqn{(N \times P)} matrix of predictive variables.
#' @param folds Integer number indicating number of folds.
#' @param seed Seed for the random generator number.
#' @return A list containing two inside lists, one for training and one for testing. The length of the inside lists is equal to the number of folds
kfold_cv_points <- function(Y, folds = 3, seed = NULL)
{
  value <- NULL # Required to avoid warning when compiling package
  if(!is.null(seed)){set.seed(seed)}
  Y.long <- pivot_long_x(Y)
  Y.long <- tibble::as_tibble(Y.long)
  Y.long <- dplyr::mutate(.data = Y.long, is_na = is.na(value))
  test_idx <- caret::createFolds(Y.long$is_na, k = folds, list = FALSE)
  Y.train.list <- list()
  Y.test.list <- list()
  for(idx in 1:folds)
  {
    Y.train.list[[idx]] <- Y.long[test_idx != idx, c('row', 'column', 'value')]
    Y.train.list[[idx]] <- pivot_wide_x(x.long = Y.train.list[[idx]], dimensions = c(nrow(Y), ncol(Y)))

    Y.test.list[[idx]] <- Y.long[test_idx == idx, c('row', 'column', 'value')]
    Y.test.list[[idx]] <- pivot_wide_x(x.long = Y.test.list[[idx]], dimensions = c(nrow(Y), ncol(Y)))
  }
  results <- list(Y.train.list = Y.train.list, Y.test.list = Y.test.list)
  return(results)
}

#' @title Split a given matrix Y into k-folds
#' @description Splits a given matrix Y into k-folds. based on two possible criteria: either based on the total number of rows or the total number of data points.
#' @param Y \eqn{(N \times P)} matrix of predictive variables.
#' @param criteria Criteria used to divide the data. Valid values are \code{'rows'}, which considers the division based on full rows, or \code{'points'}, which considers the division based on points within the matrix.
#' @param folds Integer number indicating number of folds
#' @param seed Seed for the random generator number.
#' @return A list containing two inside lists, one for training and one for testing. The length of the inside lists is equal to the number of folds
#' @export
#' @examples
#' # Generate a small matrix
#' Y <- matrix(rnorm(50), nrow = 10)
#'
#' # Divide based on full rows
#' kfolds_rows <- create_folds(Y, criteria = 'rows', folds = 3, seed = 1)
#' Y.train.list <- kfolds_rows$Y.train.list
#' Y.test.list <- kfolds_rows$Y.test.list
#'
#' kfolds_points <- create_folds(Y, criteria<-'points', folds<-3, seed<-1)
#' Y.train.list <- kfolds_rows$Y.train.list
#' Y.test.list <- kfolds_rows$Y.test.list
#'
create_folds <- function(Y, criteria = 'points', folds = 3, seed = NULL)
{
  Y <- base::unname(base::as.matrix(Y))
  if(criteria == 'points')
  {
    return(kfold_cv_points(Y = Y, folds = folds, seed = seed))
  } else if(criteria == 'rows')
  {
    return(kfold_cv_rows(Y = Y, folds = folds, seed = seed))
  } else{stop('Invalid criteria')}
}

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
#' @param method Method used in the resolution of the quantile regression model. It currently accepts the methods \code{c('br', 'fn', 'pfn', 'sfn')} from \code{quantreg} package along with any available solver in \code{CVXR} package.
#' @param penalized Boolean indicating if the smoothness should be controlled using a second derivative penalty. This functionality is experimental and is much slower than the control of the smoothness using the degrees of freedom.
#' @param lambda.grid An array containing the list of possible lambda values (these should be always positive numbers).
#' @param weights Optional array of weights used in penalized fqpca models. It must have dimension equal to (npc+1)
#' @param n.folds Number of folds to be used on cross validation.
#' @param return.models Should the list of all the models built be returned?
#' @param criteria Criteria used to divide the data. Valid values are \code{'rows'}, which considers the division based on full rows, or \code{'points'}, which considers the division based on points within the matrix.
#' @param tol Tolerance on the convergence of the algorithm.
#' @param n.iters Maximum number of iterations.
#' @param verbose.fqpca Boolean indicating verbosity of the fqpca function.
#' @param verbose.cv Boolean indicating verbosity of the cross-validation process.
#' @param seed Seed for the random generator number.
#' @param parallelized.scores Should the scores be computed in parallel? Experimantal component.
#' @param num.cores Number of cores to use in parallelized executions.
#' @return A list containing the matrix of scores, the matrix of loadings, a list with all the trained models (if the return_models param is TRUE) and a secondary list with extra information.
#' @export
#' @examples
#' # Generate fake dataset with 200 observations and 144 time points
#'
#' Y <- matrix(rep(sin(seq(0, 2*pi, length.out = 144)), 200), byrow = TRUE, nrow = 200)
#' Y <- Y + matrix(rnorm(200*144, 0, 0.4), nrow = 200)
#'
#' # Add missing observations
#' Y[sample(200*144, as.integer(0.2*200*144))] <- NA
#'
#' cv_result <- cross_validation_lambda(data=Y, lambda.grid = c(0, 1e-6), n.folds = 2)
cross_validation_lambda <- function(data, colname=NULL, npc = 2,  pve=NULL, quantile.value = 0.5, lambda.grid  =  c(0, 1e-10, 1e-5), n.folds = 3, return.models=TRUE, criteria = 'points', periodic = TRUE, splines.df = 10, tol = 1e-3, n.iters = 20, method = 'conquer', penalized=TRUE, weights=NULL, parallelized.scores=FALSE, num.cores=NULL, verbose.fqpca = FALSE, verbose.cv = TRUE, seed = NULL)
{
  start_time <- Sys.time()
  if(!base::is.null(seed))
  {
    base::set.seed(seed)
  } else{seed <- sample(1e9, 1)}

  if(method != 'conquer'){stop('For a penalized second derivative approach, the method must be conquer')}
  if(!penalized){stop('For a penalized second derivative approach, the penalized parameter must be set to TRUE.')}
  if(!n.folds == floor(n.folds)){stop('n.folds must be an integer number. Value provided: ', n.folds)}
  if(!(criteria %in% c('rows', 'points'))){stop('Invalid criteria. Valid criterias are c("rows", "points". Value provided: ', criteria)}

  # Check the input parameters except Y and colname
  check_fqpca_params(npc=npc, quantile.value=quantile.value, periodic=periodic,
                     splines.df=splines.df, method=method, penalized=penalized,
                     lambda.ridge=lambda.grid[1], tol=tol, n.iters=n.iters,
                     verbose=verbose.fqpca, seed=seed, parallelized.scores=parallelized.scores,
                     num.cores=num.cores)

  # Check Y and colname and return an unnamed matrix
  Y <- check_input_data(data=data, colname=colname)

  # KFOLDS
  Y.folds <- create_folds(Y, criteria = criteria, folds = n.folds, seed = seed)

  # Initialize error storage
  error.matrix <- matrix(-1, nrow = length(lambda.grid), ncol = n.folds)
  colnames(error.matrix) <- paste('Fold', 1:n.folds)
  list.models = list()

  # CROSS VALIDATION
  for(i in seq_along(lambda.grid))
  {
    start_loop_time <- Sys.time()
    if(verbose.cv){message('lambda=', lambda.grid[i], ' ---------------------')}
    lambda.ridge <- lambda.grid[i]
    for(j in 1:n.folds)
    {
      if(verbose.cv){message('Fold: ', j)}

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
      fqpca_results <- fqpca(data = Y.train, npc = npc,  quantile.value = quantile.value,  periodic = periodic, splines.df = splines.df, method = method, penalized=TRUE, lambda.ridge = lambda.ridge, weights=weights, tol = tol, n.iters = n.iters, verbose = verbose.fqpca, seed = seed)
      if(return.models)
      {
        name.model <- paste0('lambda_idx=', i, '.fold=', j)
        list.models[[name.model]] <- fqpca_results
      }
      if(is.null(pve))
      {
        npc.reconstruction <- fqpca_results$npc+1 # Add 1 to take the intercept into account
      } else{
        npc.reconstruction <- which(cumsum(fqpca_results$pve) >= pve)[1]+1
      }
      # Obtain reconstruction
      if(criteria == 'points')
      {
        loadings <- matrix(fqpca_results$loadings[, 1:npc.reconstruction], ncol=npc.reconstruction)
        scores <- matrix(fqpca_results$scores[, 1:npc.reconstruction], ncol=npc.reconstruction)
        Y.predicted <- scores %*% t(loadings)
      }else if(criteria=='rows')
      {
        test.scores <- predict.fqpca_object(fqpca_results, Y.test)
        loadings <- matrix(fqpca_results$loadings[, 1:npc.reconstruction], ncol=npc.reconstruction)
        scores <- matrix(test.scores[, 1:npc.reconstruction], ncol=npc.reconstruction)
        Y.predicted <- scores %*% t(loadings)
      }
      error.matrix[i, j] <- quantile_error(Y = Y.test, Y.pred = Y.predicted, quantile.value = quantile.value)
    }
    end_loop_time <-  Sys.time()
    if(verbose.cv){message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ' lambda: ', lambda.ridge, '. Execution completed in: ', round(difftime(end_loop_time, start_loop_time, units = "secs"), 3), ' seconds.')}
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
#' @param method Method used in the resolution of the quantile regression model. It currently accepts the methods \code{c('br', 'fn', 'pfn', 'sfn')} from \code{quantreg} package along with any available solver in \code{CVXR} package.
#' @param penalized Boolean indicating if the smoothness should be controlled using a second derivative penalty. This functionality is experimental and is much slower than the control of the smoothness using the degrees of freedom.
#' @param lambda.ridge  Hyper parameter controlling the penalization on the second derivative of the splines. It has effect only with \code{penalized=TRUE}. Experimantal component.
#' @param weights Optional array of weights used in penalized fqpca models. It must have dimension equal to (npc+1)
#' @param n.folds Number of folds to be used on cross validation.
#' @param return.models Should the list of all the models built be returned?
#' @param criteria Criteria used to divide the data. Valid values are \code{'rows'}, which considers the division based on full rows, or \code{'points'}, which considers the division based on points within the matrix.
#' @param tol Tolerance on the convergence of the algorithm.
#' @param n.iters Maximum number of iterations.
#' @param verbose.fqpca Boolean indicating verbosity of the fqpca function.
#' @param verbose.cv Boolean indicating verbosity of the cross-validation process.
#' @param seed Seed for the random generator number.
#' @param parallelized.scores Should the scores be computed in parallel? Experimantal component.
#' @param num.cores Number of cores to use in parallelized executions.
#' @return A list containing the matrix of scores, the matrix of loadings, and a secondary list with extra information.
#' @export
#' @examples
#' # Generate fake dataset with 200 observations and 144 time points
#'
#' Y <- matrix(rep(sin(seq(0, 2*pi, length.out = 144)), 200), byrow = TRUE, nrow = 200)
#' Y <- Y + matrix(rnorm(200*144, 0, 0.4), nrow = 200)
#'
#' # Add missing observations
#' Y[sample(200*144, as.integer(0.2*200*144))] <- NA
#'
#' cv_result <- cross_validation_df(data=Y, splines.df.grid = c(5, 10, 15), n.folds = 2)
cross_validation_df <- function(data, colname=NULL, npc = 2,  pve=NULL, quantile.value = 0.5,  lambda.ridge = 0, n.folds = 3, return.models = TRUE, criteria = 'points', periodic = TRUE, splines.df.grid = c(5, 10, 15, 20), tol = 1e-3, n.iters = 20, method = 'conquer', penalized=FALSE,  weights=NULL, parallelized.scores=FALSE, num.cores=NULL, verbose.fqpca = FALSE, verbose.cv = TRUE, seed = NULL)
{
  start_time <- Sys.time()
  if(!base::is.null(seed)){base::set.seed(seed)}

  if(!n.folds == floor(n.folds)){stop('n.folds must be an integer number. Value provided: ', n.folds)}
  if(!(criteria %in% c('rows', 'points'))){stop('Invalid criteria. Valid criterias are c("rows", "points". Value provided: ', criteria)}

  # Check the input parameters except Y and colname
  check_fqpca_params(npc=npc, quantile.value=quantile.value, periodic=periodic,
                     splines.df=splines.df.grid[1], method=method, penalized=penalized,
                     lambda.ridge=lambda.ridge, tol=tol, n.iters=n.iters,
                     verbose=verbose.fqpca, seed=seed, parallelized.scores=parallelized.scores,
                     num.cores=num.cores)

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

    # Detect if the number of components is larger than the degrees of freedom.
    npc_splines_warning <- FALSE
    if(npc > splines.df)
    {
      npc_splines_warning <- TRUE
      warning('The npc is larger than splines.df\nCurrent splines.df value: ', splines.df, '\nTaking npc=splines.df for this iteration of the Cross-Validation process.')
      npc_original <- npc
      npc <- splines.df
    }

    for(j in 1:n.folds)
    {
      if(verbose.cv){message('Fold: ', j)}

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
      fqpca_results <- fqpca(data = Y.train, npc = npc,  quantile.value = quantile.value,  periodic = periodic, splines.df = splines.df, method = method, penalized=penalized, lambda.ridge = lambda.ridge, weights = weights, tol = tol, n.iters = n.iters, verbose = verbose.fqpca, seed = seed)
      if(return.models)
      {
        name.model <- paste0('df_idx=', i, '.fold=', j)
        list.models[[name.model]] <- fqpca_results
      }
      if(is.null(pve))
      {
        npc.reconstruction <- fqpca_results$npc+1 # Add 1 to take the intercept into account
      } else{
        npc.reconstruction <- which(cumsum(fqpca_results$pve) >= pve)[1]+1
      }
      # Obtain reconstruction
      if(criteria == 'points')
      {
        loadings <- matrix(fqpca_results$loadings[, 1:npc.reconstruction], ncol=npc.reconstruction)
        scores <- matrix(fqpca_results$scores[, 1:npc.reconstruction], ncol=npc.reconstruction)
        Y.predicted <- scores %*% t(loadings)
      }else if(criteria=='rows')
      {
        test.scores <- predict.fqpca_object(fqpca_results, Y.test)
        loadings <- matrix(fqpca_results$loadings[, 1:npc.reconstruction], ncol=npc.reconstruction)
        scores <- matrix(test.scores[, 1:npc.reconstruction], ncol=npc.reconstruction)
        Y.predicted <- scores %*% t(loadings)
      }
      error.matrix[i, j] <- quantile_error(Y = Y.test, Y.pred = Y.predicted, quantile.value = quantile.value)
    }

    end_loop_time <-  Sys.time()
    if(verbose.cv){message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ' Degrees of freedom: ', splines.df, '. Execution completed in: ', round(difftime(end_loop_time, start_loop_time, units = "secs"), 3), ' seconds.')}
  }

  # Recover the original npc value
  if(npc_splines_warning)
  {
    npc <- npc_original
  }

  end_time <- Sys.time()
  execution.time <- difftime(end_time, start_time, units = "secs")
  return(list(error.matrix = error.matrix, execution.time = execution.time, splines.df.grid = splines.df.grid, criteria = criteria, list.models = list.models))
}
