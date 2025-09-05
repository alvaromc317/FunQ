
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
  mean(Y <= Y.pred, na.rm = T)
}

# VARIABILITY OF SCORES -------------------------------------------------------

#' @title Explained variability of the fqpca scores computation
#' @description Computes the percentage of explained variability based on the variance of the scores matrix
#' @param scores Matrix of scores.
#' @return The percentage of variability each component is explaining.
compute_explained_variability <- function(scores)
{
  variability <- base::diag(stats::var(scores))
  pve <- variability / base::sum(variability)
  return(pve)
}

#' @title Obtain npc
#' @description Computes the required number of components to achieve the desired pve value
#' @param scores matrix of scores
#' @param pve If smaller than 1, taken as percentage of explained variability. If greater than 1, taken as number of components.
obtain_npc <- function(scores, pve)
{
  if(is.null(pve))
  {
    npc <- ncol(scores)
  }else if(pve < 1)
  {
    score.variance <- cumsum(compute_explained_variability(scores))
    npc <- min(which(score.variance > pve))
  }else
  {
    if(!(pve == floor(pve))){stop('pve must be either a floating point number smaller than 1 or an integer number smaller than the number of pc. Value provided: ', pve)}
    if(pve > ncol(scores)){stop('pve must be either a floating point number smaller than 1 or an integer number smaller than the number of pc. Value provided: ', pve)}
    npc <- pve
  }
  return(npc)
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
