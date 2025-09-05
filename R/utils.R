
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
  if(is.null(pve)){
    npc <- ncol(scores)
  }else if(pve == 0){
    npc <- 0
  }else if(pve < 1){
    score.variance <- cumsum(compute_explained_variability(scores))
    npc <- min(which(score.variance >= pve))
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
  # Validate inputs
  if (is.null(train.pct) && is.null(train.size)) {
    stop("Either train.pct or train.size must be provided.")
  }
  if (!is.null(train.pct)) {
    if (!is.numeric(train.pct) || length(train.pct) != 1L || is.na(train.pct) || train.pct < 0 || train.pct > 1) {
      stop("train.pct must be a numeric scalar in [0, 1].")
    }
  }
  if (!is.null(train.size)) {
    if (!is.numeric(train.size) || length(train.size) != 1L || is.na(train.size)) {
      stop("train.size must be a numeric scalar.")
    }
  }
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) stop("seed must be a numeric scalar.")
    set.seed(seed)
  }

  Y <- as.matrix(Y)
  n <- nrow(Y); m <- ncol(Y)
  if (is.null(n) || is.null(m)) stop("Y must be a 2D object (matrix/data.frame) coercible to a matrix.")
  if (n == 0L) stop("Y must have at least one row.")

  # Determine train proportion
  if (is.null(train.pct)) {
    if (train.size < 0 || train.size > n) stop("train.size must be in [0, nrow(Y)].")
    p <- if (n == 0) 0 else train.size / n
  } else {
    p <- train.pct
  }

  # Stratify rows: all-NA vs not-all-NA to mirror original intent
  all_na <- rowSums(is.na(Y)) == m
  idx_true  <- which(all_na)
  idx_false <- which(!all_na)

  total_target <- as.integer(round(p * n))
  ns <- c(length(idx_true), length(idx_false))
  raw <- p * ns
  base_alloc <- floor(raw)
  frac <- raw - base_alloc

  alloc <- base_alloc
  remainder <- total_target - sum(alloc)
  if (remainder > 0L) {
    order_add <- order(frac, decreasing = TRUE)
    for (j in order_add) {
      if (remainder <= 0L) break
      if (alloc[j] < ns[j]) {
        alloc[j] <- alloc[j] + 1L
        remainder <- remainder - 1L
      }
    }
  } else if (remainder < 0L) {
    order_sub <- order(frac, decreasing = FALSE)
    for (j in order_sub) {
      if (remainder >= 0L) break
      if (alloc[j] > 0L) {
        alloc[j] <- alloc[j] - 1L
        remainder <- remainder + 1L
      }
    }
  }
  alloc <- pmin(pmax(alloc, 0L), ns)

  pick_from <- function(ix, k) {
    if (k <= 0L || length(ix) == 0L) return(integer(0L))
    if (k >= length(ix)) return(ix)
    sample(ix, size = k, replace = FALSE)
  }

  train_idx <- c(
    pick_from(idx_true,  alloc[1L]),
    pick_from(idx_false, alloc[2L])
  )
  train_idx <- sort(train_idx)
  test_idx <- setdiff(seq_len(n), train_idx)

  Y.train <- Y[train_idx, , drop = FALSE]
  Y.test  <- Y[test_idx,  , drop = FALSE]

  list(Y.train = Y.train, Y.test = Y.test, train.idx=train_idx)
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
  # Validate inputs
  if (is.null(train.pct) && is.null(train.size)) {
    stop("Either train.pct or train.size must be provided.")
  }
  if (!is.null(train.pct)) {
    if (!is.numeric(train.pct) || length(train.pct) != 1L || is.na(train.pct) || train.pct < 0 || train.pct > 1) {
      stop("train.pct must be a numeric scalar in [0, 1].")
    }
  }
  if (!is.null(train.size)) {
    if (!is.numeric(train.size) || length(train.size) != 1L || is.na(train.size) || train.size < 0) {
      stop("train.size must be a non-negative numeric scalar.")
    }
  }
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) stop("seed must be a numeric scalar.")
    set.seed(seed)
  }

  Y <- as.matrix(Y)
  n <- nrow(Y); m <- ncol(Y)
  if (is.null(n) || is.null(m)) stop("Y must be a 2D object (matrix/data.frame) coercible to a matrix.")
  if (n == 0L || m == 0L) stop("Y must have at least one row and one column.")

  total_points <- n * m

  # Determine train proportion
  if (is.null(train.pct)) {
    if (train.size > total_points) stop("train.size cannot exceed length(Y) = nrow(Y) * ncol(Y).")
    p <- if (total_points == 0) 0 else train.size / total_points
  } else {
    p <- train.pct
  }

  # Per-row allocation to keep the proportion within each row
  raw_row <- rep(p * m, n)
  base_k <- floor(raw_row)
  frac <- raw_row - base_k
  target_total <- as.integer(round(p * total_points))
  remainder <- target_total - sum(base_k)

  k_row <- base_k
  if (remainder > 0L) {
    order_add <- order(frac, decreasing = TRUE)
    for (i in order_add) {
      if (remainder <= 0L) break
      if (k_row[i] < m) {
        k_row[i] <- k_row[i] + 1L
        remainder <- remainder - 1L
      }
    }
  } else if (remainder < 0L) {
    order_sub <- order(frac, decreasing = FALSE)
    for (i in order_sub) {
      if (remainder >= 0L) break
      if (k_row[i] > 0L) {
        k_row[i] <- k_row[i] - 1L
        remainder <- remainder + 1L
      }
    }
  }
  k_row <- pmin(pmax(k_row, 0L), m)

  # Prepare outputs
  Y.train <- matrix(NA, nrow = n, ncol = m)
  Y.test  <- matrix(NA, nrow = n, ncol = m)

  cols <- seq_len(m)
  for (i in seq_len(n)) {
    ki <- k_row[i]
    if (ki <= 0L) {
      Y.test[i, ] <- Y[i, ]
    } else if (ki >= m) {
      Y.train[i, ] <- Y[i, ]
    } else {
      trc <- sample(cols, size = ki, replace = FALSE)
      tec <- cols[!(cols %in% trc)]
      if (length(trc)) Y.train[i, trc] <- Y[i, trc]
      if (length(tec)) Y.test[i, tec] <- Y[i, tec]
    }
  }

  list(Y.train = Y.train, Y.test = Y.test)
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
  # Validate inputs
  if (!is.numeric(folds) || length(folds) != 1L || is.na(folds)) stop("folds must be a numeric scalar.")
  folds <- as.integer(folds)
  if (folds < 2L) stop("folds must be at least 2.")
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) stop("seed must be a numeric scalar.")
    set.seed(seed)
  }

  Y <- as.matrix(Y)
  n <- nrow(Y); m <- ncol(Y)
  if (is.null(n) || is.null(m)) stop("Y must be a 2D object (matrix/data.frame) coercible to a matrix.")
  if (n == 0L) stop("Y must have at least one row.")
  if (folds > n) stop("folds cannot exceed the number of rows in Y.")

  all_na <- rowSums(is.na(Y)) == m
  idx_true  <- which(all_na)
  idx_false <- which(!all_na)

  # Assign folds within each stratum; always return a data.frame (even if empty)
  assign_folds <- function(ix, k) {
    if (length(ix) == 0L) return(data.frame(index = integer(0), fold = integer(0)))
    perm <- sample(ix, length(ix), replace = FALSE)
    fold_seq <- rep_len(seq_len(k), length(ix))
    data.frame(index = perm, fold = fold_seq)
  }

  df_true  <- assign_folds(idx_true,  folds)
  df_false <- assign_folds(idx_false, folds)

  fold_id <- integer(n)
  if (nrow(df_true) > 0L)  fold_id[df_true$index]  <- df_true$fold
  if (nrow(df_false) > 0L) fold_id[df_false$index] <- df_false$fold

  Y.train.list <- vector("list", folds)
  Y.test.list  <- vector("list", folds)
  train.idx.list <- vector("list", folds)
  for (k in seq_len(folds)) {
    test_rows  <- which(fold_id == k)
    train_rows <- which(fold_id != k)
    Y.train.list[[k]] <- Y[train_rows, , drop = FALSE]
    Y.test.list[[k]]  <- Y[test_rows,  , drop = FALSE]
    train.idx.list[[k]] <- train_rows
  }

  list(Y.train.list = Y.train.list, Y.test.list = Y.test.list)
}

#' @title Split a given matrix Y into k-folds
#' @description Splits a given matrix Y into k-folds. This function keeps the dimensions of the original Y in both train and test and substitutes the values in both split by NAs
#' @param Y \eqn{(N \times P)} matrix of predictive variables.
#' @param folds Integer number indicating number of folds.
#' @param seed Seed for the random generator number.
#' @return A list containing two inside lists, one for training and one for testing. The length of the inside lists is equal to the number of folds
kfold_cv_points <- function(Y, folds = 3, seed = NULL)
{
  # Validate inputs
  if (!is.numeric(folds) || length(folds) != 1L || is.na(folds)) stop("folds must be a numeric scalar.")
  folds <- as.integer(folds)
  if (folds < 2L) stop("folds must be at least 2.")
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) stop("seed must be a numeric scalar.")
    set.seed(seed)
  }

  Y <- as.matrix(Y)
  n <- nrow(Y); m <- ncol(Y)
  if (is.null(n) || is.null(m)) stop("Y must be a 2D object (matrix/data.frame) coercible to a matrix.")
  if (n == 0L || m == 0L) stop("Y must have at least one row and one column.")

  # Per-row fold assignment for columns: balanced and randomized
  fold_ids <- matrix(0L, nrow = n, ncol = m)
  base_pattern <- seq_len(folds)
  for (i in seq_len(n)) {
    perm <- sample.int(m, m, replace = FALSE)
    fold_seq <- rep_len(base_pattern, m)
    fold_ids[i, perm] <- fold_seq
  }

  Y.train.list <- vector("list", folds)
  Y.test.list  <- vector("list", folds)
  for (k in seq_len(folds)) {
    Ytrain <- matrix(NA, nrow = n, ncol = m)
    Ytest  <- matrix(NA, nrow = n, ncol = m)
    test_mask <- (fold_ids == k)
    train_mask <- !test_mask
    if (any(train_mask)) Ytrain[train_mask] <- Y[train_mask]
    if (any(test_mask))  Ytest[test_mask]   <- Y[test_mask]
    Y.train.list[[k]] <- Ytrain
    Y.test.list[[k]]  <- Ytest
  }

  list(Y.train.list = Y.train.list, Y.test.list = Y.test.list)
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
#' kfolds_points <- create_folds(Y, criteria='points', folds=3, seed=1)
#' Y.train.list <- kfolds_points$Y.train.list
#' Y.test.list <- kfolds_points$Y.test.list
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
