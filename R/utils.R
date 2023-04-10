# MATRIX OPERATIONS -----------------------------------------------------------
#' Pivot long a matrix
#'
#' Given a new matrix x of dimensions (n, m) pivot into a (n*m, 3) matrix indexing (row, column, value)
#'
#' @param x Input matrix
#' @return The pivoted matrix
#'
pivot_long_x = function(x)
{
  rowid = column = row = NULL # Required to avoid warning when compiling package
  x_df = data.frame(x)
  colnames(x_df) = 1:ncol(x_df)
  x_df = tibble::rowid_to_column(x_df)
  x_df = tidyr::pivot_longer(data=x_df, cols=-rowid, names_to='column')
  x_df = dplyr::mutate(.data=x_df, column = as.integer(column))
  x_df = dplyr::rename(.data=x_df, row = rowid)
  x_df = as.matrix(x_df)
  return(x_df)
}

#' Pivot short a matrix
#'
#' Given a (n*m, 3) matrix indexing (row, column, value) pivots it into an (n, m) matrix.
#'
#' @param x_long Input matrix
#' @param dimensions array of two elements indicating the dimensions of the original data.
#' @return The pivoted matrix as a tibble
#'
pivot_wide_x = function(x_long, dimensions=NULL)
{
  x_long = as.matrix(x_long)
  if(!is.null(dimensions))
  {
    dimensions = c(max(x_long[,1], dimensions[1]), max(x_long[,2], dimensions[2]))
  } else {
    dimensions = c(max(x_long[,1]), max(x_long[,2]))
  }
  x = matrix(NA, nrow=dimensions[1], ncol=dimensions[2])
  x = base::replace(x, x_long[,c(1,2)], x_long[,3])
  return(x)
}

# QUANTILE BASED FUNCTIONS ----------------------------------------------------

#' quantile_regression objective value
#'
#' Objective value function of a quantile regression model
#'
#' @param quantile_value The quantile considered. Default is the median 0.5.
#' @param x The vector.
#'
#' @return Objective value function for a quantile regression model.
quantile_function = function(x, quantile_value=0.5)
{
  # Quantile check function
  return(0.5 * base::abs(x) + (quantile_value - 0.5) * x)
}

#' Quantile regression model
#'
#' Uses the CVXR library to solve a quantile regression model
#'
#' @param x N times p matrix of predictive variables.
#' @param y N times 1 vector of response.
#' @param quantile_value The quantile considered. Default is the median 0.5.
#' @param solver Solver to be used by the CVXR package in the resolution of the penalized quantile regression models. The speed of the methodology can vary depending on the solver used. There are free alternatives available like 'SCS', 'ECOS' or 'OSQP'. It is also possible to use licensed programs like 'MOSEK' or 'GUROBI'. Default is 'SCS'.
#'
#' @return Beta coefficients of the penalized quantile regression model.
quantile_regression = function(x, y, quantile_value=0.5, solver='SCS')
{
  # Solve a quantile regression model with a ridge type penalty
  n = dim(x)[1]
  m = dim(x)[2]
  beta_var = CVXR::Variable(m)
  objective_function = (1.0 / n) * base::sum(quantile_function(quantile_value=quantile_value, x=(y - x %*% beta_var)))
  objective = CVXR::Minimize(objective_function)
  problem = CVXR::Problem(objective)
  solution = CVXR::solve(problem, solver=solver)
  beta_sol = c(solution$getValue(beta_var))
  results = beta_sol
  return(results)
}


#' Penalized quantile regression
#'
#' Uses the CVXR library to solve a penalized quantile regression model
#'
#' @param x N times p matrix of predictive variables.
#' @param y N times 1 vector of response.
#' @param R Quadratic matrix used to apply a ridge based penalty.
#' @param quantile_value The quantile considered. Default is the median 0.5.
#' @param lambda The hyperparameter controlling the penalization. Default is 1e-3.
#' @param solver Solver to be used by the CVXR package in the resolution of the penalized quantile regression models. The speed of the methodology can vary depending on the solver used. There are free alternatives available like 'SCS', 'ECOS' or 'OSQP'. It is also possible to use licensed programs like 'MOSEK' or 'GUROBI'. Default is 'SCS'.
#'
#' @return Beta coefficients of the penalized quantile regression model.
quantile_regression_ridge = function(x, y, R, quantile_value=0.5, lambda=1e-3, solver='SCS')
{
  # Solve a quantile regression model with a ridge type penalty
  n = dim(x)[1]
  m = dim(x)[2]
  beta_var = CVXR::Variable(m)
  objective_function = (1.0 / n) * base::sum(quantile_function(quantile_value=quantile_value, x=(y - x %*% beta_var)))
  ridge_penalization = lambda * CVXR::quad_form(beta_var, R)
  objective = CVXR::Minimize(objective_function + ridge_penalization)
  problem = CVXR::Problem(objective)
  solution = CVXR::solve(problem, solver=solver)
  beta_sol = c(solution$getValue(beta_var))
  results = beta_sol
  return(results)
}

# TRAIN TEST SPLIT ------------------------------------------------------------

#' Split rows of a given matrix X into train / test
#'
#' Split rows of a given matrix X into train / test based on parameters
#' train_pct and train_size. If both are informed train_pct takes preference
#'
#' @param x N times p matrix of predictive variables.
#' @param train_pct Float number indicating the % of rows used for training
#' @param train_size Integer number indicating number of rows used for training. train_size is superseded by train_pct
#' @param seed Seed for the random generator number. Parameter included for reproducibility purposes. Default is NULL (meaning no seed is assigned).
#'
#' @return a list containing a matrix x_train and a matrix x_test
train_test_split_rows = function(x, train_pct=NULL, train_size=NULL, seed=NULL)
{
  # train_pct takes preference over train_size
  if(is.null(train_pct) & is.null(train_size)){stop('Either train_pct or train_size must be filled.')}
  if(!is.null(seed)){set.seed(seed)}
  if(is.null(train_pct)){train_pct = train_size / nrow(x) }
  x = as.matrix(x)
  mask_x = base::is.na(x)
  curve_NAs = base::rowSums(mask_x) == ncol(x)
  train_idx = caret::createDataPartition(curve_NAs, p = train_pct, list = FALSE)
  x_train = x[train_idx, ]
  x_test = x[-train_idx,]
  results = list(x_train=x_train, x_test=x_test)
  return(results)
}

#' Split points of a given matrix X into train / test
#'
#' Split points of a given matrix X into train / test based on parameters
#' train_pct and train_size. If both are informed train_pct takes preference.
#' This function keeps the dimensions of the original X in both train and test
#' and substitutes the values on either division by NAs
#'
#' @param x N times p matrix of predictive variables.
#' @param train_pct Float number indicating the % of rows used for training
#' @param train_size Integer number indicating number of rows used for training. train_size is superseded by train_pct
#' @param seed Seed for the random generator number. Parameter included for reproducibility purposes. Default is NULL (meaning no seed is assigned).
#'
#' @return a list containing a matrix x_train and a matrix x_test
train_test_split_points = function(x, train_pct=NULL, train_size=NULL, seed=NULL)
{
  value = NULL
  if(is.null(train_pct) & is.null(train_size)){stop('Either train_pct or train_size must be filled.')}
  if(!is.null(seed)){set.seed(seed)}
  if(is.null(train_pct)){train_pct = train_size / length(x) }
  x_long = pivot_long_x(x)
  x_long = tibble::as_tibble(x_long)
  x_long = dplyr::mutate(.data=x_long, is_na = is.na(value))
  train_idx = caret::createDataPartition(x_long$is_na, p = train_pct, list = FALSE)
  train_idx = as.vector(train_idx)
  x_train = x_long[train_idx, c('row', 'column', 'value')]
  x_train = pivot_wide_x(x_long=x_train, dimensions = c(nrow(x), ncol(x)))

  x_test = x_long[-train_idx, c('row', 'column', 'value')]
  x_test = pivot_wide_x(x_long=x_test, dimensions = c(nrow(x), ncol(x)))

  results = list(x_train=x_train, x_test=x_test)
  return(results)
}


#' Split a given matrix X into train / test
#'
#' Splits a given matrix into a train / test split based on two posible
#' criterias: either based on the total number of rows or the total number of
#' data points.
#'
#' @param x N times p matrix of predictive variables.
#' @param criteria Criteria used to divide the data. Valid values are 'rows', which considers the division based on full rows, or 'points', which considers the division based on points within the matrix.
#' @param train_pct Float number indicating the % of rows used for training
#' @param train_size Integer number indicating number of rows used for training. train_size is superseded by train_pct
#' @param seed Seed for the random generator number. Parameter included for reproducibility purposes. Default is NULL (meaning no seed is assigned).
#'
#' @return a list containing a matrix x_train and a matrix x_test
#' @export
#'
#' @examples
#' # Generate a small matrix
#' x = matrix(rnorm(50), nrow=10)
#'
#' # Divide based on full rows
#' train_test_rows = train_test_split(x, criteria='rows', train_pct=0.5)
#' x_train = train_test_rows$x_train
#' x_test = train_test_rows$x_test
#'
#' train_test_points = train_test_split(x, criteria='points', train_pct=0.5)
#'
train_test_split = function(x, criteria='points', train_pct=NULL, train_size=NULL, seed=NULL)
{
  if(criteria == 'points')
  {
    return(train_test_split_points(x=x, train_pct=train_pct, train_size=train_size, seed=seed))
  } else if(criteria == 'rows')
  {
    return(train_test_split_rows(x=x, train_pct=train_pct, train_size=train_size, seed=seed))
  } else{stop('Invalid criteria')}
}

# K-FOLD DIVISION -------------------------------------------------------------

#' Split a given matrix X into k-folds
#'
#' Splits full rows of a given matrix X into k-folds
#'
#' @param x N times p matrix of predictive variables.
#' @param folds Integer number indicating number of folds
#' @param seed Seed for the random generator number. Parameter included for reproducibility purposes. Default is NULL (meaning no seed is assigned).
#'
#' @return A list containing two inside lists, one for training and one for testing. The length of the inside lists is equal to the number of folds
#'
kfold_cv_rows = function(x, folds=3, seed=NULL)
{
  if(!is.null(seed)){set.seed(seed)}
  x = as.matrix(x)
  mask_x = base::is.na(x)
  curve_NAs = base::rowSums(mask_x) == ncol(x)
  test_idx = caret::createFolds(curve_NAs, k=folds, list = FALSE)
  x_train_list = list()
  x_test_list = list()
  for(idx in 1:folds)
  {
    x_train_list[[idx]] = x[test_idx != idx,]
    x_test_list[[idx]] = x[test_idx == idx,]
  }
  results = list(x_train_list=x_train_list, x_test_list=x_test_list)
  return(results)
}

#' Split a given matrix X into k-folds
#'
#' Splits a given matrix X into k-folds. This function keeps the dimensions of
#' the original X in both train and test and substitutes the values on either
#' division by NAs
#'
#' @param x N times p matrix of predictive variables.
#' @param folds Integer number indicating number of folds
#' @param seed Seed for the random generator number. Parameter included for reproducibility purposes. Default is NULL (meaning no seed is assigned).
#'
#' @return A list containing two inside lists, one for training and one for testing. The length of the inside lists is equal to the number of folds
#'
kfold_cv_points = function(x, folds=3, seed=NULL)
{
  value = NULL # Required to avoid warning when compiling package
  if(!is.null(seed)){set.seed(seed)}
  x_long = pivot_long_x(x)
  x_long = tibble::as_tibble(x_long)
  x_long = dplyr::mutate(.data=x_long, is_na = is.na(value))
  test_idx = caret::createFolds(x_long$is_na, k=folds, list = FALSE)
  x_train_list = list()
  x_test_list = list()
  for(idx in 1:folds)
  {
    x_train_list[[idx]] = x_long[test_idx != idx, c('row', 'column', 'value')]
    x_train_list[[idx]] = pivot_wide_x(x_long=x_train_list[[idx]], dimensions = c(nrow(x), ncol(x)))

    x_test_list[[idx]] = x_long[test_idx == idx, c('row', 'column', 'value')]
    x_test_list[[idx]] = pivot_wide_x(x_long=x_test_list[[idx]], dimensions = c(nrow(x), ncol(x)))
  }
  results = list(x_train_list=x_train_list, x_test_list=x_test_list)
  return(results)
}

#' Split a given matrix X into k-folds
#'
#' Splits a given matrix X into k-folds. based on two posible
#' criterias: either based on the total number of rows or the total number of
#' data points.
#'
#' @param x N times p matrix of predictive variables.
#' @param criteria Criteria used to divide the data. Valid values are 'rows', which considers the division based on full rows, or 'points', which considers the division based on points within the matrix.
#' @param folds Integer number indicating number of folds
#' @param seed Seed for the random generator number. Parameter included for reproducibility purposes. Default is NULL (meaning no seed is assigned).
#'
#' @return A list containing two inside lists, one for training and one for testing. The length of the inside lists is equal to the number of folds
#' @export
#'
#' @examples
#' # Generate a small matrix
#' x = matrix(rnorm(50), nrow=10)
#'
#' # Divide based on full rows
#' kfolds_rows = create_folds(x, criteria='rows', folds=3, seed=1)
#' x_train_list = kfolds_rows$x_train_list
#' x_test_list = kfolds_rows$x_test_list
#'
#' kfolds_points = create_folds(x, criteria='points', folds=3, seed=1)
#'
create_folds = function(x, criteria='points', folds=3, seed=NULL)
{
  if(criteria == 'points')
  {
    return(kfold_cv_points(x=x, folds=folds, seed=seed))
  } else if(criteria == 'rows')
  {
    return(kfold_cv_rows(x=x, folds=folds, seed=seed))
  } else{stop('Invalid criteria')}
}

# CROSS VALIDATION ------------------------------------------------------------

#' CROSS VALIDATION
#'
#' Performs cross validation on alpha parameter of fqpca
#
#' @param x The N by T matrix of observed time instants
#' @param n_components Default=2. The number of estimated components.
#' @param quantile_value Default=0.5. The quantile considered.
#' @param periodic Default=TRUE. Boolean indicating if the data is expected to
#'                  be periodic (start coincides with end) or not.
#' @param splines_df Default=10. Degrees of freedom for the splines.
#' @param method Default='fn'. Method used in the resolution of the quantile
#'                regression model. It currently accepts the methods
#'                c('br', 'fn', 'pfn', 'sfn') along with any available solver in
#'                CVXR package like the free 'SCS' or the comercial 'MOSEK'.
#' @param alpha_grid Default=c(0, 1e-16, 1e-14). An array containing the list of
#'                    possible alpha values  (these should be always positive
#'                    numbers).
#' @param n_folds Default=3. Number of folds to be used on cross validation.
#' @param criteria Default='points'. Criteria used to divide the data.
#'                 Valid values are 'rows', which considers the division based
#'                 on full rows, or 'points', which considers the division based
#'                 on points within the matrix. Default is points
#' @param tol Default=1e-3. Tolerance on the convergence of the algorithm.
#'             Smaller values can speed up computation but may affect the
#'             quality of the estimations.
#' @param n_iters Default=30. Maximum number of iterations.
#' @param verbose_fqpca Default=FALSE. Boolean indicating verbosity of the
#'                       fqpca function.
#' @param verbose_cv Default=TRUE. Boolean indicating verbosity of the cross
#'                    validation process.
#' @param seed Default=NULL. Seed for the random generator number.
#'              Parameter included for reproducibility purposes.
#' @return A list containing the matrix of scores, the matrix of loadings, and a secondary list with extra information.
#' @export
#'
#' @examples
#' # Generate fake dataset with 200 observations and 144 time points
#'
#' x = matrix(rep(sin(seq(0, 2*pi, length.out=144)), 200), byrow=TRUE, nrow=200)
#' x = x + matrix(rnorm(200*144, 0, 0.4), nrow=200)
#'
#' # Add missing observations
#' x[sample(200*144, as.integer(0.2*200*144))] = NA
#'
#' cv_result = cross_validation_alpha(x, alpha_grid=c(0, 1e-15), n_folds=2)
#'
cross_validation_alpha = function(x, n_components=2,  quantile_value=0.5, alpha_grid = c(0, 1e-16, 1e-14), n_folds=3, criteria='points', periodic=TRUE, splines_df=10, tol=1e-3, n_iters=50, method='fn', verbose_fqpca=FALSE, verbose_cv=TRUE, seed=NULL)
{
  start_time = Sys.time()
  if(!base::is.null(seed)){base::set.seed(seed)}

  # KFOLDS
  x_folds = create_folds(x, criteria=criteria, folds=n_folds, seed=seed)

  # Initialize error storage
  error_matrix = matrix(-1, nrow=length(alpha_grid), ncol=n_folds)

  # CROSS VALIDATION
  for(i in seq_along(alpha_grid))
  {
    start_loop_time = Sys.time()
    if(verbose_cv){message('Alpha: ', alpha_grid[i], ' ---------------------')}
    alpha_ridge = alpha_grid[i]
    for(j in 1:n_folds)
    {
      if(verbose_cv){message('Fold: ', j)}

      # Access data depending on criteria
      if(criteria=='rows')
      {
        x_train = x_folds$x_train_list[[j]]
        x_test = x_folds$x_test_list[[j]]
      } else if(criteria == 'points')
      {
        x_train = x_folds$x_train_list[[j]]
        x_test =  x_folds$x_test_list[[j]]
      } else{stop('Invalid value for criteria. Valid values are observations or curves')}

      # Execute model
      fqpca_results = fqpca(x=x_train, n_components=n_components,  quantile_value=quantile_value,  periodic=periodic, splines_df=splines_df, method=method, alpha_ridge=alpha_ridge, tol=tol, n_iters=n_iters, verbose=verbose_fqpca, seed=seed)

      # Obtain reconstruction
      if(criteria == 'points')
      {
        x_predicted = fqpca_results$scores %*% t(fqpca_results$loadings)
      }else if(criteria=='rows')
      {
        test_scores = predict.fqpca_object(fqpca_results, x_test)
        x_predicted = test_scores %*% t(fqpca_results$loadings)
      }
      error_matrix[i, j] = mean(quantile_function(quantile_value = quantile_value, x=(x_test - x_predicted)), na.rm=TRUE)
    }
    end_loop_time =  Sys.time()
    if(verbose_cv){message('Alpha ', alpha_ridge, ' execution completed in: ', round(difftime(end_loop_time, start_loop_time, units = "secs"), 2), ' seconds')}
  }
  end_time = Sys.time()
  execution_time = difftime(end_time, start_time, units = "secs")
  return(list(error_matrix = error_matrix, execution_time=execution_time, alpha_grid=alpha_grid, criteria=criteria))
}
