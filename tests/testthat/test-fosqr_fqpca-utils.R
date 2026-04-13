library(testthat)
library(FunQ)

test_that('FOSQR - FQPCA CV df functions correctly', {
  data <- test_data_fosqr_fqpca()
  Y <- data$Y
  regressors <- data$regressors
  
  # Run a lightweight cross-validation to not block tests
  cv_res <- fosqr_fqpca_cv_df(
    Y = Y[1:40, ],
    regressors = regressors[1:40, , drop=FALSE],
    npc = 1,
    splines.df.grid = c(4, 5),
    n.folds = 2,
    criteria = 'points',
    max.iters = 2,
    verbose.cv = FALSE,
    return.models = TRUE,
    seed = 1
  )
  
  expect_type(cv_res, "list")
  expect_true("error.matrix" %in% names(cv_res))
  expect_true("list.models" %in% names(cv_res))
  expect_equal(dim(cv_res$error.matrix), c(2, 2))
  
  # Check criteria = 'rows'
  cv_res_rows <- fosqr_fqpca_cv_df(
    Y = Y[1:40, ],
    regressors = regressors[1:40, , drop=FALSE],
    npc = 1,
    splines.df.grid = c(5),
    n.folds = 2,
    criteria = 'rows',
    max.iters = 2,
    verbose.cv = FALSE,
    return.models = FALSE,
    seed = 1
  )
  
  expect_type(cv_res_rows, "list")
  expect_equal(dim(cv_res_rows$error.matrix), c(1, 2))
  expect_equal(length(cv_res_rows$list.models), 0)
})

test_that('FOSQR - FQPCA CV constraints and validation', {
  data <- test_data_fosqr_fqpca()
  Y <- data$Y
  regressors <- data$regressors

  expect_error(fosqr_fqpca_cv_df(Y = Y, regressors = regressors, n.folds = 2.5), "integer number")
  expect_error(fosqr_fqpca_cv_df(Y = Y, regressors = regressors, criteria = "invalid"), "Invalid criteria")
  expect_error(fosqr_fqpca_cv_df(Y = Y, regressors = regressors, splines.df.grid = c(2.5, 4)), "positive integer array")
})
