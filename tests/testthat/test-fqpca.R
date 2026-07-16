test_that('FQPCA accepts a tf-object as data input', {
  Y = test_data_fqpca(proportion_na = 0)
  Y.tf = tf::tfd(Y, arg = seq(0, 1, length.out = ncol(Y)))
  results = fqpca(data=Y.tf, npc=1, quantile.value=0.5, seed=5, verbose=FALSE)
  expect_s3_class(results, "fqpca_object")
  expect_equal(dim(results$loadings), c(ncol(Y), 1))
})

test_that('FQPCA algorithm for quantile.value=0.1 works', {
  Y = test_data_fqpca()
  results = fqpca(data=Y, npc=1, quantile.value=0.1, seed=5, verbose=FALSE)
  expected_result = readRDS(test_path("fixtures", "fqpca_01.rds"))
  expect_equal(round(results$loadings, 4), round(expected_result$loadings, 4))
  expect_equal(round(results$scores, 4), round(expected_result$scores, 4))
  expect_equal(results$n.iters, expected_result$n.iters)
})

test_that('FQPCA algorithm for quantile.value=0.1 with quantreg works', {
  Y = test_data_fqpca()
  results = fqpca(data=Y, npc=1, quantile.value=0.1, seed=5, verbose=FALSE, splines.method='quantreg')
  expected_result = readRDS(test_path("fixtures", "fqpca_01_quantreg.rds"))
  expect_equal(round(results$loadings, 4), round(expected_result$loadings, 4))
  expect_equal(round(results$scores, 4), round(expected_result$scores, 4))
  expect_equal(results$n.iters, expected_result$n.iters)
})

test_that('FQPCA algorithm for quantile.value=0.5 works', {
  Y = test_data_fqpca()
  results = fqpca(data=Y, npc=1, quantile.value=0.5, seed=5, verbose=FALSE)
  expected_result = readRDS(test_path("fixtures", "fqpca_05.rds"))
  expect_equal(round(results$loadings, 4), round(expected_result$loadings, 4))
  expect_equal(round(results$scores, 4), round(expected_result$scores, 4))
  expect_equal(results$n.iters, expected_result$n.iters)
})

test_that('FQPCA algorithm for quantile.value=0.5 with conquer works', {
  Y = test_data_fqpca()
  results = fqpca(data=Y, npc=1, quantile.value=0.5, seed=5, verbose=FALSE, splines.method='conquer')
  expected_result = readRDS(test_path("fixtures", "fqpca_05_conquer.rds"))
  expect_equal(round(results$loadings, 4), round(expected_result$loadings, 4))
  expect_equal(round(results$scores, 4), round(expected_result$scores, 4))
  expect_equal(results$n.iters, expected_result$n.iters)
})


test_that('FQPCA Fitted values works', {
  Y = test_data_fqpca()
  results = fqpca(data=Y, npc=1, quantile.value=0.5, seed=5, verbose=FALSE)
  fitted.values = fitted(results, pve=0.95)
  fitted2 = sweep(results$scores %*% t(results$loadings), 2, results$intercept, FUN = "+")
  expect_equal(round(fitted.values, 4), round(fitted2, 4))
})

test_that('FQPCA error checking works', {
  Y = test_data_fqpca()
  # check_fqpca_params
  expect_error(fqpca(data=Y, npc=-1), "Expected a positive integer number")
  expect_error(fqpca(data=Y, quantile.value=-0.1), "Expected a float number in \\(0, 1\\)")
  expect_error(fqpca(data=Y, periodic="yes"), "Expected a Boolean value")
  expect_error(fqpca(data=Y, npc=3, splines.df=2), "Expected an integer number larger or equal than 'npc'")
  expect_error(fqpca(data=Y, splines.method='invalid'), "Expected 'conquer' or 'quantreg'")
  expect_error(fqpca(data=Y, penalized="yes"), "Expected a Boolean value")
  expect_error(fqpca(data=Y, lambda.ridge=-1), "Expected a positive number or zero")
  expect_error(fqpca(data=Y, tol=-1), "Expected a positive float number")
  expect_error(fqpca(data=Y, max.iters=-1), "Expected an integer number >= 2")
  expect_error(fqpca(data=Y, seed=-1), "Expected a positive integer number or NULL")

  # Mismatched params
  expect_error(fqpca(data=Y, penalized=TRUE, splines.method='quantreg'), "If penalized is set to  TRUE then splines.method must be 'conquer'")

  # Data check
  df <- data.frame(a=1:5)
  expect_error(fqpca(data=df, colname="a"), "The column 'a' in the data.frame is not a 'tf' object")
  expect_error(fqpca(data=df), "Data is a data.frame, but 'colname' is not provided")
  expect_error(fqpca(data=df, colname="b"), "does not correspond to any column in the data")
  expect_error(fqpca(data="string"), "Data is not a matrix, 'tf' object, or data.frame")
})

test_that('FQPCA algorithm works for extreme quantiles', {
  Y = test_data_fqpca()
  results_01 = fqpca(data=Y, npc=1, quantile.value=0.01, seed=5, verbose=FALSE)
  expect_true(is.matrix(results_01$loadings))

  results_99 = fqpca(data=Y, npc=1, quantile.value=0.99, seed=5, verbose=FALSE)
  expect_true(is.matrix(results_99$loadings))
})

test_that('FQPCA algorithm predicts new data and errors on invalid newdata', {
  Y = test_data_fqpca()
  results = fqpca(data=Y, npc=2, quantile.value=0.5, seed=5, verbose=FALSE)

  preds = predict(results, newdata=Y)
  expect_true(is.matrix(preds))
  expect_equal(dim(preds), c(nrow(Y), 2))

  # Invalid predict calls
  expect_error(predict(results, newdata=matrix(NA, nrow=5, ncol=ncol(Y))), "newdata contains no information")
  expect_error(predict.fqpca_object("not_an_object", newdata=Y), "must be of class fqpca_object")
})

test_that('FQPCA Fitted values works for pve=0 edge case', {
  Y = test_data_fqpca()
  results = fqpca(data=Y, npc=1, quantile.value=0.5, seed=5, verbose=FALSE)

  # Test with pve such that n.components = 0
  fitted.values.zero = fitted(results, pve=0) # Using pve=0 to trigger n.components=0
  expect_equal(fitted.values.zero[1, 1], results$intercept[1])

  expect_error(fitted.fqpca_object("not_an_object"), "must be of class fqpca_object")
})

test_that('FQPCA plot method returns ggplot object', {
  Y = test_data_fqpca()
  results = fqpca(data=Y, npc=2, quantile.value=0.5, seed=5, verbose=FALSE)

  p = plot(results, pve=0.99)
  expect_true(inherits(p, "ggplot"))

  p_zero = plot(results, pve=0)
  expect_true(inherits(p_zero, "ggplot"))

  expect_error(plot.fqpca_object("not_an_object"), "must be of class fqpca_object")
})

# CROSS VALIDATION PROCESS ----------------------------------------------------

test_that("fqpca_cv_lambda incorrect inputs detection works", {
  Y = test_data_fqpca()
  expect_error(fqpca_cv_lambda(data=Y, lambda.grid=c(0, 1e-15), n.folds=3, verbose.cv=FALSE, splines.method='aaa')) # method
  expect_error(fqpca_cv_lambda(data=Y, lambda.grid=c(0, 1e-15), n.folds=2.5, verbose.cv=FALSE, splines.method='conquer')) # n.folds
})

test_that("fqpca_cv_lambda based on points function works", {
  Y = test_data_fqpca()
  cv_result = fqpca_cv_lambda(data=Y, lambda.grid=c(0, 1e-3), n.folds=3, verbose.cv=FALSE, splines.method='conquer', seed=5, return.models=FALSE)
  expected_result = readRDS(test_path("fixtures", "cv_lambda_05.rds"))
  expect_equal(cv_result$error.matrix, expected_result$error.matrix, tolerance=0.01)
})

test_that("fqpca_cv_df based on points function works", {
  Y = test_data_fqpca()
  cv_result = fqpca_cv_df(data=Y, quantile.value=0.9, splines.df.grid=c(5, 8), n.folds=3, verbose.cv=FALSE, splines.method='conquer', seed=5, return.models=FALSE)
  expected_result = readRDS(test_path("fixtures", "cv_df_09.rds"))
  expect_equal(cv_result$error.matrix, expected_result$error.matrix, tolerance=0.01)
})

test_that("fqpca_cv_lambda errors on wrong inputs", {
  Y = test_data_fqpca()
  expect_error(fqpca_cv_lambda(data=Y, lambda.grid=c(0, 1e-15), n.folds=3.5, verbose.cv=FALSE), "must be an integer number")
  expect_error(fqpca_cv_lambda(data=Y, lambda.grid=c(0, 1e-15), n.folds=2, criteria='invalid'), "Invalid criteria")
  expect_error(fqpca_cv_lambda(data=Y, lambda.grid=c(0, 1e-15), n.folds=2, splines.method='quantreg', penalized=TRUE), "splines.method must be conquer")
  expect_error(fqpca_cv_lambda(data=Y, lambda.grid=c(0, 1e-15), n.folds=2, penalized=FALSE), "penalized parameter must be set to TRUE")
})

test_that("fqpca_cv_lambda with criteria rows works and retains models", {
  Y = test_data_fqpca()
  cv_result = fqpca_cv_lambda(data=Y, lambda.grid=c(0, 1e-3), n.folds=2, verbose.cv=TRUE, splines.method='conquer', criteria='rows', seed=5, return.models=TRUE)
  expect_true(!is.null(cv_result$error.matrix))
  expect_true(length(cv_result$list.models) > 0)
})

test_that("fqpca_cv_df errors on wrong inputs", {
  Y = test_data_fqpca()
  expect_error(fqpca_cv_df(data=Y, splines.df.grid=c(5, 8), n.folds=3.5, verbose.cv=FALSE), "must be an integer number")
  expect_error(fqpca_cv_df(data=Y, splines.df.grid=c(5, 8), n.folds=2, criteria='invalid'), "Invalid criteria")
  expect_error(fqpca_cv_df(data=Y, splines.df.grid=c(5.5, 8), n.folds=2), "must be a positive integer array")
})

test_that("fqpca_cv_df with criteria rows works and retains models", {
  Y = test_data_fqpca()
  cv_result = fqpca_cv_df(data=Y, quantile.value=0.9, splines.df.grid=c(5, 8), n.folds=2, verbose.cv=TRUE, splines.method='conquer', criteria='rows', seed=5, return.models=TRUE)
  expect_true(!is.null(cv_result$error.matrix))
  expect_true(length(cv_result$list.models) > 0)
})
