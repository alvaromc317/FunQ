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
  cv_result = fqpca_cv_df(data=Y, quantile.value=0.9, splines.df.grid=c(5, 10), n.folds=3, verbose.cv=FALSE, splines.method='conquer', seed=5, return.models=FALSE)
  expected_result = readRDS(test_path("fixtures", "cv_df_09.rds"))
  expect_equal(cv_result$error.matrix, expected_result$error.matrix, tolerance=0.01)
})

# ERROR FUNCTIONS -------------------------------------------------------------

test_that("quantile error 0.1 works", {
  Y = test_data_fqpca()
  Y.pred = test_data_fqpca(seed=10)
  expect_equal(round(quantile_error(Y=Y, Y.pred=Y.pred, quantile.value=0.1), 5), round(0.20384035, 5))
})

test_that("quantile error 0.5 works", {
  Y = test_data_fqpca()
  Y.pred = test_data_fqpca(seed=10)
  expect_equal(round(quantile_error(Y=Y, Y.pred=Y.pred, quantile.value=0.5), 5), round(0.22865256, 5))
})

test_that("proportion under quantile works", {
  Y = test_data_fqpca()
  Y.pred = test_data_fqpca(seed=10)
  expect_equal(round(proportion_under_quantile(Y=Y, Y.pred=Y.pred), 5), round(0.46753247, 5))
})



