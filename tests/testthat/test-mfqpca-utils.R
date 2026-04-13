
test_that('MFQPCA fitted function works', {
  data <- test_data_mfqpca()
  Y <- data[[1]]
  group <- data[[2]]
  cv_result <- mfqpca_cv_df(data=Y, group=group, npc.between = 5, npc.within=5, quantile.value = 0.5, return.models=F, verbose.mfqpca = FALSE, verbose.cv=F, seed=1)
  expected_result = readRDS(test_path("fixtures", "cv_mfqpca_05.rds"))
  expect_equal(cv_result$error.matrix, expected_result$error.matrix, tolerance=0.01)
})

test_that("mfqpca_cv_df errors on wrong inputs", {
  data <- test_data_mfqpca()
  Y <- data[[1]]
  group <- data[[2]]
  
  expect_error(mfqpca_cv_df(data=Y, group=group, splines.df.grid=c(5, 8), n.folds=3.5, verbose.cv=FALSE), "must be an integer number")
  expect_error(mfqpca_cv_df(data=Y, group=group, splines.df.grid=c(5, 8), n.folds=2, criteria='invalid'), "Invalid criteria")
  expect_error(mfqpca_cv_df(data=Y, group=group, splines.df.grid=c(5.5, 8), n.folds=2), "must be a positive integer array")
})

test_that("mfqpca_cv_df with criteria rows works, output verbosely, and retains models", {
  data <- test_data_mfqpca()
  Y <- data[[1]]
  group <- data[[2]]
  
  cv_result = mfqpca_cv_df(data=Y, group=group, quantile.value=0.5, npc.between=1, npc.within=1, splines.df.grid=c(5, 8), n.folds=2, verbose.cv=TRUE, splines.method='conquer', criteria='rows', seed=5, return.models=TRUE)
  expect_true(!is.null(cv_result$error.matrix))
  expect_true(length(cv_result$list.models) > 0)
})
