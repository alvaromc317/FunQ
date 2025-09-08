
test_that('MFQPCA fitted function works', {
  data <- test_data_fosqr_fqpca()
  Y <- data$Y
  regressors <- data$regressors
  cv_result <- fosqr_fqpca_cv_df(Y=Y, regressors=regressors, npc=5, quantile.value = 0.5, return.models=FALSE, verbose.fosqr_fqpca = FALSE, verbose.cv=F, seed=1)
  expected_result = readRDS(test_path("fixtures", "cv_fosqr_fqpca_05.rds"))
  expect_equal(cv_result$error.matrix, expected_result$error.matrix, tolerance=0.01)
})
