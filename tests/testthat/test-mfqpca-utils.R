
test_that('MFQPCA fitted function works', {
  data <- test_data_mfqpca()
  Y <- data[[1]]
  group <- data[[2]]
  cv_result <- mfqpca_cv_df(data=Y, group=group, npc.between = 5, npc.within=5, quantile.value = 0.5, return.models=F, verbose.mfqpca = FALSE, verbose.cv=F, seed=1)
  expected_result = readRDS(test_path("fixtures", "cv_mfqpca_05.rds"))
  expect_equal(cv_result$error.matrix, expected_result$error.matrix, tolerance=0.01)
})
