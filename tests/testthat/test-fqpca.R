
test_that('FQFA incorrect inputs detection works', {
  expect_error(fqpca(x=5, n_components=1, quantile_value=0.5)) # Check x
  expect_error(fqpca(x='a', n_components=1, quantile_value=0.5)) # check x
  expect_error(fqpca(x=matrix(0, 5, 5), n_components=1.5, quantile_value=0.5, verbose=FALSE)) # check n_components
  expect_error(fqpca(x=matrix(0, 5, 5), n_components=1, quantile_value=-1, verbose=FALSE)) # check quantile_value
  expect_error(fqpca(x=matrix(0, 5, 5), n_components=1, quantile_value=2, verbose=FALSE)) # check quantile_value
  expect_error(fqpca(x=matrix(0, 5, 5), n_components=1, quantile_value=2, verbose=FALSE, method='AAA')) # check quantile_value
})

test_that('rotate_scores_and_loadings function works', {
  set.seed(5)
  loadings = matrix(rnorm(18), nrow=6)
  scores = matrix(rnorm(15), nrow=5)
  results = rotate_scores_and_loadings(loadings, scores)
  expected_result = readRDS(test_path("fixtures", "rotation_result.rds"))
  expect_equal(round(results$loadings, 4), round(expected_result$loadings, 4))
  expect_equal(round(results$scores, 4), round(expected_result$scores, 4))
  expect_equal(round(results$rotation.matrix, 4), round(expected_result$rotation.matrix, 4))
})

test_that('FQPCA algorithm works for quantile.value=0.1', {
  Y = test_data_fqpca()
  results = fqpca(Y=Y, npc=1, quantile.value=0.1, seed=5, verbose=FALSE)
  expected_result = readRDS(test_path("fixtures", "fqpca_01.rds"))
  expect_equal(round(results$loadings, 4), round(expected_result$loadings, 4))
  expect_equal(round(results$scores, 4), round(expected_result$scores, 4))
  expect_equal(results$n.iters, expected_result$n.iters)
})

test_that('FQPCA algorithm works for quantile.value=0.5', {
  Y = test_data_fqpca()
  results = fqpca(Y=Y, npc=1, quantile.value=0.5, seed=5, verbose=FALSE)
  expected_result = readRDS(test_path("fixtures", "fqpca_05.rds"))
  expect_equal(round(results$loadings, 4), round(expected_result$loadings, 4))
  expect_equal(round(results$scores, 4), round(expected_result$scores, 4))
  expect_equal(results$n.iters, expected_result$n.iters)
})

test_that('FQPCA Predictions works', {
  Y = test_data_fqpca()
  results = fqpca(Y=Y, npc=1, quantile.value=0.5, seed=5, verbose=FALSE)
  expected_result = readRDS(test_path("fixtures", "fqpca_05.rds"))
  newdata = Y[5:9, ]
  predictions = predict(results, newdata=newdata)
  expected_predictions = expected_result$scores[5:9,]
  expect_equal(round(predictions, 4), round(expected_predictions, 4))
})
