
test_that('rotate_scores_and_loadings function works', {
  set.seed(5)
  scores = matrix(rnorm(15), nrow=5)
  intercept = rnorm(6)
  loadings = matrix(rnorm(18), nrow=6)
  results = rotate_scores_and_loadings(scores, intercept, loadings)
  expected_result = readRDS(test_path("fixtures", "rotation_result.rds"))
  expect_equal(round(results$loadings, 4), round(expected_result$loadings, 4))
  expect_equal(round(results$scores, 4), round(expected_result$scores, 4))
  expect_equal(round(results$rotation.matrix, 4), round(expected_result$rotation.matrix, 4))
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

test_that('FQPCA algorithm for quantile.value=0.5 with penalized conquer works', {
  Y = test_data_fqpca()
  results = fqpca(data=Y, npc=1, quantile.value=0.5, seed=5, verbose=FALSE, splines.method='conquer', penalized=T, lambda.ridge=0.1)
  expected_result = readRDS(test_path("fixtures", "fqpca_05_conquer_penalized.rds"))
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

