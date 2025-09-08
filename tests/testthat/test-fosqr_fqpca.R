
test_that('FOSQR - FQPCA fitted function works', {
  data <- test_data_fosqr_fqpca()
  Y <- data$Y
  regressors <- data$regressors
  results <- fosqr_fqpca(Y = Y, regressors=regressors, npc = 1, quantile.value = 0.5, seed=1)

  Yhat.fosqr <- sweep(regressors %*% t(results$fosqr.loadings), MARGIN = 2, STATS = results$fosqr.intercept, FUN = "+")
  Yhat.fqpca <- results$fqpca.scores %*% t(results$fqpca.loadings)
  Y.hat <- Yhat.fqpca + Yhat.fosqr
  Y.fitted <- fitted(object = results, pve=0.95)
  expect_equal(Y.fitted, Y.hat)
})
