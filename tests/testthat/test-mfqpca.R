
test_that('MFQPCA fitted function works', {
  data <- test_data_mfqpca()
  Y <- data[[1]]
  group <- data[[2]]
  results <- mfqpca(data=Y, group=group, npc.between = 5, npc.within=5, quantile.value = 0.5)

  variability.between = diag(var(results$scores.between[,-1]))
  pve.between <- cumsum(variability.between) / base::sum(variability.between)
  npc.between.fitted =  min(which(pve.between > 0.1))
  variability.within = diag(var(results$scores.within))
  pve.within <- cumsum(variability.within) / base::sum(variability.within)
  npc.within.fitted =  min(which(pve.within > 0.1))
  scores.between = matrix(results$scores.between.full[,1:npc.between.fitted], ncol=npc.between.fitted)
  scores.within = matrix(results$scores.within[,1:npc.within.fitted], ncol=npc.within.fitted)
  loadings.between = matrix(results$loadings.between[,1:npc.between.fitted], ncol=npc.between.fitted)
  loadings.within = matrix(results$loadings.within[,1:npc.within.fitted], ncol=npc.within.fitted)
  Y.fitted = scores.between %*% t(loadings.between) + scores.within %*% t(loadings.within)
  Y.fitted = sweep(Y.fitted, 2, results$intercept, FUN = "+")
  Yhat <- fitted(object = results, pve.between=0.1, pve.within=0.1)

  expect_equal(Y.fitted, Yhat)
})
