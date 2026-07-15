
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

test_that('MFQPCA Fitted values works for pve.between=0 and pve.within=0', {
  data <- test_data_mfqpca()
  Y <- data[[1]]
  group <- data[[2]]
  results <- mfqpca(data=Y, group=group, npc.between=1, npc.within=1, quantile.value=0.5, seed=1, verbose=FALSE)
  
  # if both are 0, it should just be intercept
  Yhat <- fitted(results, pve.between=0, pve.within=0)
  expect_equal(Yhat[1,1], results$intercept[1])
  
  # Invalid predict calls
  expect_error(fitted.mfqpca_object("not_an_object"), "must be of class mfqpca_object")
})

test_that('MFQPCA error checking works', {
  data <- test_data_mfqpca()
  Y <- data[[1]]
  group <- data[[2]]
  
  # check_params
  expect_error(mfqpca(data=Y, group=group, npc.between=-1), "Expected a positive integer number")
  expect_error(mfqpca(data=Y, group=group, npc.within=-1), "Expected a positive integer number")
  expect_error(mfqpca(data=Y, group=group, quantile.value=-0.1), "Expected a float number in \\(0, 1\\)")
  expect_error(mfqpca(data=Y, group=group, periodic="yes"), "Expected a Boolean value")
  expect_error(mfqpca(data=Y, group=group, npc.between=3, npc.within=3, splines.df=2), "Expected an integer number larger or equal than")
  expect_error(mfqpca(data=Y, group=group, splines.method='invalid'), "Expected 'conquer' or 'quantreg'")
  expect_error(mfqpca(data=Y, group=group, tol=-1), "Expected a positive float number")
  expect_error(mfqpca(data=Y, group=group, max.iters=-1), "Expected a positive integer number")
  expect_error(mfqpca(data=Y, group=group, seed=-1), "Expected a positive integer number or NULL")

  # Data structure checking
  df <- data.frame(a=1:5, b=1:5)
  expect_error(mfqpca(data=df, group="a", colname="b"), "not a 'tf' object")
  expect_error(mfqpca(data=df, group="a"), "not provided as a single string")
  expect_error(mfqpca(data=df, group="a", colname="x"), "does not match any column")
  
  # Group checking
  expect_error(mfqpca(data=Y, group="a string"), "For matrix input, 'group' cannot be a single string")
  expect_error(mfqpca(data=Y, group=rep(1, 1000)), "Length of 'group' \\(1000\\) does not match the number of rows")
})

test_that('MFQPCA algorithm works for extreme quantiles', {
  data <- test_data_mfqpca()
  Y <- data[[1]]
  group <- data[[2]]
  
  results_01 = mfqpca(data=Y, group=group, npc.between=1, npc.within=1, quantile.value=0.01, seed=5, verbose=FALSE)
  expect_true(is.matrix(results_01$loadings.between))
  
  results_99 = mfqpca(data=Y, group=group, npc.between=1, npc.within=1, quantile.value=0.99, seed=5, verbose=FALSE)
  expect_true(is.matrix(results_99$loadings.within))
})

test_that('MFQPCA algorithm predicts new data and errors on invalid newdata', {
  data <- test_data_mfqpca()
  Y <- data[[1]]
  group <- data[[2]]
  results <- mfqpca(data=Y, group=group, npc.between=1, npc.within=1, quantile.value=0.5, seed=5, verbose=FALSE)
  
  preds = predict(results, newdata=Y, newdata.group=group)
  expect_true(is.list(preds))
  expect_true(is.matrix(preds$scores.between))
  
  # Invalid predict calls
  expect_error(predict(results, newdata=matrix(NA, nrow=5, ncol=ncol(Y)), newdata.group=group[1:5]), "newdata contains no information")
  expect_error(predict.mfqpca_object("not_an_object", newdata=Y, newdata.group=group), "must be of class mfqpca_object")
})

test_that('MFQPCA plot method returns ggplot object', {
  data <- test_data_mfqpca()
  Y <- data[[1]]
  group <- data[[2]]
  results <- mfqpca(data=Y, group=group, npc.between=1, npc.within=1, quantile.value=0.5, seed=5, verbose=FALSE)
  
  p = plot(results, pve.between=0.99, pve.within=0.99)
  expect_true(inherits(p, "ggplot"))
  
  p_zero = plot(results, pve.between=0, pve.within=0)
  expect_true(inherits(p_zero, "ggplot"))
  
  expect_error(plot.mfqpca_object("not_an_object"), "must be of class mfqpca_object")
})
