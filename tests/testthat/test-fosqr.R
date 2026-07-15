library(testthat)
library(FunQ)

test_that('FOSQR fitted function works', {
  data <- test_data_fosqr_fqpca(seed = 1)
  Y <- data$Y
  regressors <- data$regressors
  
  # Run fosqr model
  results <- fosqr(Y = Y, as.matrix(regressors), quantile.value = 0.5, splines.df = 10, seed = 1)
  
  # Manual calculation
  Yhat.fosqr <- sweep(regressors %*% t(results$fosqr.loadings), MARGIN = 2, STATS = results$fosqr.intercept, FUN = "+")
  
  # Fitted method
  Y.fitted <- fitted(object = results)
  
  expect_equal(Y.fitted, Yhat.fosqr)
  expect_s3_class(results, "fosqr_object")
})

test_that('FOSQR input validation functions', {
  data <- test_data_fosqr_fqpca(seed = 2)
  Y <- data$Y
  regressors <- data$regressors
  
  # invalid Y
  expect_error(fosqr(Y = "string", regressors=regressors))
  
  # invalid regressors
  expect_error(fosqr(Y = Y, regressors="string"))
  reg_na <- regressors
  reg_na[1, 1] <- NA
  expect_error(fosqr(Y = Y, regressors=reg_na))
  
  reg_const <- cbind(regressors, rep(1, nrow(regressors)))
  expect_error(fosqr(Y = Y, regressors=reg_const))
  
  # invalid params
  expect_error(fosqr(Y = Y, regressors=regressors, quantile.value="0.5"))
  expect_error(fosqr(Y = Y, regressors=regressors, quantile.value=1.5))
  expect_error(fosqr(Y = Y, regressors=regressors, periodic=1))
  expect_error(fosqr(Y = Y, regressors=regressors, splines.df=0)) # < 1
  expect_error(fosqr(Y = Y, regressors=regressors, splines.df=1.5))
  expect_error(fosqr(Y = Y, regressors=regressors, splines.method="invalid"))
  expect_error(fosqr(Y = Y, regressors=regressors, verbose=1))
  expect_error(fosqr(Y = Y, regressors=regressors, seed="1"))
})

test_that('FOSQR S3 Predict', {
  data <- test_data_fosqr_fqpca(seed = 3)
  Y <- data$Y
  regressors <- data$regressors
  # Fit with reduced sample
  fit <- fosqr(Y = Y[1:100, ], regressors=regressors[1:100, , drop=FALSE], quantile.value = 0.5, seed=1)
  
  # Predict requires the length to match newregressors columns
  expect_error(predict(fit, newregressors=matrix(1, nrow=50, ncol=1)), "constant column")
  
  # Predict dimension checks
  preds <- predict(fit, newregressors=regressors[101:150, , drop=FALSE])
  expect_equal(dim(preds), c(50, ncol(Y)))
  expect_true(is.numeric(preds))
  
  # Invalid model dimension for regressors matching
  expect_error(predict(fit, newregressors=cbind(regressors[101:150, , drop=FALSE], regressors[101:150, , drop=FALSE])))
})

test_that('FOSQR S3 Variance', {
  data <- test_data_fosqr_fqpca(seed = 4)
  Y <- data$Y
  regressors <- data$regressors
  
  # Use small spline df down to run fast
  fit <- fosqr(Y = Y[1:20,], regressors=regressors[1:20,,drop=FALSE], splines.df=5, quantile.value = 0.5, seed=1)
  
  variances <- compute_variance(fit)
  expect_type(variances, "list")
  expect_true("variance" %in% names(variances))
  expect_true("var.analytical" %in% names(variances))
  
  expect_equal(dim(variances$variance), c(ncol(Y), ncol(regressors) + 1))
})

test_that('FOSQR S3 Plot', {
  data <- test_data_fosqr_fqpca(seed = 5)
  Y <- data$Y
  regressors <- data$regressors
  fit <- fosqr(Y = Y[1:20,], regressors=regressors[1:20,,drop=FALSE], splines.df=5, quantile.value = 0.5, seed=1)
  
  p <- plot(fit)
  expect_s3_class(p, "ggplot")
})

test_that('FOSQR Edge Cases', {
  data <- test_data_fosqr_fqpca(seed = 6)
  Y <- data$Y
  regressors <- data$regressors
  
  # periodic = FALSE and method = quantreg
  fit_edge <- fosqr(Y = Y[1:20,], regressors=regressors[1:20,,drop=FALSE], periodic=FALSE, splines.method="quantreg", splines.df=5, seed=1)
  expect_s3_class(fit_edge, "fosqr_object")
  
  # Extreme quantile value
  fit_edge2 <- fosqr(Y = Y[1:20,], regressors=regressors[1:20,,drop=FALSE], periodic=TRUE, splines.method="conquer", splines.df=5, quantile.value=0.9, seed=1)
  expect_s3_class(fit_edge2, "fosqr_object")
})
