library(testthat)
library(FunQ)

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

test_that('FOSQR - FQPCA input validation functions', {
  data <- test_data_fosqr_fqpca()
  Y <- data$Y
  regressors <- data$regressors
  
  # check_Y logic
  expect_error(fosqr_fqpca(Y = "string", regressors=regressors, npc=1))
  expect_error(fosqr_fqpca(Y = array(1:24, dim=c(2,3,4)), regressors=regressors, npc=1))
  
  # check_regressors logic
  expect_error(fosqr_fqpca(Y = Y, regressors="string", npc=1))
  
  reg_na <- regressors
  reg_na[1, 1] <- NA
  expect_error(fosqr_fqpca(Y = Y, regressors=reg_na, npc=1))
  
  reg_const <- cbind(regressors, rep(1, nrow(regressors)))
  expect_error(fosqr_fqpca(Y = Y, regressors=reg_const, npc=1))
  
  # check_fosqr_fqpca_params
  expect_error(fosqr_fqpca(Y = Y, regressors=regressors, npc="1"))
  expect_error(fosqr_fqpca(Y = Y, regressors=regressors, npc=0))
  expect_error(fosqr_fqpca(Y = Y, regressors=regressors, quantile.value="0.5"))
  expect_error(fosqr_fqpca(Y = Y, regressors=regressors, quantile.value=1.5))
  expect_error(fosqr_fqpca(Y = Y, regressors=regressors, periodic=1))
  expect_error(fosqr_fqpca(Y = Y, regressors=regressors, splines.df=1)) # < npc or not matching
  expect_error(fosqr_fqpca(Y = Y, regressors=regressors, npc=2, splines.df=1))
  expect_error(fosqr_fqpca(Y = Y, regressors=regressors, splines.method="invalid"))
  expect_error(fosqr_fqpca(Y = Y, regressors=regressors, tol="1e-3"))
  expect_error(fosqr_fqpca(Y = Y, regressors=regressors, max.iters=-1))
  expect_error(fosqr_fqpca(Y = Y, regressors=regressors, verbose=1))
  expect_error(fosqr_fqpca(Y = Y, regressors=regressors, seed="1"))
})

test_that('FOSQR - FQPCA S3 Predict', {
  data <- test_data_fosqr_fqpca()
  Y <- data$Y
  regressors <- data$regressors
  # Fit with reduced sample
  fit <- fosqr_fqpca(Y = Y[1:100, ], regressors=regressors[1:100, , drop=FALSE], npc = 1, quantile.value = 0.5, seed=1)
  
  expect_error(predict(fit, newdata=Y[101:150, ], newregressors=matrix(1, nrow=50, ncol=1)), "constant column")
  
  # Prediction with correct sizes
  preds <- predict(fit, newdata=Y[101:150, ], newregressors=regressors[101:150, , drop=FALSE])
  expect_equal(dim(preds), c(50, 1))
  expect_true(is.numeric(preds))
  
  # Invalid model dimension
  expect_error(predict(fit, newdata=Y[101:150, 1:10], newregressors=regressors[101:150, , drop=FALSE]))
})

test_that('FOSQR - FQPCA S3 Variance', {
  data <- test_data_fosqr_fqpca()
  Y <- data$Y
  regressors <- data$regressors
  fit <- fosqr_fqpca(Y = Y[1:20,], regressors=regressors[1:20,,drop=FALSE], npc = 1, splines.df=5, quantile.value = 0.5, seed=1)
  
  # Use n.bootstrap=5 to run quickly
  variances <- compute_variance(fit, pve=0.95, n.bootstrap=5)
  expect_type(variances, "list")
  expect_true("variance" %in% names(variances))
  expect_true("var.correction" %in% names(variances))
  expect_true("var.analytical" %in% names(variances))
  
  # Variance computation with pve=0
  variances_0 <- compute_variance(fit, pve=0, n.bootstrap=5)
  expect_equal(sum(variances_0$var.correction), 0)
})

test_that('FOSQR - FQPCA S3 Variance matches numeric reference', {
  data <- test_data_fosqr_fqpca()
  Y <- data$Y
  regressors <- data$regressors
  fit <- fosqr_fqpca(Y = Y[1:20,], regressors=regressors[1:20,,drop=FALSE], npc = 1, splines.df=5, quantile.value = 0.5, seed=1)

  set.seed(1)
  variances <- compute_variance(fit, pve=0.95, n.bootstrap=5)
  expected_result <- readRDS(test_path("fixtures", "fosqr_fqpca_variance_05.rds"))
  expect_equal(round(variances$variance, 4), round(expected_result$variance, 4))
  expect_equal(round(variances$var.analytical, 4), round(expected_result$var.analytical, 4))
  expect_equal(round(variances$var.correction, 4), round(expected_result$var.correction, 4))
})

test_that('FOSQR - FQPCA S3 methods reject objects of the wrong class', {
  expect_error(predict.fosqr_fqpca_object("not_an_object", newdata=matrix(1), newregressors=matrix(1)), "must be of class fosqr_fqpca_object")
  expect_error(fitted.fosqr_fqpca_object("not_an_object"), "must be of class fosqr_fqpca_object")
  expect_error(compute_variance.fosqr_fqpca_object("not_an_object"), "must be of class fosqr_fqpca_object")
  expect_error(plot.fosqr_fqpca_object("not_an_object"), "must be of class fosqr_fqpca_object")
})

test_that('FOSQR - FQPCA S3 Plot', {
  data <- test_data_fosqr_fqpca()
  Y <- data$Y
  regressors <- data$regressors
  fit <- fosqr_fqpca(Y = Y[1:20,], regressors=regressors[1:20,,drop=FALSE], npc = 1, splines.df=5, quantile.value = 0.5, seed=1)
  
  p <- plot(fit)
  expect_s3_class(p, "ggplot")
  
  p0 <- plot(fit, pve=0)
  expect_s3_class(p0, "ggplot")
})

test_that('FOSQR - FQPCA Edge Cases', {
  data <- test_data_fosqr_fqpca()
  Y <- data$Y
  regressors <- data$regressors
  # periodic = FALSE and method = quantreg
  fit_edge <- fosqr_fqpca(Y = Y[1:20,], regressors=regressors[1:20,,drop=FALSE], npc=1, periodic=FALSE, splines.method="quantreg", splines.df=5, max.iters=2, seed=1)
  expect_s3_class(fit_edge, "fosqr_fqpca_object")
  
  # Check differing quantile values
  fit_edge2 <- fosqr_fqpca(Y = Y[1:20,], regressors=regressors[1:20,,drop=FALSE], npc=1, periodic=TRUE, splines.method="conquer", splines.df=5, quantile.value=0.9, max.iters=2, seed=1)
  expect_s3_class(fit_edge2, "fosqr_fqpca_object")
})

# CROSS VALIDATION PROCESS ----------------------------------------------------

test_that('FOSQR - FQPCA CV df functions correctly', {
  data <- test_data_fosqr_fqpca()
  Y <- data$Y
  regressors <- data$regressors

  # Run a lightweight cross-validation to not block tests
  cv_res <- fosqr_fqpca_cv_df(
    Y = Y[1:40, ],
    regressors = regressors[1:40, , drop=FALSE],
    npc = 1,
    splines.df.grid = c(4, 5),
    n.folds = 2,
    criteria = 'points',
    max.iters = 2,
    verbose.cv = FALSE,
    return.models = TRUE,
    seed = 1
  )

  expect_type(cv_res, "list")
  expect_true("error.matrix" %in% names(cv_res))
  expect_true("list.models" %in% names(cv_res))
  expect_equal(dim(cv_res$error.matrix), c(2, 2))

  # Check criteria = 'rows'
  cv_res_rows <- fosqr_fqpca_cv_df(
    Y = Y[1:40, ],
    regressors = regressors[1:40, , drop=FALSE],
    npc = 1,
    splines.df.grid = c(5),
    n.folds = 2,
    criteria = 'rows',
    max.iters = 2,
    verbose.cv = FALSE,
    return.models = FALSE,
    seed = 1
  )

  expect_type(cv_res_rows, "list")
  expect_equal(dim(cv_res_rows$error.matrix), c(1, 2))
  expect_equal(length(cv_res_rows$list.models), 0)
})

test_that('FOSQR - FQPCA CV df matches numeric reference', {
  data <- test_data_fosqr_fqpca()
  Y <- data$Y
  regressors <- data$regressors

  cv_res <- fosqr_fqpca_cv_df(
    Y = Y[1:40, ],
    regressors = regressors[1:40, , drop=FALSE],
    npc = 1,
    splines.df.grid = c(4, 5),
    n.folds = 2,
    criteria = 'points',
    max.iters = 2,
    verbose.cv = FALSE,
    return.models = FALSE,
    seed = 1
  )
  expected_result <- readRDS(test_path("fixtures", "cv_fosqr_fqpca_df_05.rds"))
  expect_equal(round(cv_res$error.matrix, 4), round(expected_result$error.matrix, 4))
})

test_that('FOSQR - FQPCA CV constraints and validation', {
  data <- test_data_fosqr_fqpca()
  Y <- data$Y
  regressors <- data$regressors

  expect_error(fosqr_fqpca_cv_df(Y = Y, regressors = regressors, n.folds = 2.5), "integer number")
  expect_error(fosqr_fqpca_cv_df(Y = Y, regressors = regressors, criteria = "invalid"), "Invalid criteria")
  expect_error(fosqr_fqpca_cv_df(Y = Y, regressors = regressors, splines.df.grid = c(2.5, 4)), "positive integer array")
})
