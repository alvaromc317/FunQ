test_that("quantile function works", {
  expect_equal(quantile_function(matrix(c(1, -2, 7), nrow=3), quantile_value=0.2), matrix(c(0.2, 1.6, 1.4), nrow=3))
})

test_that("CVXR Quantile ridge regression function works", {
  set.seed(5)
  x = matrix(rnorm(20), nrow=10)
  y = x %*% matrix(c(1,2), nrow=2) +  rnorm(10)
  solution = qr_ridge(x=x, y=y, R=diag(2), quantile_value=0.5, lambda=10)
  expect_equal(round(solution, 5), round(c(0.01354876, 0.00802739), 5))
})

test_that('algorithm_normalization function works', {
  set.seed(5)
  Fhat = matrix(rnorm(12), nrow=6)
  Lhat = matrix(rnorm(10), nrow=5)
  results = algorithm_normalization(Fhat, Lhat)
  expected_loadings = matrix(c(-0.84085548, 1.38435934, -1.25549186, 0.07014277, 1.71144087, -0.60290798, -0.68183259, -0.91750891, -0.41267185, 0.19943539, 1.77276147, -1.15781084), nrow=6)
  expected_scores = matrix(c(-1.0803926, -0.1575344, -1.0717600, -0.1389861, -0.5973131, -1.5123884, 0.1667650, -0.1796026, 0.6236010, 0.6522408), nrow=5)
  expected_normalization = as.matrix(1.444052)
  expect_equal(round(results$Fhat, 4), round(expected_loadings, 4))
  expect_equal(round(results$Lhat, 4), round(expected_scores, 4))
  expect_equal(round(results$normalization_matrix, 4), round(expected_normalization, 4))
})

test_that('FQFA incorrect inputs', {
  expect_error(fqpca(x=5, n_components=1, quantile_value=0.5)) # Check x
  expect_error(fqpca(x='a', n_components=1, quantile_value=0.5)) # check x
  expect_error(fqpca(x=matrix(0, 5, 5), n_components=1.5, quantile_value=0.5, verbose=FALSE)) # check n_components
  expect_error(fqpca(x=matrix(0, 5, 5), n_components=1, quantile_value=-1, verbose=FALSE)) # check quantile_value
  expect_error(fqpca(x=matrix(0, 5, 5), n_components=1, quantile_value=2, verbose=FALSE)) # check quantile_value
})

test_that('QFA algorithm works', {
  set.seed(5)
  x = matrix(rep(sin(seq(0, 2*pi, length.out=12)), 20), byrow=TRUE, nrow=20)
  x = x + matrix(rnorm(20*12, 0, 0.4), nrow=20)
  x[sample(20*12, as.integer(0.2*20*12))] = NA

  results = fqpca(x=x, n_components=1, quantile_value=0.5, seed=5, verbose=FALSE)
  expected_loadings = matrix(c(-0.12518748, 0.91245330, 0.92836363, 0.82455808, 0.93730186, 0.30849490, -0.15466429, -0.83787132, -1.11032493, -0.84120435, -0.34595214, -0.12518748, 0.11003081, -0.03325746, 1.02132530, -0.71262706, -0.45000976, 0.13329688, -0.10406764, -0.47427733, 1.81574647, 2.30423128, -1.16673678, 0.11003081), nrow=12)
  expected_scores = matrix(c(1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 0.34225282, -0.14231325, 0.09258555, 0.41611763, -0.18981294, 0.02760994, -0.03931977, -0.09162334, 0.09855796, -0.07074595, -0.24216121, 0.05600382, -0.11376266, -0.30919321, -0.05842662, 0.13815427, 0.11134802, -0.05858253, -0.04007066, 0.07338212), nrow=20)
  expect_equal(round(results$loadings, 4), round(expected_loadings, 4))
  expect_equal(round(results$scores, 4), round(expected_scores, 4))
  expect_equal(results$n_iters, 3)
})

test_that('FQPCA Predictions works', {
  set.seed(5)
  x = matrix(rep(sin(seq(0, 2*pi, length.out=12)), 20), byrow=TRUE, nrow=20)
  x = x + matrix(rnorm(20*12, 0, 0.4), nrow=20)
  x[sample(20*12, as.integer(0.2*20*12))] = NA
  results = fqpca(x=x, n_components=1, quantile_value=0.5, seed=5, verbose=FALSE)
  newdata = x[5:9, ]
  predictions = predict(results, newdata=newdata)
  expected_predictions = matrix(c(1, 1, 1, 1, 1, -0.18981294, 0.02760994, -0.03931977, -0.09162334, 0.09855796), ncol=2)
  expect_equal(round(predictions, 4), round(expected_predictions, 4))

  newdata = matrix(x[5, ], nrow=1)
  predictions = predict(results, newdata=newdata)
  expected_predictions = matrix(c(1, -0.1898129), ncol=2)
  expect_equal(round(predictions, 4), round(expected_predictions, 4))
})
