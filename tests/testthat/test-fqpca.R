
test_that('FQFA incorrect inputs', {
  expect_error(fqpca(x=5, n_components=1, quantile_value=0.5)) # Check x
  expect_error(fqpca(x='a', n_components=1, quantile_value=0.5)) # check x
  expect_error(fqpca(x=matrix(0, 5, 5), n_components=1.5, quantile_value=0.5, verbose=FALSE)) # check n_components
  expect_error(fqpca(x=matrix(0, 5, 5), n_components=1, quantile_value=-1, verbose=FALSE)) # check quantile_value
  expect_error(fqpca(x=matrix(0, 5, 5), n_components=1, quantile_value=2, verbose=FALSE)) # check quantile_value
})

test_that('algorithm_normalization function works', {
  set.seed(5)
  Fhat = matrix(rnorm(12), nrow=6)
  Lhat = matrix(rnorm(10), nrow=5)
  results = algorithm_normalization(Fhat, Lhat)
  expected_loadings = matrix(c(-0.84085548, 1.38435934, -1.25549186, 0.07014277, 1.71144087, -0.60290798, -0.68183259, -0.91750891, -0.41267185, 0.19943539, 1.77276147, -1.15781084), nrow=6)
  expected_scores = matrix(c(-1.0803926, -0.1575344, -1.0717600, -0.1389861, -0.5973131, -1.5123884, 0.1667650, -0.1796026, 0.6236010, 0.6522408), nrow=5)
  expected_normalization = as.matrix(1.444052)
  expect_equal(round(results$loadings, 4), round(expected_loadings, 4))
  expect_equal(round(results$scores, 4), round(expected_scores, 4))
  expect_equal(round(results$normalization_matrix, 4), round(expected_normalization, 4))
})

test_that('FQPCA algorithm works', {
  set.seed(5)
  x = matrix(rep(sin(seq(0, 2*pi, length.out=12)), 20), byrow=TRUE, nrow=20)
  x = x + matrix(rnorm(20*12, 0, 0.4), nrow=20)
  x[sample(20*12, as.integer(0.2*20*12))] = NA

  results = fqpca(x=x, n_components=1, quantile_value=0.5, alpha_ridge=1e-12, seed=5, verbose=FALSE)
  expected_loadings = matrix(c(-0.10375779,0.90597597,1.12715232,0.68514343,0.84965556,0.33445475,-0.17404908,-0.93024271,-0.75668444,-0.39242729,-0.57319216,-0.10375779,0.11003609,-0.03325922,1.02089589,-0.71372588,-0.45002830,0.13331065,-0.09538148,-0.47429959,1.81584091,2.30433939,-1.16680408,0.11003609), nrow=12)
  expected_scores = matrix(c(1.00000000,1.00000000,1.00000000,1.00000000,1.00000000,1.00000000,1.00000000,1.00000000,1.00000000,1.00000000,1.00000000,1.00000000 ,1.00000000,1.00000000,1.00000000,1.00000000,1.00000000,1.00000000,1.00000000,1.00000000,0.15425417,-0.34314130,-0.10274230,0.23007249 ,-0.39195760,-0.16872070,-0.23742170,-0.29110932,-0.09589277,-0.27014792,-0.44563142,-0.13951461,-0.31383480,-0.51443796,-0.25703369,-0.05524880 ,-0.08276430,-0.25718975,-0.23816828,-0.12166785), nrow=20)
  expect_equal(round(results$loadings, 4), round(expected_loadings, 4))
  expect_equal(round(results$scores, 4), round(expected_scores, 4))
  expect_equal(results$n_iters, 3)
})

test_that('FQPCA Predictions works', {
  set.seed(5)
  x = matrix(rep(sin(seq(0, 2*pi, length.out=12)), 20), byrow=TRUE, nrow=20)
  x = x + matrix(rnorm(20*12, 0, 0.4), nrow=20)
  x[sample(20*12, as.integer(0.2*20*12))] = NA
  results = fqpca(x=x, n_components=1, quantile_value=0.5, alpha_ridge=1e-12, seed=5, verbose=FALSE)
  newdata = x[5:9, ]
  predictions = predict(results, newdata=newdata)
  expected_predictions = matrix(c(1.00000000,1.00000000,1.00000000,1.00000000,1.00000000,-0.39195760,-0.16872070,-0.23742170,-0.29110932,-0.09589277), ncol=2)
  expect_equal(round(predictions, 4), round(expected_predictions, 4))
})
