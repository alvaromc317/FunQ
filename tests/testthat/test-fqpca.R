
test_that('FQFA incorrect inputs', {
  expect_error(fqpca(x=5, n_components=1, quantile_value=0.5)) # Check x
  expect_error(fqpca(x='a', n_components=1, quantile_value=0.5)) # check x
  expect_error(fqpca(x=matrix(0, 5, 5), n_components=1.5, quantile_value=0.5, verbose=FALSE)) # check n_components
  expect_error(fqpca(x=matrix(0, 5, 5), n_components=1, quantile_value=-1, verbose=FALSE)) # check quantile_value
  expect_error(fqpca(x=matrix(0, 5, 5), n_components=1, quantile_value=2, verbose=FALSE)) # check quantile_value
  expect_error(fqpca(x=matrix(0, 5, 5), n_components=1, quantile_value=2, verbose=FALSE, method='AAA')) # check quantile_value
})

test_that('algorithm_normalization function works', {
  set.seed(5)
  Fhat = matrix(rnorm(12), nrow=6)
  Lhat = matrix(rnorm(10), nrow=5)
  results = algorithm_normalization(Fhat, Lhat)
  expected_loadings = matrix(c(-0.80684783,1.43012178,-1.23490910,0.06019556,1.62302115,-0.54516005,-0.68183259,-0.91750891,-0.41267185,0.19943539,1.77276147,-1.15781084), nrow=6)
  expected_scores = matrix(c(-1.0803926,-0.1575344,-1.0717600,-0.1389861,-0.5973131,-1.4625116,0.2166418,-0.1297258,0.6734778,0.7021177), nrow=5)
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
  results = fqpca(x=x, n_components=1, quantile_value=0.5, seed=5, verbose=FALSE)
  expected_loadings = matrix(c(-0.14922418,0.83067279,0.89341446,0.89453901,0.81433555,0.35614450,-0.19893906,-0.74567344,-1.12742675,-0.92568372, -0.45669027,-0.14922418,-0.16433113,1.32398882,0.15210724,0.51990480,1.44931314,-0.07888021,-0.59762701,-0.48200042, -0.79115868,-0.25820286,-2.55162550,-0.16433113), nrow=12)
  expected_scores = matrix(c(1.000000000,1.000000000,1.000000000,1.000000000,1.000000000,1.000000000,1.000000000,1.000000000,1.000000000, 1.000000000,1.000000000,1.000000000,1.000000000,1.000000000,1.000000000,1.000000000,1.000000000,1.000000000, 1.000000000,1.000000000,0.050721321,0.065495862,-0.122164026,-0.005534437,-0.010349780,-0.031494479,-0.062044514, -0.113918655,-0.163990574,-0.083162737,-0.123615812,-0.065724076,0.450269860,0.208171892,0.228643569,0.065449802,0.055262073,-0.189511617,-0.090807301,-0.061696372), nrow=20)
  expect_equal(round(results$loadings, 4), round(expected_loadings, 4))
  expect_equal(round(results$scores, 4), round(expected_scores, 4))
  expect_equal(results$n_iters, 4)
})

test_that('FQPCA Predictions works', {
  set.seed(5)
  x = matrix(rep(sin(seq(0, 2*pi, length.out=12)), 20), byrow=TRUE, nrow=20)
  x = x + matrix(rnorm(20*12, 0, 0.4), nrow=20)
  x[sample(20*12, as.integer(0.2*20*12))] = NA
  results = fqpca(x=x, n_components=1, quantile_value=0.5, seed=5, verbose=FALSE)
  newdata = x[5:9, ]
  predictions = predict(results, newdata=newdata)
  expected_predictions = matrix(c(1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, -0.01034978, -0.03149448, -0.06204451, -0.11391865, -0.16399057), ncol=2)
  expect_equal(round(predictions, 4), round(expected_predictions, 4))
})
