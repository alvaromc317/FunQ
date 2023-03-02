
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
  results = fqpca(x=x, n_components=1, quantile_value=0.5, alpha_ridge=1e-12, seed=5, verbose=FALSE)
  expected_loadings = matrix(c(-0.12544755,0.91253186,0.92591851,0.82582945,0.93836285,0.30817723,-0.15524797,-0.83675118,-1.11461378,-0.84664697,-0.34319767,-0.12544755,0.11003609,-0.03325922,1.02089589,-0.71372588,-0.45002830,0.13331065,-0.09538148,-0.47429959,1.81584091,2.30433939,-1.16680408,0.11003609), nrow=12)
  expected_scores = matrix(c(1.00000000,1.00000000,1.00000000,1.00000000,1.00000000,1.00000000,1.00000000,1.00000000,1.00000000,1.00000000,1.00000000,1.00000000,1.00000000,1.00000000,1.00000000,1.00000000,1.00000000,1.00000000,1.00000000,1.00000000,0.35136909,-0.14602638,0.09437262,0.42718741,-0.19484268,0.02839422,-0.04030678,-0.09399440,0.10122215,-0.07303300,-0.24851650,0.05760031,-0.11671988,-0.31732303,-0.05991877,0.14186612,0.11435062,-0.06007483,-0.04105336,0.07544707), nrow=20)
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
  expected_predictions = matrix(c(1.00000000,1.00000000,1.00000000,1.00000000,1.00000000,-0.19484268,0.02839422,-0.04030678,-0.09399440,0.10122215), ncol=2)
  expect_equal(round(predictions, 4), round(expected_predictions, 4))
})
