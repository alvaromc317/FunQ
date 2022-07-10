test_that('QFA incorrect inputs', {
  expect_error(qfa(x=5, n_components=1, quantile_value=0.5)) # Check x
  expect_error(qfa(x='a', n_components=1, quantile_value=0.5)) # check x
  expect_error(qfa(x=matrix(0, 5, 5), n_components=1.5, quantile_value=0.5, verbose=FALSE)) # check n_components
  expect_error(qfa(x=matrix(0, 5, 5), n_components=1, quantile_value=-1, verbose=FALSE)) # check quantile_value
  expect_error(qfa(x=matrix(0, 5, 5), n_components=1, quantile_value=2, verbose=FALSE)) # check quantile_value
})

test_that('QFA algorithm works', {
  set.seed(5)
  x = matrix(rep(sin(seq(0, 2*pi, length.out=12)), 20), byrow=TRUE, nrow=20)
  x = x + matrix(rnorm(20*12, 0, 0.4), nrow=20)
  results = qfa(x=x, n_components=1, quantile_value=0.5, seed=5, verbose=FALSE)
  expected_loadings = matrix(c(-0.2778827, 1.2412438, 1.3241667, 1.2113886, 1.1596838, 0.4078259, -0.2567481, -1.0531221, -1.5924436, -1.2874327, -0.4780583, -0.2305246), nrow=12)
  expected_scores = matrix(c(0.8236004, 0.7287768, 0.6261189, 0.6433792, 0.7499094, 0.7872014, 0.6796629, 0.7174746, 0.4764496, 0.7247509, 0.7497037, 0.5927589, 0.9316852, 0.8277844, 0.7648711, 0.5381713, 0.5321416, 0.6810075, 0.6675863, 0.6312321), nrow=20)
  expect_equal(round(results$loadings, 4), round(expected_loadings, 4))
  expect_equal(round(results$scores, 4), round(expected_scores, 4))
  expect_equal(results$n_iters, 3)
})

test_that('QFA Predictions works', {
  set.seed(5)
  x = matrix(rep(sin(seq(0, 2*pi, length.out=12)), 20), byrow=TRUE, nrow=20)
  x = x + matrix(rnorm(20*12, 0, 0.4), nrow=20)
  results = qfa(x=x, n_components=1, quantile_value=0.5, seed=5, verbose=FALSE)
  newdata = x[5:7, ]
  predictions = predict(results, newdata=newdata)
  expected_predictions = matrix(c(0.7499094, 0.7872014, 0.6796629), ncol=1)
  expect_equal(round(predictions, 4), round(expected_predictions, 4))

  newdata = matrix(x[5, ], nrow=1)
  predictions = predict(results, newdata=newdata)
  expected_predictions = matrix(c(0.7499094), ncol=1)
  expect_equal(round(predictions, 4), round(expected_predictions, 4))
})
