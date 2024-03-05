# QUANTILE RELATED TESTS ------------------------------------------------------

test_that("quantile function works", {
  expect_equal(quantile_function(matrix(c(1, -2, 7), nrow=3), quantile.value=0.2), matrix(c(0.2, 1.6, 1.4), nrow=3))
})

test_that("CVXR Quantile regression function works", {
  set.seed(5)
  data = test_data_cvxr()
  x = data$x
  y = data$y
  solution = quantile_regression(x=x, y=y, quantile.value=0.5, solver='SCS')
  expect_equal(round(solution, 4), round(c(0.8278325, 1.8020721), 4), tolerance=0.01)
})

test_that("CVXR Quantile ridge regression function works", {
  set.seed(5)
  data = test_data_cvxr()
  x = data$x
  y = data$y
  solution = quantile_regression_ridge(x=x, y=y, R=diag(2), quantile.value=0.5, lambda=10, solver='SCS')
  expect_equal(round(solution, 4), round(c(0.001840445, 0.027260196), 4), tolerance=0.01)
})

# TRAIN TEST SPLIT TESTS ------------------------------------------------------

test_that("train_test_split based on rows function; train portion works", {
  set.seed(5)
  x = matrix(rnorm(50), nrow=10)
  train_test_rows = train_test_split(x, criteria='rows', train.pct=0.5, seed=1)
  Y.train = train_test_rows$Y.train
  expected_output_train = matrix(c(-0.84085548,1.2276303,0.9005119,0.315915,1.5500604, 1.38435934,-0.8017795,0.9418694,1.109694,-0.8024232, 0.07014277,-0.1575344,0.7067611,1.217104,1.8956680, -0.47216639,-0.5973131,1.4185891,-1.009533,-0.8870085, -0.28577363,0.2408173,-0.6570821,-1.762186,-0.7243285), nrow=5, byrow=T)
  expect_equal(round(Y.train, 4), round(expected_output_train, 4))
})

test_that("train_test_split based on rows function; test portion works", {
  set.seed(5)
  x = matrix(rnorm(50), nrow=10)
  train_test_rows = train_test_split(x, criteria='rows', train.pct=0.5, seed=1)
  Y.test = train_test_rows$Y.test
  expected_output_test = matrix(c(-1.2554919,-1.0803926,1.4679619,2.2154606,-0.07457892, 1.7114409,-1.0717600,0.8190089,1.4792218,-0.45656894, -0.6029080,-0.1389861,-0.2934818,0.9515738,0.56222336, -0.6353713,-2.1839668,1.4987738,-2.0004727,-0.46024458, 0.1381082,-0.2593554,-0.8527954,-0.1426081,-0.06921116), nrow=5, byrow=T)
  expect_equal(round(Y.test, 4), round(expected_output_test, 4))
})

test_that("train_test_split based on points function; train portion works", {
  set.seed(5)
  x = matrix(rnorm(50), nrow=10)
  train_test_rows = train_test_split(x, criteria='points', train.pct=0.5, seed=1)
  Y.train = train_test_rows$Y.train
  expected_output_train = matrix(c(-0.8408555,NA,0.9005119,0.3159150,NA, 1.3843593,-0.8017795,NA,1.1096942,-0.80242318, NA,-1.0803926,NA,2.2154606,-0.07457892), nrow=3, byrow=T)
  expect_equal(round(Y.train[1:3,], 4), round(expected_output_train, 4))
})

test_that("train_test_split based on points function; test portion works", {
  set.seed(5)
  x = matrix(rnorm(50), nrow=10)
  train_test_rows = train_test_split(x, criteria='points', train.pct=0.5, seed=1)
  Y.test = train_test_rows$Y.test
  expected_output_test = matrix(c(NA, NA, -1.2554919, 1.2276303, NA, NA, NA, 0.9418694, 1.4679619, NA, NA, NA, 1.5500604, NA, NA), nrow=3, byrow=F)
  expect_equal(round(Y.test[1:3,], 4), round(expected_output_test, 4))
})

# KFOLDS TESTS ----------------------------------------------------------------

test_that("create_folds based on rows function; train portion works", {
  set.seed(5)
  x = matrix(rnorm(50), nrow=10)
  kfolds = create_folds(x, criteria='rows',folds=3, seed=1)
  true_partial_output = list(train1=round(kfolds$Y.train.list[[1]], 4),
                             test1=round(kfolds$Y.test.list[[1]], 4))
  expected_partial_output = list(train1=round(matrix(c(0.07014277,-0.1575344,0.7067611,1.2171036,1.89566795, 1.71144087,-1.0717600,0.8190089,1.4792218,-0.45656894, -0.60290798,-0.1389861,-0.2934818,0.9515738,0.56222336, -0.47216639,-0.5973131,1.4185891,-1.0095326,-0.88700851, -0.28577363,0.2408173,-0.6570821,-1.7621859,-0.72432849, 0.13810822,-0.2593554,-0.8527954,-0.1426081,-0.06921116), nrow=6, byrow=T), 4),
                                 test1=round(matrix(c(-0.8408555,1.2276303,0.9005119,0.315915,1.55006037, 1.3843593,-0.8017795,0.9418694,1.109694,-0.80242318, -1.2554919,-1.0803926,1.4679619,2.215461,-0.07457892, -0.6353713,-2.1839668,1.4987738,-2.000473,-0.46024458), nrow=4, byrow=T), 4))
  expect_equal(expected_partial_output, true_partial_output)
})

test_that("create_folds based on points function; train portion works", {
  set.seed(5)
  x = matrix(rnorm(50), nrow=10)
  kfolds = create_folds(x, criteria='points',folds=3, seed=1)
  true_partial_output = list(train1=round(kfolds$Y.train.list[[1]][1:3,], 4),
                             test1=round(kfolds$Y.test.list[[1]][1:3,], 4))
  expected_partial_output = list(train1=round(matrix(c(-0.8408555,NA,NA,0.315915,NA, 1.3843593,-0.8017795,0.9418694,1.109694,-0.80242318, NA,NA,1.4679619,2.215461,-0.07457892), nrow=3, byrow=T), 4),
                                 test1=round(matrix(c(NA,1.227630,0.9005119,NA,1.55006, NA,NA,NA,NA,NA, -1.255492,-1.080393,NA,NA,NA), nrow=3, byrow=T), 4))
  expect_equal(expected_partial_output, true_partial_output)
})

# CROSS VALIDATION PROCESS ----------------------------------------------------

test_that("cross_validation_alpha incorrect inputs detection works", {
  Y = test_data_fqpca()
  expect_error(cross_validation_alpha(Y=Y, alpha.grid=c(0, 1e-15), n.folds=3, verbose.cv=FALSE, method='fn')) # method
  expect_error(cross_validation_alpha(Y=Y, alpha.grid=c(0, 1e-15), n.folds=2.5, verbose.cv=FALSE, method='SCS')) # n.folds
})

test_that("cross_validation_alpha based on points function works", {
  Y = test_data_fqpca()
  cv_result = cross_validation_alpha(Y=Y, alpha.grid=c(0, 1e-15), n.folds=3, verbose.cv=FALSE, method='SCS', seed=5)
  expected_result = readRDS(test_path("fixtures", "cv_alpha_05.rds"))
  expect_equal(cv_result$error.matrix, expected_result$error.matrix, tolerance=0.01)
})

test_that("cross_validation_df based on points function works", {
  Y = test_data_fqpca()
  cv_result = cross_validation_df(Y=Y, quantile.value=0.9, splines.df.grid=c(5, 10), n.folds=3, verbose.cv=FALSE, method='SCS', seed=5)
  expected_result = readRDS(test_path("fixtures", "cv_df_09.rds"))
  expect_equal(cv_result$error.matrix, expected_result$error.matrix, tolerance=0.01)
})

# ERROR FUNCTIONS -------------------------------------------------------------

test_that("quantile error 0.1 works", {
  Y = test_data_fqpca()
  Y.pred = test_data_fqpca(seed=10)
  expect_equal(round(quantile_error(Y=Y, Y.pred=Y.pred, quantile.value=0.1), 5), round(0.20384035, 5))
})

test_that("quantile error 0.5 works", {
  Y = test_data_fqpca()
  Y.pred = test_data_fqpca(seed=10)
  expect_equal(round(quantile_error(Y=Y, Y.pred=Y.pred, quantile.value=0.5), 5), round(0.22865256, 5))
})

test_that("proportion under quantile works", {
  Y = test_data_fqpca()
  Y.pred = test_data_fqpca(seed=10)
  expect_equal(round(proportion_under_quantile(Y=Y, Y.pred=Y.pred), 5), round(0.46753247, 5))
})



