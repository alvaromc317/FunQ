# QUANTILE RELATED TESTS ------------------------------------------------------

test_that("quantile function works", {
  expect_equal(quantile_function(matrix(c(1, -2, 7), nrow=3), quantile.value=0.2), matrix(c(0.2, 1.6, 1.4), nrow=3))
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
  x = matrix(1:50, nrow=10)
  train_test_rows = train_test_split(x, criteria='points', train.pct=0.5, seed=1)
  Y.train = train_test_rows$Y.train
  expected_output_train = matrix(
    c(1, 2, NA, 4, 5, NA, 7, 8, 9, 10, NA, 12, 13, NA, 15, 16, NA, NA, NA, NA, 21, NA, 23, 24, NA, 26, 27, NA, NA, NA, 31, NA, 33, NA, NA, NA, NA, NA, NA, NA, NA, 42, NA, 44, 45, NA, NA, 48, 49, 50),
    nrow=10)
  expect_equal(round(Y.train, 4), round(expected_output_train, 4))
})

test_that("train_test_split based on points function; test portion works", {
  set.seed(5)
  x = matrix(1:50, nrow=10)
  train_test_rows = train_test_split(x, criteria='points', train.pct=0.5, seed=1)
  Y.test = train_test_rows$Y.test
  expected_output_test = matrix(
    c(NA, NA, 3, NA, NA, 6, NA, NA, NA, NA, 11, NA, NA, 14, NA, NA, 17, 18, 19, 20, NA, 22,  NA, NA, 25, NA, NA, 28, 29, 30, NA, 32, NA, 34, 35, 36, 37, 38, 39, 40, 41, NA, 43, NA,  NA, 46, 47, NA, NA, NA),
    nrow=10)
  expect_equal(round(Y.test, 4), round(expected_output_test, 4))
})

# KFOLDS TESTS ----------------------------------------------------------------

test_that("create_folds based on rows function; train portion works", {
  set.seed(5)
  x = matrix(1:50, nrow=10)
  kfolds = create_folds(x, criteria='rows',folds=3, seed=1)
  true_partial_output = list(train1=round(kfolds$Y.train.list[[1]], 4),
                             test1=round(kfolds$Y.test.list[[1]], 4))
  expected_partial_output = list(
    train1 = matrix(c(2, 4, 5, 6, 7, 10, 12, 14, 15, 16, 17, 20, 22, 24, 25, 26, 27, 30, 32, 34, 35, 36, 37, 40, 42, 44, 45, 46, 47, 50), nrow=6),
    test1 = matrix(c(1, 3, 8, 9, 11, 13, 18, 19, 21, 23, 28, 29, 31, 33, 38, 39, 41, 43, 48, 49), nrow=4)
  )
  expect_equal(expected_partial_output, true_partial_output)
})

test_that("create_folds based on points function; train portion works", {
  set.seed(5)
  x = matrix(1:50, nrow=10)
  kfolds = create_folds(x, criteria='points',folds=3, seed=1)
  true_partial_output = list(train1=round(kfolds$Y.train.list[[1]], 4),
                             test1=round(kfolds$Y.test.list[[1]], 4))
  expected_partial_output = list(
    train1 = matrix(c(NA, 2, 3, 4, NA, NA, 7, NA, 9, 10, 11, NA, 13, NA, 15, 16, NA, NA, NA, 20, 21, 22,  NA, NA, 25, NA, 27, 28, 29, NA, 31, 32, NA, 34, NA, 36, NA, 38, NA, NA, NA, NA, 43, 44,  45, 46, 47, 48, 49, 50),
                    nrow=10),
    test1 = matrix(c(1, NA, NA, NA, 5, 6, NA, 8, NA, NA, NA, 12, NA, 14, NA, NA, 17, 18, 19, NA, NA, NA, 23, 24, NA, 26, NA, NA, NA, 30, NA, NA, 33, NA, 35, NA, 37, NA, 39, 40, 41, 42, NA, NA, NA, NA, NA, NA, NA, NA),
                   nrow=10)
  )
  expect_equal(expected_partial_output, true_partial_output)
})

