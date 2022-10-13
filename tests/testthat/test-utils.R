test_that("quantile function works", {
  expect_equal(quantile_function(matrix(c(1, -2, 7), nrow=3), quantile_value=0.2), matrix(c(0.2, 1.6, 1.4), nrow=3))
})

test_that("CVXR Quantile regression function works", {
  set.seed(5)
  x = matrix(rnorm(20), nrow=10)
  y = x %*% matrix(c(1,2), nrow=2) +  rnorm(10)
  solution = qr(x=x, y=y, quantile_value=0.5)
  expect_equal(round(solution, 4), round(c(1.041265, 1.301729), 4))
})

test_that("CVXR Quantile ridge regression function works", {
  set.seed(5)
  x = matrix(rnorm(20), nrow=10)
  y = x %*% matrix(c(1,2), nrow=2) +  rnorm(10)
  solution = qr_ridge(x=x, y=y, R=diag(2), quantile_value=0.5, lambda=10)
  expect_equal(round(solution, 4), round(c(0.01354876, 0.00802739), 4))
})
