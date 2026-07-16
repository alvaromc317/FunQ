# QUANTILE RELATED TESTS ------------------------------------------------------

test_that("quantile function works", {
  expect_equal(check_loss(matrix(c(1, -2, 7), nrow=3), quantile.value=0.2), matrix(c(0.2, 1.6, 1.4), nrow=3))
})

test_that("proportion_under_quantile works", {
  Y <- matrix(c(1, 2, 3, 4), nrow = 2)
  Y.pred <- matrix(c(1, 1, 3, 5), nrow = 2)
  # Y <= Y.pred: TRUE, FALSE, TRUE, TRUE -> 3/4
  expect_equal(proportion_under_quantile(Y, Y.pred), 0.75)

  # NAs are dropped from the proportion (na.rm = TRUE)
  Y.na <- Y
  Y.na[1, 1] <- NA
  expect_equal(proportion_under_quantile(Y.na, Y.pred), 2/3)
})
