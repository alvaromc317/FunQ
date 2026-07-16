test_that('build_tensor_matrix_sparse matches the dense tensor matrix construction', {
  set.seed(42)
  n.obs <- 4
  n.time <- 6
  npc <- 2
  Y.axis <- seq(0, 1, length.out = n.time)
  spline.basis <- pbs::pbs(Y.axis, degree = 3, df = 4, intercept = TRUE, periodic = TRUE, Boundary.knots = c(0, 1))
  scores <- matrix(rnorm(n.obs * npc), nrow = n.obs, ncol = npc)
  Y.mask <- matrix(TRUE, n.obs, n.time)
  Y.mask[1, 2] <- FALSE
  Y.mask[3, c(1, 4)] <- FALSE

  dense <- cbind(1, build_tensor_matrix(scores = scores, Y.mask = Y.mask, spline.basis = spline.basis, intercept = TRUE))
  sparse <- build_tensor_matrix_sparse(scores = scores, Y.mask = Y.mask, spline.basis = spline.basis)

  # Custom helper to convert CSR format to a dense matrix without S4 method dispatch issues
  csr_to_dense <- function(x) {
    nrow <- x@dimension[1]
    ncol <- x@dimension[2]
    res <- matrix(0, nrow = nrow, ncol = ncol)
    for (i in 1:nrow) {
      start <- x@ia[i]
      end <- x@ia[i + 1] - 1
      if (start <= end) {
        cols <- x@ja[start:end]
        res[i, cols] <- x@ra[start:end]
      }
    }
    res
  }

  sparse.dense <- csr_to_dense(sparse)

  expect_equal(dim(sparse.dense), dim(dense))
  expect_equal(unname(sparse.dense), unname(dense), tolerance = 1e-10)
})

test_that('fit_tensor_quantile_regression uses the sparse solver path on large problems', {
  # Cross the use.sparse threshold (length(Y.vector) * ((npc+1) * n.basis) > 5e6) with
  # method='quantreg' so fit_tensor_quantile_regression routes through
  # build_tensor_matrix_sparse() / quantreg::rq.fit.sfn(), which is otherwise never
  # exercised on the small datasets used by the rest of the suite.
  set.seed(1)
  n.obs <- 60
  n.time <- 5000
  npc <- 1
  n.basis <- 10

  Y.axis <- seq(0, 1, length.out = n.time)
  spline.basis <- pbs::pbs(Y.axis, degree = 3, df = n.basis, intercept = TRUE, periodic = TRUE, Boundary.knots = c(0, 1))
  scores <- matrix(rnorm(n.obs * npc), nrow = n.obs, ncol = npc)
  true.loadings <- spline.basis %*% matrix(rnorm(ncol(spline.basis)), ncol = 1)
  Y <- outer(scores[, 1], true.loadings[, 1]) + matrix(rnorm(n.obs * n.time, sd = 0.1), n.obs, n.time)
  Y.mask <- matrix(TRUE, n.obs, n.time)
  Y.list <- lapply(seq_len(n.obs), function(j) Y[j, Y.mask[j, ]])
  Y.vector <- unlist(Y.list, use.names = FALSE)

  expect_true(length(Y.vector) * ((npc + 1) * ncol(spline.basis)) > 5e6)

  result <- fit_tensor_quantile_regression(
    Y.vector = Y.vector,
    Y.mask = Y.mask,
    scores = scores,
    spline.basis = spline.basis,
    quantile.value = 0.5,
    method = 'quantreg')

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(ncol(spline.basis), npc + 1))
  expect_true(all(is.finite(result)))
})
