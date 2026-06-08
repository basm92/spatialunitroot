test_that("sigma_lbm is symmetric", {
  n <- 20
  coords <- matrix(rnorm(n * 2), ncol = 2)
  dmat <- euclidean_distances(coords)
  dmat <- normalize_distances(dmat)

  slbm <- sigma_lbm(dmat)
  expect_equal(dim(slbm), c(n, n))
  expect_true(isSymmetric(slbm))
})

test_that("sigma_lbm_dm is double-centered and symmetric", {
  n <- 20
  coords <- matrix(rnorm(n * 2), ncol = 2)
  dmat <- euclidean_distances(coords)
  dmat <- normalize_distances(dmat)

  slbm_dm <- sigma_lbm_dm(dmat)
  expect_equal(dim(slbm_dm), c(n, n))
  expect_true(isSymmetric(slbm_dm))

  # Row and column means should be near zero
  row_means <- rowMeans(slbm_dm)
  col_means <- colMeans(slbm_dm)
  expect_true(all(abs(row_means) < 1e-12))
  expect_true(all(abs(col_means) < 1e-12))
})

test_that("sigma_dm works", {
  n <- 20
  coords <- matrix(rnorm(n * 2), ncol = 2)
  dmat <- euclidean_distances(coords)
  dmat <- normalize_distances(dmat)

  sig <- sigma_dm(dmat, 1.0)
  expect_equal(dim(sig), c(n, n))
  expect_true(isSymmetric(sig))
})

test_that("get_R returns top q eigenvectors", {
  n <- 30; q <- 15
  coords <- matrix(rnorm(n * 2), ncol = 2)
  dmat <- euclidean_distances(coords)
  dmat <- normalize_distances(dmat)

  sig <- sigma_lbm_dm(dmat)
  R <- get_R(sig, q)
  expect_equal(dim(R), c(n, q))

  # Should be approximately orthonormal
  expect_equal(t(R) %*% R, diag(q), tolerance = 1e-10)
})

test_that("getcbar returns reasonable values", {
  n <- 20
  coords <- matrix(rnorm(n * 2), ncol = 2)
  dmat <- euclidean_distances(coords)
  dmat <- normalize_distances(dmat)

  c_small <- getcbar(0.95, dmat)
  c_large <- getcbar(0.01, dmat)

  expect_true(c_small > 0)
  expect_true(c_large > 0)
  expect_true(c_small < c_large)
})
