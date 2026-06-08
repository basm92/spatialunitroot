test_that("euclidean_distances works", {
  n <- 10
  coords <- matrix(rnorm(n * 2), ncol = 2)
  dmat <- euclidean_distances(coords)

  expect_equal(dim(dmat), c(n, n))
  expect_equal(diag(dmat), rep(0, n))
  expect_true(isSymmetric(dmat))
  expect_true(all(dmat >= 0))
})

test_that("haversine_distances works", {
  n <- 10
  coords <- cbind(runif(n, -90, 90), runif(n, -180, 180))
  dmat <- haversine_distances(coords)

  expect_equal(dim(dmat), c(n, n))
  expect_equal(diag(dmat), rep(0, n))
  expect_true(all(dmat >= 0 & dmat <= 1))
})

test_that("normalize_distances works", {
  dmat <- matrix(runif(100, 0, 5), 10, 10)
  diag(dmat) <- 0
  dmat <- (dmat + t(dmat)) / 2

  normed <- normalize_distances(dmat)
  expect_equal(dim(normed), dim(dmat))
  expect_equal(max(normed), 1)
})
