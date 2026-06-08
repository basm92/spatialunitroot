test_that("spurtransform lbmgls works", {
  n <- 30
  y <- rnorm(n)
  x <- rnorm(n)
  coords <- matrix(rnorm(n * 2), ncol = 2)
  dat <- data.frame(y = y, x = x, s_1 = coords[,1], s_2 = coords[,2])

  transformed <- spurtransform(~ y + x, data = dat, coords = ~ s_1 + s_2,
                                prefix = "h_", transformation = "lbmgls")

  expect_true("h_y" %in% names(transformed))
  expect_true("h_x" %in% names(transformed))
  expect_equal(length(transformed$h_y), n)
})

test_that("spurtransform nn works", {
  n <- 30
  y <- rnorm(n)
  coords <- matrix(rnorm(n * 2), ncol = 2)
  dat <- data.frame(y = y, s_1 = coords[,1], s_2 = coords[,2])

  transformed <- spurtransform(~ y, data = dat, coords = ~ s_1 + s_2,
                                prefix = "h_", transformation = "nn")

  expect_true("h_y" %in% names(transformed))
  expect_equal(sum(is.na(transformed$h_y)), 0)
})

test_that("spurtransform cluster works", {
  n <- 30
  y <- rnorm(n)
  cl <- rep(1:3, each = 10)
  coords <- matrix(rnorm(n * 2), ncol = 2)
  dat <- data.frame(y = y, cluster = cl, s_1 = coords[,1], s_2 = coords[,2])

  transformed <- spurtransform(~ y, data = dat, coords = ~ s_1 + s_2,
                                prefix = "h_", transformation = "cluster",
                                cluster = "cluster")

  expect_true("h_y" %in% names(transformed))
})

test_that("spurtransform error for missing radius with iso", {
  n <- 20
  y <- rnorm(n)
  coords <- matrix(rnorm(n * 2), ncol = 2)
  dat <- data.frame(y = y, s_1 = coords[,1], s_2 = coords[,2])

  expect_error(
    spurtransform(~ y, data = dat, coords = ~ s_1 + s_2,
                  prefix = "h_", transformation = "iso"),
    "adius"
  )
})
