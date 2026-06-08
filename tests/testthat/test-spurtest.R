test_that("spurtest i1 returns valid structure", {
  set.seed(123)
  n <- 30; q <- 10; nrep <- 1000

  y <- rnorm(n)
  coords <- matrix(rnorm(n * 2), ncol = 2)
  dat <- data.frame(y = y, s_1 = coords[,1], s_2 = coords[,2])

  result <- spurtest(y ~ 1, data = dat, coords = ~ s_1 + s_2,
                     type = "i1", latlong = FALSE, q = q, nrep = nrep)

  expect_s3_class(result, c("spur_test_i1", "spur_test"))
  expect_named(result$critical_values, c("1%", "5%", "10%"))
  expect_true(result$p_value >= 0 && result$p_value <= 1)
  expect_true(result$statistic > 0)
  expect_true(result$ha_parm > 0)
})

test_that("spurtest i0 returns valid structure", {
  set.seed(123)
  n <- 30; q <- 10; nrep <- 1000

  y <- rnorm(n)
  coords <- matrix(rnorm(n * 2), ncol = 2)
  dat <- data.frame(y = y, s_1 = coords[,1], s_2 = coords[,2])

  result <- spurtest(y ~ 1, data = dat, coords = ~ s_1 + s_2,
                     type = "i0", latlong = FALSE, q = q, nrep = nrep)

  expect_s3_class(result, c("spur_test_i0", "spur_test"))
  expect_true(result$p_value >= 0 && result$p_value <= 1)
  expect_true(result$statistic > 0)
})

test_that("spurtest i1resid works with covariates", {
  set.seed(123)
  n <- 30; q <- 10; nrep <- 1000

  y <- rnorm(n)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  coords <- matrix(rnorm(n * 2), ncol = 2)
  dat <- data.frame(y = y, x1 = x1, x2 = x2,
                     s_1 = coords[,1], s_2 = coords[,2])

  result <- spurtest(y ~ x1 + x2, data = dat, coords = ~ s_1 + s_2,
                     type = "i1resid", latlong = FALSE, q = q, nrep = nrep)

  expect_s3_class(result, c("spur_test_i1resid", "spur_test"))
  expect_true(result$p_value >= 0 && result$p_value <= 1)
})

test_that("print and summary methods work", {
  set.seed(123)
  n <- 20; q <- 8

  y <- rnorm(n)
  coords <- matrix(rnorm(n * 2), ncol = 2)
  dat <- data.frame(y = y, s_1 = coords[,1], s_2 = coords[,2])

  result <- spurtest(y ~ 1, data = dat, coords = ~ s_1 + s_2,
                     type = "i1", latlong = FALSE, q = q, nrep = 500)

  expect_output(print(result), "Test Results")
  expect_output(summary(result), "Test Summary")
})

test_that("spurtest on Chetty data gives reasonable results", {
  data(chetty, package = "spatialunitroot")

  set.seed(42)
  t1 <- spurtest(am ~ 1, data = chetty, coords = ~ s_1 + s_2,
                 type = "i1", latlong = TRUE, q = 15, nrep = 2000)

  # AM is known to be spatially persistent
  # I(1) test should NOT reject (high p-value)
  expect_true(t1$p_value > 0.05)

  set.seed(42)
  t2 <- spurtest(am ~ 1, data = chetty, coords = ~ s_1 + s_2,
                 type = "i0", latlong = TRUE, q = 15, nrep = 2000)

  # I(0) test should REJECT (low p-value)
  expect_true(t2$p_value < 0.05)
})
