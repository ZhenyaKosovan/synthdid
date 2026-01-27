test_that("timesteps falls back to labels when not parseable", {
  Y <- matrix(1:4, nrow = 2)
  colnames(Y) <- c("foo", "bar")
  expect_equal(timesteps(Y), c("foo", "bar"))

  Y2 <- matrix(1:4, nrow = 2)
  colnames(Y2) <- c("2020-01-01", "2020-01-02")
  expect_equal(timesteps(Y2), as.Date(c("2020-01-01", "2020-01-02")))
})

test_that("panel.matrices errors when treatment starts in first period", {
  panel <- data.frame(
    unit = rep(c("A", "B"), each = 2),
    time = rep(1:2, times = 2),
    outcome = c(1, 2, 3, 4),
    treated = c(0, 0, 1, 1)
  )
  expect_error(
    panel.matrices(panel),
    "Treatment starts in the first period"
  )
})

test_that("synthdid_estimate errors out on single pre-treatment period", {
  Y <- matrix(c(1, 2, 1, 2), nrow = 2, byrow = TRUE)
  expect_error({
    est <- synthdid_estimate(Y, N0 = 1, T0 = 1, sparsify = NULL, max.iter = 5)
  })
})

test_that("fit_ar2 errors on singular system", {
  E <- matrix(0, nrow = 2, ncol = 3)
  expect_error(synthdid:::fit_ar2(E), "singular")
})
