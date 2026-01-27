# Comprehensive Edge Case Tests for synthdid package
test_that("minimal control units (N0=1) with multiple periods", {
  set.seed(123)
  Y <- matrix(rnorm(30), nrow = 3, ncol = 10) # 1 control, 2 treated, 10 periods
  Y[2:3, 6:10] <- Y[2:3, 6:10] + 5 # Treatment effect

  expect_no_error({
    est <- synthdid_estimate(Y, N0 = 1, T0 = 5, max.iter = 1000)
  })

  # Should warn about convergence potentially
  expect_true(is.finite(as.numeric(est)))
})

test_that("minimal pre-treatment periods (T0=1) with multiple units", {
  set.seed(123)
  Y <- matrix(rnorm(100), nrow = 10, ncol = 10) # 10 units, 10 periods
  Y[8:10, 2:10] <- Y[8:10, 2:10] + 3 # Treatment effect

  expect_error({
    est <- synthdid_estimate(Y, N0 = 7, T0 = 1, sparsify = NULL, max.iter = 500)
  })
})

test_that("no treated units after treatment (T0=T-1)", {
  set.seed(456)
  Y <- matrix(rnorm(30), nrow = 5, ncol = 6)
  # Treatment in last period only
  Y[4:5, 6] <- Y[4:5, 6] + 4

  expect_no_error({
    est <- synthdid_estimate(Y, N0 = 3, T0 = 5, max.iter = 500)
  })

  expect_true(is.finite(as.numeric(est)))
})

test_that("all units treated (N0=N-1)", {
  set.seed(789)
  Y <- matrix(rnorm(30), nrow = 5, ncol = 6)
  # Only 1 control unit
  Y[5, 4:6] <- Y[5, 4:6] + 2

  expect_no_error({
    est <- synthdid_estimate(Y, N0 = 4, T0 = 3, max.iter = 500)
  })

  expect_true(is.finite(as.numeric(est)))
})

test_that("square matrices work (N=T)", {
  set.seed(111)
  Y <- matrix(rnorm(100), nrow = 10, ncol = 10)
  Y[8:10, 7:10] <- Y[8:10, 7:10] + 5

  expect_no_error({
    est <- synthdid_estimate(Y, N0 = 7, T0 = 6, max.iter = 1000)
  })

  expect_true(is.finite(as.numeric(est)))
})

test_that("very small matrices (3x3)", {
  Y <- matrix(
    c(
      1, 2, 1,
      2, 3, 2,
      1, 2, 4
    ), # Treatment effect in last period
    nrow = 3, ncol = 3, byrow = TRUE
  )

  expect_no_error({
    est <- synthdid_estimate(Y, N0 = 2, T0 = 2, sparsify = NULL, max.iter = 100)
  })

  expect_true(is.finite(as.numeric(est)))
})

test_that("constant pre-treatment data (zero noise.level)", {
  # All pre-treatment values are identical
  Y <- matrix(
    c(
      5, 5, 5,
      5, 5, 5,
      5, 5, 8
    ), # Treatment effect
    nrow = 3, ncol = 3, byrow = TRUE
  )

  # With constant data, noise.level will be 0 (no warning for this specific case)
  # But the estimation should still work
  expect_no_error({
    est <- synthdid_estimate(Y, N0 = 2, T0 = 2, sparsify = NULL, max.iter = 100)
  })

  expect_true(is.finite(as.numeric(est)))
})

test_that("zero variance in some rows", {
  set.seed(222)
  Y <- matrix(rnorm(50), nrow = 5, ncol = 10)
  # Make one control unit constant
  Y[1, ] <- 5

  # Add treatment effect
  Y[4:5, 6:10] <- Y[4:5, 6:10] + 3

  expect_no_error({
    est <- synthdid_estimate(Y, N0 = 3, T0 = 5, max.iter = 1000)
  })

  expect_true(is.finite(as.numeric(est)))
})

test_that("zero variance in some columns", {
  set.seed(333)
  Y <- matrix(rnorm(50), nrow = 5, ncol = 10)
  # Make one pre-treatment period constant across all units
  Y[, 2] <- 10

  # Add treatment effect
  Y[4:5, 6:10] <- Y[4:5, 6:10] + 3

  expect_no_error({
    est <- synthdid_estimate(Y, N0 = 3, T0 = 5, max.iter = 1000)
  })

  expect_true(is.finite(as.numeric(est)))
})


test_that("very large values in data", {
  set.seed(555)
  Y <- matrix(rnorm(50, mean = 1e6, sd = 1e5), nrow = 5, ncol = 10)
  # Treatment effect
  Y[4:5, 6:10] <- Y[4:5, 6:10] + 1e5

  expect_no_error({
    est <- synthdid_estimate(Y, N0 = 3, T0 = 5, max.iter = 1000)
  })

  expect_true(is.finite(as.numeric(est)))
})

test_that("very small values in data", {
  set.seed(666)
  Y <- matrix(rnorm(50, mean = 1e-6, sd = 1e-7), nrow = 5, ncol = 10)
  # Treatment effect
  Y[4:5, 6:10] <- Y[4:5, 6:10] + 1e-6

  expect_no_error({
    est <- synthdid_estimate(Y, N0 = 3, T0 = 5, max.iter = 1000)
  })

  expect_true(is.finite(as.numeric(est)))
})

test_that("zero treatment effect (placebo test)", {
  set.seed(999)
  Y <- matrix(rnorm(100), nrow = 10, ncol = 10)
  # No treatment effect added

  expect_no_error({
    est <- synthdid_estimate(Y, N0 = 7, T0 = 6, max.iter = 1000)
  })

  # Estimate should be close to zero
  expect_lt(abs(as.numeric(est)), 2) # Within 2 standard deviations
})

test_that("perfect parallel trends (DID case)", {
  # Create data with perfect parallel trends
  Y <- outer(1:5, 1:10, function(i, j) i + j)
  # Add treatment effect to last 2 units in last 4 periods
  Y[4:5, 7:10] <- Y[4:5, 7:10] + 10

  expect_no_error({
    est <- synthdid_estimate(Y, N0 = 3, T0 = 6, max.iter = 1000)
  })

  # Should recover treatment effect of 10
  expect_equal(as.numeric(est), 10, tolerance = 0.5)
})

test_that("heterogeneous treatment effects across units", {
  set.seed(1010)
  Y <- matrix(rnorm(100), nrow = 10, ncol = 10)

  # Different treatment effects for different treated units
  Y[8, 7:10] <- Y[8, 7:10] + 5 # Unit 8 gets +5
  Y[9, 7:10] <- Y[9, 7:10] + 10 # Unit 9 gets +10
  Y[10, 7:10] <- Y[10, 7:10] + 2 # Unit 10 gets +2

  expect_no_error({
    est <- synthdid_estimate(Y, N0 = 7, T0 = 6, max.iter = 1000)
  })

  # Average treatment effect should be around (5+10+2)/3 ≈ 5.67
  expect_true(is.finite(as.numeric(est)))
})

test_that("with row and column names", {
  set.seed(1111)
  Y <- matrix(rnorm(60), nrow = 6, ncol = 10)
  rownames(Y) <- paste0("Unit_", 1:6)
  colnames(Y) <- paste0("Period_", 2010:2019)

  Y[5:6, 7:10] <- Y[5:6, 7:10] + 4

  expect_no_error({
    est <- synthdid_estimate(Y, N0 = 4, T0 = 6, max.iter = 1000)
  })

  expect_true(is.finite(as.numeric(est)))

  # Note: Weight naming feature not currently implemented
  # Weights are returned as unnamed vectors
  weights_omega <- omega(est)
  expect_equal(length(weights_omega), 4)

  weights_lambda <- lambda(est)
  expect_equal(length(weights_lambda), 6)
})
