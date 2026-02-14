# Error Condition Tests for synthdid package
# Addresses SUGGESTIONS.md Issue #8: Insufficient tests for error conditions

test_that("synthdid_estimate errors on non-numeric Y", {
  Y <- matrix(c("a", "b", "c", "d"), 2, 2)

  expect_error({
    expect_warning({
      synthdid_estimate(Y, N0 = 1, T0 = 1)
    })
  })
})

test_that("synthdid_estimate errors on non-matrix Y", {
  Y <- c(1, 2, 3, 4)

  expect_error({
    synthdid_estimate(Y, N0 = 1, T0 = 1)
  })
})

test_that("synthdid_estimate errors on N0 >= nrow(Y)", {
  Y <- matrix(1:10, 5, 2)

  expect_error({
    synthdid_estimate(Y, N0 = 5, T0 = 1)
  })
})

test_that("synthdid_estimate errors on T0 >= ncol(Y)", {
  Y <- matrix(1:10, 2, 5)

  expect_error({
    synthdid_estimate(Y, N0 = 1, T0 = 5)
  })
})

test_that("synthdid_estimate errors on negative N0", {
  Y <- matrix(1:10, 5, 2)

  expect_error({
    synthdid_estimate(Y, N0 = -1, T0 = 1)
  })
})

test_that("synthdid_estimate errors on negative T0", {
  Y <- matrix(1:10, 2, 5)

  expect_error({
    synthdid_estimate(Y, N0 = 1, T0 = -1)
  })
})

test_that("synthdid_estimate errors on zero N0", {
  Y <- matrix(1:10, 5, 2)

  expect_error({
    synthdid_estimate(Y, N0 = 0, T0 = 1)
  })
})

test_that("synthdid_estimate errors on zero T0", {
  Y <- matrix(1:10, 2, 5)

  expect_error({
    synthdid_estimate(Y, N0 = 1, T0 = 0)
  })
})

test_that("synthdid_estimate handles NA values", {
  Y <- matrix(1:10, 2, 5)
  Y[1, 1] <- NA

  expect_error({
    synthdid_estimate(Y, N0 = 1, T0 = 2)
  })
})

test_that("synthdid_estimate handles Inf values", {
  Y <- matrix(1:10, 2, 5)
  Y[1, 1] <- Inf

  # Should error
  expect_error({
    synthdid_estimate(Y, N0 = 1, T0 = 2)
  })
})

test_that("synthdid_estimate errors on mismatched covariate dimensions", {
  Y <- matrix(rnorm(20), 4, 5)
  X <- array(rnorm(15), dim = c(3, 5, 1)) # Wrong number of rows

  expect_error({
    synthdid_estimate(Y, N0 = 2, T0 = 3, X = X)
  })
})

test_that("synthdid_estimate errors on mismatched covariate periods", {
  Y <- matrix(rnorm(20), 4, 5)
  X <- array(rnorm(16), dim = c(4, 4, 1)) # Wrong number of columns

  expect_error({
    synthdid_estimate(Y, N0 = 2, T0 = 3, X = X)
  })
})


test_that("panel.matrices errors on non-data.frame", {
  panel <- matrix(1:20, 5, 4)

  expect_error({
    panel.matrices(panel)
  })
})

test_that("panel.matrices errors on unbalanced panel (missing observation)", {
  data(california_prop99)
  panel_unbalanced <- california_prop99[-10, ] # Remove one observation

  expect_error({
    panel.matrices(panel_unbalanced)
  })
})

test_that("panel.matrices errors on non-simultaneous treatment", {
  data(california_prop99)
  panel_mod <- california_prop99

  # Make Kansas treated one year earlier than California
  panel_mod[panel_mod$State == "Kansas" & panel_mod$Year >= 1988, "treated"] <- 1L

  expect_error({
    panel.matrices(panel_mod)
  })
})

test_that("panel.matrices errors when treatment starts in first period", {
  panel <- data.frame(
    unit = rep(c("A", "B"), each = 3),
    time = rep(1:3, times = 2),
    outcome = rnorm(6),
    treated = c(1, 1, 1, 0, 0, 0) # Treatment starts in period 1
  )

  expect_error(
    {
      panel.matrices(panel)
    },
    "first period"
  )
})

test_that("vcov errors on unknown method", {
  set.seed(123)
  Y <- matrix(rnorm(30), 5, 6)
  Y[4:5, 4:6] <- Y[4:5, 4:6] + 2
  estimate <- synthdid_estimate(Y, N0 = 3, T0 = 3)

  expect_error({
    vcov(estimate, method = "unknown_method")
  })
})

test_that("placebo_se errors when N0 <= N1", {
  set.seed(123)
  Y <- matrix(rnorm(30), 5, 6)
  Y[4:5, 4:6] <- Y[4:5, 4:6] + 2
  estimate <- synthdid_estimate(Y, N0 = 2, T0 = 3) # N0=2, N1=3

  expect_error(
    {
      vcov(estimate, method = "placebo", replications = 10)
    },
    "more controls than treated"
  )
})

test_that("sc_weight_fw errors on non-matrix Y", {
  Y <- c(1, 2, 3, 4)

  expect_error({
    synthdid:::sc_weight_fw(Y, zeta = 1, intercept = TRUE)
  })
})

test_that("contract3_cpp errors on non-3D array", {
  X <- matrix(1:12, 3, 4)
  v <- c(0.5, 0.5)

  expect_error(
    {
      contract3_cpp(X, v)
    },
    "3D array"
  )
})

test_that("contract3_cpp errors on mismatched vector length", {
  X <- array(1:24, dim = c(3, 4, 2))
  v <- c(0.5, 0.3, 0.2) # Length 3, should be 2

  expect_error(
    {
      contract3_cpp(X, v)
    },
    "Length of v must match third dimension"
  )
})

test_that("synthdid_controls errors on non-estimate input", {
  expect_error({
    synthdid_controls("not an estimate")
  })
})


test_that("negative noise.level parameter errors", {
  Y <- matrix(rnorm(20), 4, 5)

  expect_error({
    synthdid_estimate(Y, N0 = 2, T0 = 3, noise.level = -1)
  })
})

test_that("negative max.iter errors", {
  Y <- matrix(rnorm(20), 4, 5)

  # Should handle gracefully or error
  expect_error({
    synthdid_estimate(Y, N0 = 2, T0 = 3, max.iter = -100)
  })
})


test_that("formula interface errors on missing treatment variable", {
  data(california_prop99)
  df <- california_prop99[, c("State", "Year", "PacksPerCapita")] # No 'treated'

  expect_error(
    {
      synthdid(PacksPerCapita ~ treated,
        data = df,
        index = c("State", "Year")
      )
    },
    "not found"
  )
})

test_that("formula interface errors on non-binary treatment", {
  data(california_prop99)
  df <- california_prop99
  df$treated <- df$treated + 0.5 # Make it non-binary
  expect_error({
    synthdid(PacksPerCapita ~ treated,
      data = df,
      index = c("State", "Year")
    )
  })
})

test_that("synthdid_memory_estimate errors on negative inputs", {
  expect_error({
    synthdid_memory_estimate(N = -10, T = 50)
  })

  expect_error({
    synthdid_memory_estimate(N = 10, T = -50)
  })

  expect_error({
    synthdid_memory_estimate(N = 10, T = 50, K = -5)
  })

  expect_error({
    synthdid_memory_estimate(N = 10, T = 50, replications = -100)
  })
})

test_that("synthdid_convergence_info handles missing convergence attribute", {
  # Create estimate without running through normal flow
  Y <- matrix(rnorm(20), 4, 5)
  estimate <- synthdid_estimate(Y, N0 = 2, T0 = 3)

  # Remove convergence attribute
  attr(estimate, "convergence") <- NULL

  expect_message(
    {
      info <- synthdid_convergence_info(estimate)
    },
    "does not have convergence information"
  )

  expect_null(info)
})

test_that("confint errors when SE not available", {
  Y <- matrix(rnorm(20), 4, 5)
  estimate <- synthdid_estimate(Y, N0 = 2, T0 = 3)

  # No SE computed
  attr(estimate, "se") <- NA
  expect_error(confint(estimate))
})

test_that("sparsify function errors on bad input", {
  expect_error(synthdid:::sparsify_function("not a number"))
})

test_that("errors are informative", {
  Y <- matrix(1:10, 5, 2)

  # Check that error messages are helpful
  expect_error(
    {
      synthdid_estimate(Y, N0 = 10, T0 = 1)
    },
    ".*"
  ) # Should contain some message

  expect_error(
    {
      synthdid_estimate(Y, N0 = 1, T0 = 10)
    },
    ".*"
  ) # Should contain some message
})
