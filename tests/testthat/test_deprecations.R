# Tests for lifecycle deprecation warnings on old API functions

test_that("synthdid_estimate() is soft-deprecated", {
  withr::local_options(lifecycle_verbosity = "warning")
  setup <- random.low.rank()
  expect_snapshot({
    tau.hat <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
  })
  expect_true(is.finite(c(tau.hat)))
})

test_that("sc_estimate() is soft-deprecated", {
  withr::local_options(lifecycle_verbosity = "warning")
  setup <- random.low.rank()
  expect_snapshot({
    tau.hat <- sc_estimate(setup$Y, setup$N0, setup$T0)
  })
  expect_true(is.finite(c(tau.hat)))
})

test_that("did_estimate() is soft-deprecated", {
  withr::local_options(lifecycle_verbosity = "warning")
  setup <- random.low.rank()
  expect_snapshot({
    tau.hat <- did_estimate(setup$Y, setup$N0, setup$T0)
  })
  expect_true(is.finite(c(tau.hat)))
})

test_that("synthdid_se() is soft-deprecated", {
  withr::local_options(lifecycle_verbosity = "warning")
  setup <- random.low.rank()
  tau.hat <- withr::with_options(
    list(lifecycle_verbosity = "quiet"),
    synthdid_estimate(setup$Y, setup$N0, setup$T0)
  )
  expect_snapshot({
    se <- synthdid_se(tau.hat, method = "jackknife")
  })
})

test_that("panel.matrices() is soft-deprecated", {
  withr::local_options(lifecycle_verbosity = "warning")
  data(california_prop99)
  expect_snapshot({
    setup <- panel.matrices(california_prop99)
  })
  expect_true(is.list(setup))
})

test_that("synthdid_effect_curve() is soft-deprecated", {
  withr::local_options(lifecycle_verbosity = "warning")
  setup <- random.low.rank()
  tau.hat <- withr::with_options(
    list(lifecycle_verbosity = "quiet"),
    synthdid_estimate(setup$Y, setup$N0, setup$T0)
  )
  expect_snapshot({
    curve <- synthdid_effect_curve(tau.hat)
  })
  expect_true(is.numeric(curve))
})

test_that("synthdid_converged() is soft-deprecated", {
  withr::local_options(lifecycle_verbosity = "warning")
  setup <- random.low.rank()
  tau.hat <- withr::with_options(
    list(lifecycle_verbosity = "quiet"),
    synthdid_estimate(setup$Y, setup$N0, setup$T0)
  )
  expect_snapshot({
    conv <- synthdid_converged(tau.hat)
  })
})

test_that("synthdid_convergence_info() is soft-deprecated", {
  withr::local_options(lifecycle_verbosity = "warning")
  setup <- random.low.rank()
  tau.hat <- withr::with_options(
    list(lifecycle_verbosity = "quiet"),
    synthdid_estimate(setup$Y, setup$N0, setup$T0)
  )
  expect_snapshot({
    info <- synthdid_convergence_info(tau.hat)
  })
})

test_that("synthdid_controls() is soft-deprecated", {
  withr::local_options(lifecycle_verbosity = "warning")
  setup <- random.low.rank()
  tau.hat <- withr::with_options(
    list(lifecycle_verbosity = "quiet"),
    synthdid_estimate(setup$Y, setup$N0, setup$T0)
  )
  expect_snapshot({
    ctrls <- synthdid_controls(tau.hat)
  })
})

test_that("synthdid_memory_estimate() is soft-deprecated", {
  withr::local_options(lifecycle_verbosity = "warning")
  expect_snapshot({
    mem <- synthdid_memory_estimate(N = 39, T = 31)
  })
})

test_that("timesteps() is soft-deprecated", {
  withr::local_options(lifecycle_verbosity = "warning")
  setup <- random.low.rank()
  tau.hat <- withr::with_options(
    list(lifecycle_verbosity = "quiet"),
    synthdid_estimate(setup$Y, setup$N0, setup$T0)
  )
  expect_snapshot({
    ts <- timesteps(attr(tau.hat, "setup")$Y)
  })
})

test_that("synthdid_plot() is soft-deprecated", {
  withr::local_options(lifecycle_verbosity = "warning")
  setup <- random.low.rank()
  tau.hat <- withr::with_options(
    list(lifecycle_verbosity = "quiet"),
    synthdid_estimate(setup$Y, setup$N0, setup$T0)
  )
  expect_snapshot({
    p <- synthdid_plot(tau.hat)
  })
})

test_that("synthdid_placebo() is soft-deprecated", {
  withr::local_options(lifecycle_verbosity = "warning")
  setup <- random.low.rank()
  tau.hat <- withr::with_options(
    list(lifecycle_verbosity = "quiet"),
    synthdid_estimate(setup$Y, setup$N0, setup$T0)
  )
  expect_snapshot({
    placebo <- synthdid_placebo(tau.hat)
  })
})
