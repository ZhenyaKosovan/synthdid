tol <- .03
min.decrease <- 1e-6
max.iter <- 1e6

test_that("synthdid point estimate agrees with the reference implementation", {
  skip_if_not_installed("CVXR")
  data(california_prop99)
  setup <- panel.matrices(california_prop99)
  expect_equal(
    c(synthdid_estimate(
      setup$Y,
      setup$N0,
      setup$T0,
      min.decrease = min.decrease,
      max.iter = max.iter
    )),
    c(synthdid.reference(setup$Y, setup$N0, setup$T0)),
    tolerance = tol
  )
})

test_that("sc point estimate agrees with the reference implementation", {
  skip_if_not_installed("CVXR")
  data(california_prop99)
  setup <- panel.matrices(california_prop99)
  expect_equal(
    c(sc_estimate(
      setup$Y,
      setup$N0,
      setup$T0,
      min.decrease = min.decrease,
      max.iter = max.iter
    )),
    c(sc.reference(setup$Y, setup$N0, setup$T0)),
    tolerance = tol
  )
})

test_that("did point estimate agrees with the reference implementation", {
  skip_if_not_installed("CVXR")
  data(california_prop99)
  setup <- panel.matrices(california_prop99)
  expect_equal(
    c(did_estimate(setup$Y, setup$N0, setup$T0)),
    c(did.reference(setup$Y, setup$N0, setup$T0)),
    tolerance = tol
  )
})
