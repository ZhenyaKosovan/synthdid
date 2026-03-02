# Tests for staggered synthetic difference-in-differences

# =============================================================================
# Helper: generate a staggered panel with known ATT
# =============================================================================
make_staggered_panel <- function(N_control = 20, N_treat_per_cohort = 5,
                                 cohort_times = c(11, 15), TT = 20,
                                 tau = 2, sigma = 0.5, seed = 42) {
  set.seed(seed)
  N_treat <- length(cohort_times) * N_treat_per_cohort
  N <- N_control + N_treat
  L <- outer(1:N, 1:TT, function(i, t) sin(i / 5) + t / 10)
  Y <- L + matrix(rnorm(N * TT, sd = sigma), N, TT)
  W <- matrix(0, N, TT)

  idx <- N_control
  for (g in cohort_times) {
    units <- (idx + 1):(idx + N_treat_per_cohort)
    W[units, g:TT] <- 1
    Y[units, g:TT] <- Y[units, g:TT] + tau
    idx <- idx + N_treat_per_cohort
  }

  rownames(Y) <- rownames(W) <- paste0("unit", 1:N)
  colnames(Y) <- colnames(W) <- 1:TT
  list(Y = Y, W = W, tau = tau, N_control = N_control,
       cohort_times = cohort_times)
}


# Helper: convert matrix panel to long format
panel_to_long <- function(Y, W) {
  N <- nrow(Y)
  TT <- ncol(Y)
  data.frame(
    unit = rep(rownames(Y), each = TT),
    time = rep(as.integer(colnames(Y)), N),
    outcome = c(t(Y)),
    treatment = c(t(W)),
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# Phase 1: Design layer tests
# =============================================================================

test_that("compute_first_treat_time works correctly", {
  W <- matrix(0, 5, 10)
  W[3, 6:10] <- 1
  W[4, 8:10] <- 1
  rownames(W) <- paste0("u", 1:5)
  ft <- compute_first_treat_time(W)
  expect_equal(ft[["u1"]], Inf)
  expect_equal(ft[["u3"]], 6)
  expect_equal(ft[["u4"]], 8)
  expect_length(ft, 5)
})


test_that("validate_staggered_panel catches invalid inputs", {
  W_good <- matrix(0, 5, 10)
  W_good[4, 5:10] <- 1
  W_good[5, 8:10] <- 1
  Y <- matrix(rnorm(50), 5, 10)

  # Should pass

  expect_true(validate_staggered_panel(Y, W_good))

  # Non-binary
  W_bad <- W_good
  W_bad[1, 1] <- 0.5
  expect_error(validate_staggered_panel(Y, W_bad), "binary")

  # Treatment in first period
  W_bad2 <- W_good
  W_bad2[3, 1] <- 1
  expect_error(validate_staggered_panel(Y, W_bad2), "first period")

  # Non-absorbing
  W_bad3 <- W_good
  W_bad3[4, 7] <- 0  # gap in treatment
  expect_error(validate_staggered_panel(Y, W_bad3), "absorbing")

  # Dimension mismatch
  expect_error(validate_staggered_panel(Y, W_good[1:3, ]), "dimensions")
})


test_that("identify_cohorts groups correctly", {
  ft <- c(u1 = Inf, u2 = 5, u3 = 5, u4 = 8, u5 = Inf)
  cohorts <- identify_cohorts(ft)
  expect_equal(names(cohorts), c("5", "8"))
  expect_equal(unname(cohorts[["5"]]), c(2, 3))
  expect_equal(unname(cohorts[["8"]]), 4)
})


test_that("build_staggered_subproblem constructs valid subproblems", {
  panel <- make_staggered_panel()
  Y <- panel$Y
  W <- panel$W
  ft <- compute_first_treat_time(W)

  sub <- build_staggered_subproblem(Y, W, ft, cohort_time = 11)
  expect_equal(sub$T0, 10)
  # With max_post=NULL, post_end=20. Cohort 15 adopts at 15 <= 20,

  # so only 20 never-treated are eligible as controls.
  expect_equal(sub$N0, 20)
  expect_equal(nrow(sub$Y), 25)  # 20 control + 5 treated
  expect_equal(ncol(sub$Y), 20)  # all 20 periods

  # With max_post limited, cohort 15 becomes not-yet-treated
  sub_limited <- build_staggered_subproblem(Y, W, ft, cohort_time = 11,
                                            max_post = 3)
  # post_end = 13, cohort 15 adopts at 15 > 13, so 25 controls
  expect_equal(sub_limited$N0, 25)

  # Second cohort: controls are only never-treated (cohort 11 already treated)
  sub2 <- build_staggered_subproblem(Y, W, ft, cohort_time = 15)
  expect_equal(sub2$T0, 14)
  expect_equal(sub2$N0, 20)  # only never-treated
})


test_that("build_staggered_subproblem respects control_type", {
  panel <- make_staggered_panel()
  ft <- compute_first_treat_time(panel$W)

  # never_treated: only units that are never treated
  sub_never <- build_staggered_subproblem(panel$Y, panel$W, ft, 11,
                                          control_type = "never_treated")
  expect_equal(sub_never$N0, 20)

  # not_yet_treated: includes units treated later
  sub_nyt <- build_staggered_subproblem(panel$Y, panel$W, ft, 11,
                                        control_type = "not_yet_treated")
  expect_true(sub_nyt$N0 >= sub_never$N0)
})


test_that("enumerate_staggered_subproblems builds all cohorts", {
  panel <- make_staggered_panel()
  subs <- enumerate_staggered_subproblems(panel$Y, panel$W)
  expect_equal(names(subs), c("11", "15"))
  expect_true(all(vapply(subs, function(s) s$T0 >= 2, logical(1))))
})


# =============================================================================
# Phase 2: Estimator engine tests
# =============================================================================

test_that("synthdid_staggered_estimate returns correct structure", {
  panel <- make_staggered_panel()
  est <- synthdid_staggered_estimate(panel$Y, panel$W)

  expect_s3_class(est, "synthdid_staggered")
  expect_true(is.numeric(est$att))
  expect_true(is.data.frame(est$cohort_effects))
  expect_equal(nrow(est$cohort_effects), 2)
  expect_equal(sum(est$aggregation_weights), 1)
  expect_equal(length(est$subproblem_fits), 2)
})


test_that("staggered estimate recovers known ATT", {
  panel <- make_staggered_panel(tau = 3, sigma = 0.3, seed = 123)
  est <- synthdid_staggered_estimate(panel$Y, panel$W)

  # Should be within reasonable tolerance of true tau = 3
  expect_true(abs(est$att - 3) < 1.0,
              info = sprintf("ATT = %.3f, expected ~3", est$att))
})


test_that("staggered with method = 'did' works", {
  panel <- make_staggered_panel()
  est <- synthdid_staggered_estimate(panel$Y, panel$W, method = "did")
  expect_s3_class(est, "synthdid_staggered")
  expect_equal(est$method, "did")
})


test_that("staggered with method = 'sc' works", {
  panel <- make_staggered_panel()
  est <- synthdid_staggered_estimate(panel$Y, panel$W, method = "sc")
  expect_s3_class(est, "synthdid_staggered")
  expect_equal(est$method, "sc")
})


test_that("aggregation weights work correctly", {
  panel <- make_staggered_panel()

  est_ut <- synthdid_staggered_estimate(panel$Y, panel$W,
                                         weights = "treated_unit_time")
  est_cs <- synthdid_staggered_estimate(panel$Y, panel$W,
                                         weights = "cohort_size")
  est_eq <- synthdid_staggered_estimate(panel$Y, panel$W,
                                         weights = "equal")

  # equal should give 0.5/0.5
  expect_equal(unname(est_eq$aggregation_weights), c(0.5, 0.5))

  # cohort_size: 5 and 5, so also equal
  expect_equal(unname(est_cs$aggregation_weights), c(0.5, 0.5))

  # treated_unit_time: 5*10 vs 5*6 = 50 vs 30
  expect_equal(unname(est_ut$aggregation_weights),
               c(50, 30) / 80, tolerance = 1e-10)
})


# =============================================================================
# Equivalence: staggered with one cohort == simultaneous
# =============================================================================

test_that("staggered with single cohort matches simultaneous estimator", {
  set.seed(99)
  N <- 25; TT <- 15
  Y <- matrix(rnorm(N * TT), N, TT)
  W <- matrix(0, N, TT)
  W[21:25, 8:TT] <- 1
  Y[21:25, 8:TT] <- Y[21:25, 8:TT] + 2
  rownames(Y) <- paste0("u", 1:N)
  colnames(Y) <- 1:TT

  # Staggered estimate
  est_stag <- synthdid_staggered_estimate(Y, W)

  # Simultaneous estimate (units reordered: controls first)
  unit_order <- order(W[, 8])
  suppressWarnings({
    est_sim <- synthdid_estimate(Y[unit_order, ], N0 = 20, T0 = 7)
  })

  # Should produce very similar results (same data, same method)
  expect_equal(est_stag$att, c(est_sim), tolerance = 0.1,
               info = sprintf("stag=%.3f, sim=%.3f", est_stag$att, c(est_sim)))
})


# =============================================================================
# Formula interface integration
# =============================================================================

test_that("synthdid() auto-detects staggered adoption", {
  panel <- make_staggered_panel()
  df <- panel_to_long(panel$Y, panel$W)

  est <- synthdid(outcome ~ treatment, data = df, index = c("unit", "time"))
  expect_s3_class(est, "synthdid_staggered")
  expect_true(abs(est$att - 2) < 1.0)
})


test_that("synthdid() simultaneous adoption still works", {
  data(california_prop99)
  est <- synthdid(PacksPerCapita ~ treated, data = california_prop99,
                  index = c("State", "Year"))
  expect_s3_class(est, "synthdid")
  expect_true(abs(coef(est) - (-15.6)) < 1)
})


test_that("synthdid() errors when forcing simultaneous on staggered data", {
  panel <- make_staggered_panel()
  df <- panel_to_long(panel$Y, panel$W)

  expect_error(
    synthdid(outcome ~ treatment, data = df, index = c("unit", "time"),
             adoption = "simultaneous"),
    "staggered"
  )
})


test_that("synthdid() respects adoption = 'staggered' explicitly", {
  panel <- make_staggered_panel()
  df <- panel_to_long(panel$Y, panel$W)

  est <- synthdid(outcome ~ treatment, data = df, index = c("unit", "time"),
                  adoption = "staggered")
  expect_s3_class(est, "synthdid_staggered")
})


# =============================================================================
# S3 methods for synthdid_staggered
# =============================================================================

test_that("coef.synthdid_staggered works", {
  panel <- make_staggered_panel()
  est <- synthdid_staggered_estimate(panel$Y, panel$W)

  agg <- coef(est)
  expect_named(agg, "att")

  cohort <- coef(est, type = "cohort")
  expect_length(cohort, 2)
  expect_true(all(grepl("^cohort_", names(cohort))))
})


test_that("confint.synthdid_staggered warns without SE", {
  panel <- make_staggered_panel()
  est <- synthdid_staggered_estimate(panel$Y, panel$W)

  expect_warning(confint(est), "Standard error not available")
})


test_that("predict.synthdid_staggered returns expected types", {
  panel <- make_staggered_panel()
  est <- synthdid_staggered_estimate(panel$Y, panel$W)

  expect_type(predict(est, type = "aggregate"), "double")
  expect_s3_class(predict(est, type = "cohort"), "data.frame")
  ev <- predict(est, type = "event")
  expect_s3_class(ev, "data.frame")
  expect_true(all(c("cohort_time", "relative_time", "effect") %in% names(ev)))
})


test_that("print and summary work without error", {
  panel <- make_staggered_panel()
  est <- synthdid_staggered_estimate(panel$Y, panel$W)

  expect_output(print(est), "Staggered")
  expect_output(suppressWarnings(print(summary(est))), "ATT")
})


test_that("plot.synthdid_staggered creates ggplot objects", {
  skip_if_not_installed("ggplot2")
  panel <- make_staggered_panel()
  est <- synthdid_staggered_estimate(panel$Y, panel$W)

  # New API: type = "effect" with subtype
  p1 <- plot(est, type = "effect", subtype = "cohort")
  expect_s3_class(p1, "gg")

  p2 <- plot(est, type = "effect", subtype = "event")
  expect_s3_class(p2, "gg")

  # Legacy API: type = "cohort" / "event" still works with deprecation warning
  lifecycle::expect_deprecated(p3 <- plot(est, type = "cohort"))
  expect_s3_class(p3, "gg")

  lifecycle::expect_deprecated(p4 <- plot(est, type = "event"))
  expect_s3_class(p4, "gg")
})


# =============================================================================
# Error handling
# =============================================================================

test_that("staggered estimation errors on non-absorbing treatment", {
  W <- matrix(0, 10, 10)
  W[8, 5:8] <- 1  # gap at period 9-10
  Y <- matrix(rnorm(100), 10, 10)
  rownames(Y) <- rownames(W) <- paste0("u", 1:10)
  colnames(Y) <- colnames(W) <- 1:10

  expect_error(synthdid_staggered_estimate(Y, W), "absorbing")
})


test_that("staggered estimation errors on no treated units", {
  Y <- matrix(rnorm(100), 10, 10)
  W <- matrix(0, 10, 10)
  rownames(Y) <- rownames(W) <- paste0("u", 1:10)
  colnames(Y) <- colnames(W) <- 1:10

  expect_error(synthdid_staggered_estimate(Y, W), "No treated")
})


test_that("staggered estimation errors on treatment in first period", {
  Y <- matrix(rnorm(100), 10, 10)
  W <- matrix(0, 10, 10)
  W[8, 1:10] <- 1
  rownames(Y) <- rownames(W) <- paste0("u", 1:10)
  colnames(Y) <- colnames(W) <- 1:10

  expect_error(synthdid_staggered_estimate(Y, W), "first period")
})
