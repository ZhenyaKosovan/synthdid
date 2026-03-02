# Tests for the new plot engine (R/plot_engine.R)

# =============================================================================
# Helper: create a basic estimate for testing
# =============================================================================
make_test_estimate <- function() {
  data(california_prop99, package = "synthdid")
  setup <- withCallingHandlers(
    panel.matrices(california_prop99),
    lifecycle_warning_deprecated = function(cnd) invokeRestart("muffleWarning")
  )
  withCallingHandlers(
    synthdid_estimate(setup$Y, setup$N0, setup$T0),
    lifecycle_warning_deprecated = function(cnd) invokeRestart("muffleWarning")
  )
}

# Helper: create a staggered estimate for testing
make_test_staggered <- function() {
  set.seed(42)
  N_control <- 20
  N_treat_per_cohort <- 5
  cohort_times <- c(11, 15)
  TT <- 20
  tau <- 2
  sigma <- 0.5

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

  long_df <- data.frame(
    unit = rep(rownames(Y), each = TT),
    time = rep(as.integer(colnames(Y)), N),
    outcome = c(t(Y)),
    treatment = c(t(W)),
    stringsAsFactors = FALSE
  )

  synthdid(outcome ~ treatment, data = long_df, index = c("unit", "time"),
           adoption = "staggered")
}


# =============================================================================
# Layer 1: Payload extraction tests
# =============================================================================

test_that("extract_plot_payload.synthdid_estimate returns correct structure", {
  est <- make_test_estimate()
  payload <- extract_plot_payload.synthdid_estimate(est)

  expect_type(payload, "list")
  expect_named(payload, c(
    "treated_trajectory", "synthetic_trajectory", "effect_curve",
    "omega", "lambda", "time", "time_labels", "T0", "T1",
    "N0", "N1", "tau", "se", "estimator", "unit_names"
  ), ignore.order = TRUE)

  # Dimensions match setup
  setup <- attr(est, "setup")
  expect_length(payload$treated_trajectory, ncol(setup$Y))
  expect_length(payload$synthetic_trajectory, ncol(setup$Y))
  expect_length(payload$effect_curve, ncol(setup$Y) - setup$T0)
  expect_length(payload$omega, setup$N0)
  expect_length(payload$lambda, setup$T0)
  expect_length(payload$time, ncol(setup$Y))
  expect_equal(payload$T0, setup$T0)
  expect_equal(payload$N0, setup$N0)
})

test_that("extract_plot_payload effect_curve matches synthdid_effect_curve", {
  est <- make_test_estimate()
  payload <- extract_plot_payload.synthdid_estimate(est)
  ref_effect <- withCallingHandlers(
    synthdid_effect_curve(est),
    lifecycle_warning_deprecated = function(cnd) invokeRestart("muffleWarning")
  )

  expect_equal(as.numeric(payload$effect_curve), as.numeric(ref_effect),
               tolerance = 1e-10)
})

test_that("extract_plot_payload with include_controls adds control trajectories", {
  est <- make_test_estimate()

  payload_no <- extract_plot_payload.synthdid_estimate(est, include_controls = FALSE)
  expect_null(payload_no$control_trajectories)

  payload_yes <- extract_plot_payload.synthdid_estimate(est, include_controls = TRUE)
  expect_false(is.null(payload_yes$control_trajectories))
  expect_equal(nrow(payload_yes$control_trajectories), payload_yes$N0)
  expect_equal(ncol(payload_yes$control_trajectories), length(payload_yes$time))
})

test_that("extract_plot_payload.synthdid_staggered returns correct structure for effect", {
  stag <- make_test_staggered()
  payload <- extract_plot_payload.synthdid_staggered(stag, type = "effect")

  expect_type(payload, "list")
  expect_equal(payload$type, "effect")
  expect_false(is.null(payload$cohort_effects))
  expect_false(is.null(payload$event_data))
  expect_equal(payload$att, stag$att)
})

test_that("extract_plot_payload.synthdid_staggered trajectory returns per-cohort payloads", {
  stag <- make_test_staggered()
  payload <- extract_plot_payload.synthdid_staggered(stag, type = "trajectory")

  expect_equal(payload$type, "trajectory")
  expect_length(payload$cohort_payloads, length(stag$subproblem_fits))
  for (cp in payload$cohort_payloads) {
    expect_true("treated_trajectory" %in% names(cp))
    expect_true("synthetic_trajectory" %in% names(cp))
    expect_true("cohort_time" %in% names(cp))
  }
})


# =============================================================================
# Layer 2: Spec builder tests
# =============================================================================

test_that("build_plot_spec resolves auto mode based on N0", {
  est <- make_test_estimate()
  payload <- extract_plot_payload.synthdid_estimate(est)

  # California has N0 = 38, which is < 100 threshold
  spec <- build_plot_spec(payload, type = "trajectory", mode = "auto")
  expect_equal(spec$mode, "full")
})

test_that("build_plot_spec resolves auto to top_k for large N0", {
  est <- make_test_estimate()
  payload <- extract_plot_payload.synthdid_estimate(est)
  # Override N0 to simulate large panel
  payload$N0 <- 200

  spec <- build_plot_spec(payload, type = "trajectory", mode = "auto")
  expect_equal(spec$mode, "top_k")
})

test_that("build_plot_spec trajectory creates expected data frames", {
  est <- make_test_estimate()
  payload <- extract_plot_payload.synthdid_estimate(est, include_controls = TRUE)
  spec <- build_plot_spec(payload, type = "trajectory")

  expect_true("trajectories" %in% names(spec$data))
  expect_true("lambda" %in% names(spec$data))
  expect_true("controls" %in% names(spec$data))

  traj <- spec$data$trajectories
  expect_true(all(c("time", "value", "series") %in% names(traj)))
  expect_equal(sort(unique(traj$series)), c("Synthetic Control", "Treated"))
})

test_that("build_plot_spec weights creates expected data frames", {
  est <- make_test_estimate()
  payload <- extract_plot_payload.synthdid_estimate(est)
  spec <- build_plot_spec(payload, type = "weights")

  expect_true("weights" %in% names(spec$data))
  w <- spec$data$weights
  expect_true(all(c("unit", "weight") %in% names(w)))
  # All weights non-negative
  expect_true(all(w$weight >= 0))
})

test_that("build_plot_spec effect creates expected data frames", {
  est <- make_test_estimate()
  payload <- extract_plot_payload.synthdid_estimate(est)
  spec <- build_plot_spec(payload, type = "effect")

  expect_true("effect" %in% names(spec$data))
  e <- spec$data$effect
  expect_true(all(c("time", "effect") %in% names(e)))
  expect_equal(nrow(e), payload$T1)
})


# =============================================================================
# Layer 3: Renderer tests
# =============================================================================

test_that("render_trajectory_plot returns a ggplot", {
  est <- make_test_estimate()
  payload <- extract_plot_payload.synthdid_estimate(est, include_controls = TRUE)
  spec <- build_plot_spec(payload, type = "trajectory")

  p <- render_trajectory_plot(spec)
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
})

test_that("render_weights_plot returns a ggplot", {
  est <- make_test_estimate()
  payload <- extract_plot_payload.synthdid_estimate(est)
  spec <- build_plot_spec(payload, type = "weights")

  p <- render_weights_plot(spec)
  expect_s3_class(p, "gg")
})

test_that("render_effect_plot returns a ggplot for simultaneous estimate", {
  est <- make_test_estimate()
  payload <- extract_plot_payload.synthdid_estimate(est)
  spec <- build_plot_spec(payload, type = "effect")

  p <- render_effect_plot(spec)
  expect_s3_class(p, "gg")
})

test_that("render_effect_plot returns ggplot for staggered cohort", {
  stag <- make_test_staggered()
  payload <- extract_plot_payload.synthdid_staggered(stag, type = "effect")
  spec <- build_plot_spec(payload, type = "effect")

  p <- render_effect_plot(spec, subtype = "cohort")
  expect_s3_class(p, "gg")
})

test_that("render_effect_plot returns ggplot for staggered event", {
  stag <- make_test_staggered()
  payload <- extract_plot_payload.synthdid_staggered(stag, type = "effect")
  spec <- build_plot_spec(payload, type = "effect")

  p <- render_effect_plot(spec, subtype = "event")
  expect_s3_class(p, "gg")
})


# =============================================================================
# S3 dispatch / backward compatibility tests
# =============================================================================

test_that("plot.synthdid_estimate dispatches to new engine by default", {
  est <- make_test_estimate()
  p <- plot(est)
  expect_s3_class(p, "gg")
})

test_that("plot.synthdid_estimate type='weights' works", {
  est <- make_test_estimate()
  p <- plot(est, type = "weights")
  expect_s3_class(p, "gg")
})

test_that("plot.synthdid_estimate type='effect' works", {
  est <- make_test_estimate()
  p <- plot(est, type = "effect")
  expect_s3_class(p, "gg")
})

test_that("plot.synthdid_estimate delegates to legacy when legacy args present", {
  est <- make_test_estimate()
  # spaghetti.units is a legacy arg, should delegate to synthdid_plot
  p <- plot(est, spaghetti.units = c())
  expect_s3_class(p, "gg")
})

test_that("plot.synthdid_staggered default works", {
  stag <- make_test_staggered()
  p <- plot(stag)
  expect_s3_class(p, "gg")
})

test_that("plot.synthdid_staggered legacy type='cohort' still works", {
  stag <- make_test_staggered()
  # Should produce a deprecation warning but still work
  lifecycle::expect_deprecated(
    p <- plot(stag, type = "cohort")
  )
  expect_s3_class(p, "gg")
})

test_that("plot.synthdid_staggered legacy type='event' still works", {
  stag <- make_test_staggered()
  lifecycle::expect_deprecated(
    p <- plot(stag, type = "event")
  )
  expect_s3_class(p, "gg")
})

test_that("plot.synthdid_staggered type='trajectory' works", {
  stag <- make_test_staggered()
  p <- plot(stag, type = "trajectory")
  expect_s3_class(p, "gg")
})

test_that("plot.synthdid_staggered type='weights' works", {
  stag <- make_test_staggered()
  p <- plot(stag, type = "weights")
  expect_s3_class(p, "gg")
})

test_that("trajectory plot has expected number of layers", {
  est <- make_test_estimate()
  p <- plot(est, type = "trajectory")
  # At minimum: geom_line (controls), geom_line (trajectories),
  # geom_vline (onset), geom_ribbon (lambda) = 4 layers
  # But controls only if include_controls; minimum is 2 (lines + vline)
  expect_gte(length(p$layers), 2)
})
