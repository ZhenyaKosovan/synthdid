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

test_that("render_effect_plot event handles more than eight cohorts", {
  set.seed(123)
  event_df <- do.call(rbind, lapply(seq_len(9), function(g) {
    data.frame(
      cohort_time = g,
      relative_time = seq_len(3),
      effect = stats::rnorm(3),
      stringsAsFactors = FALSE
    )
  }))
  spec <- list(
    data = list(event_data = event_df),
    annotations = list(),
    options = list()
  )

  p <- render_effect_plot(spec, subtype = "event")
  expect_s3_class(p, "gg")
  expect_no_warning(built <- ggplot2::ggplot_build(p))
  expect_false(any(is.na(built$data[[2]]$colour)))
  expect_equal(length(unique(built$data[[2]]$group)), 9)
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

test_that("compose_staggered_trajectory_plots warns without patchwork", {
  p1 <- ggplot2::ggplot(
    data.frame(x = 1:2, y = 1:2),
    ggplot2::aes(x = x, y = y)
  ) + ggplot2::geom_line()
  p2 <- ggplot2::ggplot(
    data.frame(x = 1:2, y = 2:1),
    ggplot2::aes(x = x, y = y)
  ) + ggplot2::geom_line()

  expect_warning(
    combined <- compose_staggered_trajectory_plots(
      list(p1, p2),
      has_patchwork = FALSE
    ),
    "patchwork"
  )
  expect_identical(combined, p1)
})

test_that("trajectory plot uses linetype and color for series encoding", {
  est <- make_test_estimate()
  p <- plot(est, type = "trajectory")
  expect_false(is.null(p$scales$get_scales("colour")))
  expect_false(is.null(p$scales$get_scales("linetype")))
})

test_that("apply_time_axis_labels maps non-numeric time labels", {
  base_plot <- ggplot2::ggplot(
    data.frame(x = 1:5, y = 1:5),
    ggplot2::aes(x = x, y = y)
  ) + ggplot2::geom_line()

  p <- apply_time_axis_labels(
    base_plot,
    time_values = 1:5,
    time_labels = paste0("Q", 1:5),
    time_is_numeric = FALSE
  )

  x_scale <- p$scales$get_scales("x")
  expect_false(is.null(x_scale))
  expect_equal(x_scale$get_labels(1:5), paste0("Q", 1:5))
})

# =============================================================================
# Diagnostic plot tests
# =============================================================================

test_that("extract_diagnostic_payload returns correct structure", {
  est <- make_test_estimate()
  payload <- extract_diagnostic_payload(est)

  expect_type(payload, "list")
  expect_true("convergence_vals" %in% names(payload))
  expect_true("treated_pre" %in% names(payload))
  expect_true("synthetic_pre" %in% names(payload))
  expect_true("residuals_pre" %in% names(payload))
  expect_true("rmse_pre" %in% names(payload))
  expect_true("time_pre" %in% names(payload))

  setup <- attr(est, "setup")
  expect_length(payload$treated_pre, setup$T0)
  expect_length(payload$synthetic_pre, setup$T0)
  expect_length(payload$residuals_pre, setup$T0)
  expect_true(is.numeric(payload$rmse_pre))
  expect_true(payload$rmse_pre >= 0)
})

test_that("build_diagnostic_spec builds convergence and fit data", {
  est <- make_test_estimate()
  payload <- extract_diagnostic_payload(est)
  spec <- build_plot_spec(payload, type = "diagnostic")

  expect_equal(spec$type, "diagnostic")
  # Should have convergence data if vals exist
  if (!is.null(payload$convergence_vals)) {
    expect_true("convergence" %in% names(spec$data))
    expect_true(all(c("iteration", "rmse") %in% names(spec$data$convergence)))
  }
  # Should have pretreatment fit data
  expect_true("pretreatment_fit" %in% names(spec$data))
  expect_true("residuals" %in% names(spec$data))
  expect_true(!is.null(spec$annotations$rmse_pre))
})

test_that("render_diagnostic_plot returns ggplot for convergence", {
  est <- make_test_estimate()
  payload <- extract_diagnostic_payload(est)
  spec <- build_plot_spec(payload, type = "diagnostic")

  p <- render_diagnostic_plot(spec, subtype = "convergence")
  expect_s3_class(p, "gg")
})

test_that("render_diagnostic_plot returns ggplot for fit", {
  est <- make_test_estimate()
  payload <- extract_diagnostic_payload(est)
  spec <- build_plot_spec(payload, type = "diagnostic")

  p <- render_diagnostic_plot(spec, subtype = "fit")
  expect_s3_class(p, "gg")
})

test_that("plot.synthdid_estimate type='diagnostic' works", {
  est <- make_test_estimate()

  p <- plot(est, type = "diagnostic", subtype = "convergence")
  expect_s3_class(p, "gg")

  p2 <- plot(est, type = "diagnostic", subtype = "fit")
  expect_s3_class(p2, "gg")
})

test_that("synthdid_rmse_plot routes single estimate to new engine", {
  est <- make_test_estimate()
  lifecycle::expect_deprecated(
    p <- synthdid_rmse_plot(est)
  )
  expect_s3_class(p, "gg")
})


# =============================================================================
# Scalability guardrail tests
# =============================================================================

test_that("auto mode emits message when downsampling controls to top_k", {
  est <- make_test_estimate()
  payload <- extract_plot_payload.synthdid_estimate(est, include_controls = TRUE)
  # Override N0 to trigger top_k and verify message
  payload$N0 <- 200

  expect_message(
    spec <- build_plot_spec(payload, type = "trajectory", mode = "auto",
                            top_k = 5),
    "Showing top 5 of"
  )
  expect_equal(spec$mode, "top_k")
})

test_that("max_points guardrail stops on mode='full' with too many points", {
  # Build a fully consistent fake payload to trigger the guardrail
  n_time <- 300
  n_ctrl <- 200
  payload <- list(
    treated_trajectory = rep(1, n_time),
    synthetic_trajectory = rep(1, n_time),
    effect_curve = rep(0, 10),
    omega = rep(1 / n_ctrl, n_ctrl),
    lambda = rep(1 / (n_time - 10), n_time - 10),
    time = seq_len(n_time),
    time_labels = as.character(seq_len(n_time)),
    T0 = n_time - 10,
    T1 = 10,
    N0 = n_ctrl,
    N1 = 1,
    tau = 0,
    se = NA_real_,
    estimator = "test",
    unit_names = paste0("U", seq_len(n_ctrl + 1)),
    control_trajectories = matrix(1, nrow = n_ctrl, ncol = n_time,
                                  dimnames = list(paste0("U", seq_len(n_ctrl)),
                                                  NULL))
  )

  expect_error(
    build_plot_spec(payload, type = "trajectory", mode = "full"),
    "data points"
  )
})

test_that("max_points guardrail can be overridden with force=TRUE", {
  n_time <- 300
  n_ctrl <- 200
  payload <- list(
    treated_trajectory = rep(1, n_time),
    synthetic_trajectory = rep(1, n_time),
    effect_curve = rep(0, 10),
    omega = rep(1 / n_ctrl, n_ctrl),
    lambda = rep(1 / (n_time - 10), n_time - 10),
    time = seq_len(n_time),
    time_labels = as.character(seq_len(n_time)),
    T0 = n_time - 10,
    T1 = 10,
    N0 = n_ctrl,
    N1 = 1,
    tau = 0,
    se = NA_real_,
    estimator = "test",
    unit_names = paste0("U", seq_len(n_ctrl + 1)),
    control_trajectories = matrix(1, nrow = n_ctrl, ncol = n_time,
                                  dimnames = list(paste0("U", seq_len(n_ctrl)),
                                                  NULL))
  )

  spec <- build_plot_spec(payload, type = "trajectory", mode = "full",
                          force = TRUE)
  expect_true("controls" %in% names(spec$data))
})


# =============================================================================
# Weights subtypes tests
# =============================================================================

test_that("render_weights_plot subtype='omega' returns a ggplot", {
  est <- make_test_estimate()
  payload <- extract_plot_payload.synthdid_estimate(est)
  spec <- build_plot_spec(payload, type = "weights")

  p <- render_weights_plot(spec, subtype = "omega")
  expect_s3_class(p, "gg")
})

test_that("render_weights_plot subtype='lambda' returns a ggplot", {
  est <- make_test_estimate()
  payload <- extract_plot_payload.synthdid_estimate(est)
  spec <- build_plot_spec(payload, type = "weights")

  p <- render_weights_plot(spec, subtype = "lambda")
  expect_s3_class(p, "gg")
})

test_that("render_weights_plot subtype='cumulative' returns a ggplot", {
  est <- make_test_estimate()
  payload <- extract_plot_payload.synthdid_estimate(est)
  spec <- build_plot_spec(payload, type = "weights")

  p <- render_weights_plot(spec, subtype = "cumulative")
  expect_s3_class(p, "gg")
})

test_that("build_weights_spec includes cumulative and lambda data", {
  est <- make_test_estimate()
  payload <- extract_plot_payload.synthdid_estimate(est)
  spec <- build_plot_spec(payload, type = "weights")

  expect_true("cumulative" %in% names(spec$data))
  expect_true(all(c("rank", "cumulative") %in% names(spec$data$cumulative)))
  expect_true("lambda_weights" %in% names(spec$data))
  expect_true(all(c("time", "weight") %in% names(spec$data$lambda_weights)))
  expect_true(!is.null(spec$annotations$n_effective))
  expect_true(!is.null(spec$annotations$n_90pct))
  expect_true(!is.null(spec$annotations$n_effective_time))
})

test_that("plot.synthdid_estimate type='weights' subtype='lambda' works", {
  est <- make_test_estimate()
  p <- plot(est, type = "weights", subtype = "lambda")
  expect_s3_class(p, "gg")
})


# =============================================================================
# Gap sub-panel tests
# =============================================================================

test_that("trajectory plot with show_gap=TRUE returns combined plot", {
  est <- make_test_estimate()
  p <- plot(est, type = "trajectory", show_gap = TRUE)
  # patchwork is available, so this should be a patchwork object
  expect_true(inherits(p, "patchwork") || inherits(p, "gg"))
})


test_that("staggered effect plots include uncertainty/context subtitles", {
  stag <- make_test_staggered()
  attr(stag, "se") <- 0.05

  p_cohort <- plot(stag, type = "effect", subtype = "cohort")
  expect_match(p_cohort$labels$subtitle, "Aggregate ATT")
  expect_match(p_cohort$labels$subtitle, "\\[")

  p_event <- plot(stag, type = "effect", subtype = "event")
  expect_match(p_event$labels$subtitle, "point estimates")
})
