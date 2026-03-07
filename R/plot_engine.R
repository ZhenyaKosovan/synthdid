# Plot Engine for synthdid
#
# Layered architecture: extract -> spec -> render
# This file implements the new plot engine that replaces the legacy synthdid_plot()
# for default plot() calls, while keeping the legacy code untouched.
#
# Layer 1: Extractors  — pull data from estimate objects into a standardized payload
# Layer 2: Spec builder — transform payload into tidy data frames ready for ggplot
# Layer 3: Renderers   — produce ggplot objects from specs
#
# All functions in this file are internal (not exported).


# =============================================================================
# LAYER 1: EXTRACTORS
# =============================================================================

#' Extract plot payload from a synthdid_estimate object
#' @param object A synthdid_estimate object
#' @param include_controls If TRUE, include individual control trajectories
#' @return A list with standardized plot data
#' @keywords internal
extract_plot_payload.synthdid_estimate <- function(object, include_controls = FALSE) {
  setup <- attr(object, "setup")
  weights <- attr(object, "weights")
  Y <- setup$Y
  N0 <- setup$N0
  T0 <- setup$T0
  N1 <- nrow(Y) - N0
  T1 <- ncol(Y) - T0

  X.beta <- contract3(setup$X, weights$beta)
  omega <- weights$omega
  lambda <- weights$lambda

  # Treated trajectory: average of treated units
  treated_trajectory <- colMeans(Y[(N0 + 1):nrow(Y), , drop = FALSE])

  # Synthetic control trajectory: omega-weighted control average
  synthetic_trajectory <- as.numeric(t(omega) %*% (Y[1:N0, ] - X.beta[1:N0, ]))
  if (N1 > 0) {
    # Add back the average covariate adjustment for treated units
    synthetic_trajectory <- synthetic_trajectory +
      colMeans(X.beta[(N0 + 1):nrow(Y), , drop = FALSE])
  }

  # Effect curve (same logic as synthdid_effect_curve, inline)
  tau.sc <- as.numeric(t(c(-omega, rep(1 / N1, N1))) %*% (Y - X.beta))
  effect_curve <- tau.sc[T0 + (1:T1)] - as.numeric(tau.sc[1:T0] %*% lambda)

  # Time axis
  time_labels <- colnames(Y)
  if (is.null(time_labels)) {
    time_labels <- seq_len(ncol(Y))
  }
  time_numeric <- suppressWarnings(as.numeric(time_labels))
  if (any(is.na(time_numeric))) {
    time_numeric <- seq_len(ncol(Y))
  }

  # SE and estimator info
  se_val <- attr(object, "se")
  if (is.null(se_val) || (length(se_val) == 1 && is.na(se_val))) {
    se_val <- NA_real_
  }
  estimator <- attr(object, "estimator")
  if (is.null(estimator)) estimator <- "synthdid_estimate"

  payload <- list(
    treated_trajectory = treated_trajectory,
    synthetic_trajectory = synthetic_trajectory,
    effect_curve = effect_curve,
    omega = omega,
    lambda = lambda,
    time = time_numeric,
    time_labels = time_labels,
    T0 = T0,
    T1 = T1,
    N0 = N0,
    N1 = N1,
    tau = as.numeric(object),
    se = se_val,
    estimator = estimator,
    unit_names = rownames(Y)
  )

  if (include_controls) {
    payload$control_trajectories <- Y[1:N0, , drop = FALSE]
  }

  payload
}


#' Extract plot payload from a synthdid_staggered object
#' @param object A synthdid_staggered object
#' @param type One of "trajectory", "weights", "effect"
#' @return A list with standardized plot data
#' @keywords internal
extract_plot_payload.synthdid_staggered <- function(object, type = "effect") {
  if (type == "trajectory") {
    # Per-cohort payloads from subproblem fits
    cohort_payloads <- lapply(
      names(object$subproblem_fits),
      function(g_char) {
        fit <- object$subproblem_fits[[g_char]]
        p <- extract_plot_payload.synthdid_estimate(fit, include_controls = FALSE)
        p$cohort_time <- as.integer(g_char)
        p
      }
    )
    names(cohort_payloads) <- names(object$subproblem_fits)
    return(list(
      type = "trajectory",
      cohort_payloads = cohort_payloads,
      att = object$att,
      method = object$method
    ))
  }

  if (type == "weights") {
    # Collect omega vectors and aggregation weights
    omega_list <- lapply(
      names(object$subproblem_fits),
      function(g_char) {
        fit <- object$subproblem_fits[[g_char]]
        w <- attr(fit, "weights")$omega
        setup <- attr(fit, "setup")
        names(w) <- rownames(setup$Y)[1:setup$N0]
        list(
          cohort_time = as.integer(g_char),
          omega = w,
          aggregation_weight = object$aggregation_weights[g_char]
        )
      }
    )
    names(omega_list) <- names(object$subproblem_fits)
    return(list(
      type = "weights",
      omega_list = omega_list,
      aggregation_weights = object$aggregation_weights,
      method = object$method
    ))
  }

  # type == "effect" (default)
  se_val <- attr(object, "se")
  if (is.null(se_val) || (length(se_val) == 1 && is.na(se_val))) {
    se_val <- NA_real_
  }
  event_data <- stats::predict(object, type = "event")
  list(
    type = "effect",
    cohort_effects = object$cohort_effects,
    event_data = event_data,
    att = object$att,
    se = se_val,
    aggregation_weights = object$aggregation_weights,
    method = object$method
  )
}


#' Extract diagnostic plot payload from a synthdid_estimate object
#' @param object A synthdid_estimate object
#' @return A list with convergence trace and pre-treatment fit info
#' @keywords internal
extract_diagnostic_payload <- function(object) {
  setup <- attr(object, "setup")
  weights <- attr(object, "weights")
  Y <- setup$Y
  N0 <- setup$N0
  T0 <- setup$T0
  N1 <- nrow(Y) - N0
  T1 <- ncol(Y) - T0

  # Convergence trace (RMSE values across iterations)
  vals <- weights$vals
  if (!is.null(vals)) {
    vals <- sqrt(vals[!is.na(vals)])
  }

  # Pre-treatment fit: compare treated vs intercept-adjusted synthetic
  # SDID allows a level gap between treated and synthetic (handled by intercept),

  # so we shift the synthetic by the lambda-weighted pre-treatment gap to assess
  # how well it tracks the *shape* of the treated series.
  X.beta <- contract3(setup$X, weights$beta)
  omega <- weights$omega
  lambda <- weights$lambda
  treated_pre <- colMeans(Y[(N0 + 1):nrow(Y), 1:T0, drop = FALSE])
  synthetic_pre <- as.numeric(t(omega) %*% (Y[1:N0, 1:T0, drop = FALSE] -
                                              X.beta[1:N0, 1:T0, drop = FALSE]))
  if (N1 > 0) {
    synthetic_pre <- synthetic_pre +
      colMeans(X.beta[(N0 + 1):nrow(Y), 1:T0, drop = FALSE])
  }
  # Intercept adjustment: shift synthetic to match treated level
  intercept <- sum(lambda * (treated_pre - synthetic_pre))
  synthetic_pre <- synthetic_pre + intercept
  residuals_pre <- treated_pre - synthetic_pre
  rmse_pre <- sqrt(mean(residuals_pre^2))

  time_labels <- colnames(Y)
  if (is.null(time_labels)) time_labels <- seq_len(ncol(Y))
  time_numeric <- suppressWarnings(as.numeric(time_labels))
  if (any(is.na(time_numeric))) time_numeric <- seq_len(ncol(Y))

  estimator <- attr(object, "estimator")
  if (is.null(estimator)) estimator <- "synthdid_estimate"

  list(
    convergence_vals = vals,
    treated_pre = treated_pre,
    synthetic_pre = synthetic_pre,
    residuals_pre = residuals_pre,
    rmse_pre = rmse_pre,
    time_pre = time_numeric[1:T0],
    time_labels_pre = as.character(time_labels[1:T0]),
    T0 = T0,
    N0 = N0,
    estimator = estimator
  )
}


# =============================================================================
# LAYER 2: SPEC BUILDER
# =============================================================================

#' Build a plot specification from a payload
#' @param payload Output from an extract_plot_payload function
#' @param type Plot type: "trajectory", "weights", or "effect"
#' @param mode "auto", "full", or "top_k"
#' @param top_k Number of top units to show in top_k mode
#' @param ... Additional options
#' @return A list with type, mode, data (list of data frames), annotations, options
#' @keywords internal
build_plot_spec <- function(payload, type, mode = "auto", top_k = NULL, ...) {
  dots <- list(...)

  # Resolve mode
  if (mode == "auto") {
    n0 <- payload$N0 %||% 0
    mode <- if (n0 <= SYNTHDID_PLOT_AUTO_THRESHOLD) "full" else "top_k"
  }
  if (is.null(top_k)) top_k <- SYNTHDID_PLOT_TOP_K_DEFAULT

  spec <- list(
    type = type,
    mode = mode,
    top_k = top_k,
    data = list(),
    annotations = list(),
    options = dots
  )

  if (type == "trajectory") {
    spec <- build_trajectory_spec(spec, payload)
  } else if (type == "weights") {
    spec <- build_weights_spec(spec, payload)
  } else if (type == "effect") {
    spec <- build_effect_spec(spec, payload)
  } else if (type == "diagnostic") {
    spec <- build_diagnostic_spec(spec, payload)
  }

  spec
}


#' Build trajectory spec data frames
#' @keywords internal
build_trajectory_spec <- function(spec, payload) {
  time <- payload$time
  time_labels <- as.character(payload$time_labels)
  time_labels_numeric <- suppressWarnings(as.numeric(time_labels))
  time_is_numeric <- !any(is.na(time_labels_numeric)) &&
    length(time_labels_numeric) == length(time) &&
    all(time_labels_numeric == time)

  # Main trajectories
  traj_df <- data.frame(
    time = rep(time, 2),
    value = c(payload$treated_trajectory, payload$synthetic_trajectory),
    series = rep(c("Treated", "Synthetic Control"), each = length(time)),
    stringsAsFactors = FALSE
  )
  spec$data$trajectories <- traj_df

  # Lambda weights ribbon (bottom of plot)
  lambda_full <- rep(0, length(time))
  lambda_full[seq_len(payload$T0)] <- payload$lambda
  spec$data$lambda <- data.frame(
    time = time,
    weight = lambda_full,
    stringsAsFactors = FALSE
  )

  # Treatment onset
  spec$annotations$onset_time <- time[payload$T0 + 1]

  # Effect annotation
  spec$annotations$tau <- payload$tau
  spec$annotations$se <- payload$se
  spec$annotations$estimator <- payload$estimator
  spec$annotations$time_labels <- time_labels
  spec$annotations$time_is_numeric <- time_is_numeric

  # Control trajectories (if available)
  if (!is.null(payload$control_trajectories)) {
    ctrl <- payload$control_trajectories
    omega <- payload$omega
    original_n <- nrow(ctrl)

    # In top_k mode, select top weighted controls
    if (spec$mode == "top_k" && nrow(ctrl) > spec$top_k) {
      top_idx <- order(omega, decreasing = TRUE)[seq_len(spec$top_k)]
      ctrl <- ctrl[top_idx, , drop = FALSE]
      if (original_n > spec$top_k) {
        message(
          sprintf(
            "Showing top %d of %d control units (mode='%s'). ",
            spec$top_k, original_n, spec$mode
          ),
          "Use mode='full' to show all, or set top_k to change the number shown."
        )
      }
    }

    # Guard against too many points
    n_points <- nrow(ctrl) * ncol(ctrl)
    if (n_points > SYNTHDID_PLOT_MAX_POINTS) {
      if (isTRUE(spec$options$force)) {
        # User explicitly overrode the guardrail
      } else {
        stop(
          sprintf(
            "Control trajectories would produce %s data points (limit: %s). ",
            format(n_points, big.mark = ","),
            format(SYNTHDID_PLOT_MAX_POINTS, big.mark = ",")
          ),
          "Use mode='top_k' to reduce, or pass force=TRUE to override."
        )
      }
    }

    ctrl_names <- rownames(ctrl)
    if (is.null(ctrl_names)) ctrl_names <- paste0("Control_", seq_len(nrow(ctrl)))
    spec$data$controls <- do.call(rbind, lapply(seq_len(nrow(ctrl)), function(i) {
      data.frame(
        time = time,
        value = as.numeric(ctrl[i, ]),
        unit = ctrl_names[i],
        stringsAsFactors = FALSE
      )
    }))
  }

  spec
}


#' Build weights spec data frames
#' @keywords internal
build_weights_spec <- function(spec, payload) {
  omega <- payload$omega
  unit_names <- names(omega)
  if (is.null(unit_names)) {
    if (!is.null(payload$unit_names)) {
      unit_names <- payload$unit_names[seq_len(payload$N0)]
    } else {
      unit_names <- paste0("Unit_", seq_len(length(omega)))
    }
  }

  # Sort descending by weight
  ord <- order(omega, decreasing = TRUE)
  omega_sorted <- omega[ord]
  names_sorted <- unit_names[ord]

  # Top-k with "other" bucket
  label_k <- min(SYNTHDID_PLOT_WEIGHTS_LABEL_K, length(omega_sorted))
  if (length(omega_sorted) > label_k) {
    top_weights <- omega_sorted[seq_len(label_k)]
    top_names <- names_sorted[seq_len(label_k)]
    other_weight <- sum(omega_sorted[(label_k + 1):length(omega_sorted)])

    weights_df <- data.frame(
      unit = c(top_names, "Other"),
      weight = c(top_weights, other_weight),
      stringsAsFactors = FALSE
    )
  } else {
    weights_df <- data.frame(
      unit = names_sorted,
      weight = omega_sorted,
      stringsAsFactors = FALSE
    )
  }

  # Preserve display order (top weight first)
  weights_df$unit <- factor(weights_df$unit,
    levels = rev(weights_df$unit)
  )

  spec$data$weights <- weights_df

  # Cumulative weight curve
  spec$data$cumulative <- data.frame(
    rank = seq_along(omega_sorted),
    cumulative = cumsum(omega_sorted),
    stringsAsFactors = FALSE
  )
  # Effective number of donors (1 / sum(w^2))
  spec$annotations$n_effective <- round(1 / sum(omega_sorted^2))
  spec$annotations$n_total <- length(omega_sorted)
  # Units needed for 90% weight
  spec$annotations$n_90pct <- min(which(cumsum(omega_sorted) >= 0.9),
                                  length(omega_sorted))

  # Lambda (time) weights
  if (!is.null(payload$lambda)) {
    lambda <- payload$lambda
    time_pre <- payload$time[seq_len(payload$T0)]
    time_labels <- as.character(payload$time_labels[seq_len(payload$T0)])
    time_numeric <- suppressWarnings(as.numeric(time_labels))
    spec$data$lambda_weights <- data.frame(
      time = time_pre,
      weight = lambda,
      stringsAsFactors = FALSE
    )
    spec$annotations$lambda_time_labels <- time_labels
    spec$annotations$lambda_time_is_numeric <- !any(is.na(time_numeric)) &&
      length(time_numeric) == length(time_pre) &&
      all(time_numeric == time_pre)
    spec$annotations$n_effective_time <- round(1 / sum(lambda^2))
    spec$annotations$n_pre <- payload$T0
  }

  spec
}


#' Build effect spec data frames
#' @keywords internal
build_effect_spec <- function(spec, payload) {
  # For simultaneous (non-staggered) estimates
  if (!is.null(payload$effect_curve)) {
    post_idx <- payload$T0 + seq_len(payload$T1)
    time_post <- payload$time[post_idx]
    time_labels <- as.character(payload$time_labels)
    time_labels_numeric <- suppressWarnings(as.numeric(time_labels))
    time_is_numeric <- !any(is.na(time_labels_numeric)) &&
      length(time_labels_numeric) == length(payload$time) &&
      all(time_labels_numeric == payload$time)
    effect_df <- data.frame(
      time = time_post,
      effect = payload$effect_curve,
      stringsAsFactors = FALSE
    )
    spec$data$effect <- effect_df
    spec$annotations$tau <- payload$tau
    spec$annotations$se <- payload$se
    spec$annotations$time_labels <- time_labels[post_idx]
    spec$annotations$time_is_numeric <- time_is_numeric
    return(spec)
  }

  # For staggered: cohort effects
  if (!is.null(payload$cohort_effects)) {
    ce <- payload$cohort_effects
    ce$cohort_label <- factor(
      paste0("t=", ce$cohort_time),
      levels = paste0("t=", sort(ce$cohort_time))
    )
    spec$data$cohort_effects <- ce
    spec$annotations$att <- payload$att
    spec$annotations$se <- payload$se
  }

  # For staggered: event study
  if (!is.null(payload$event_data)) {
    ed <- payload$event_data
    ed$cohort_label <- factor(paste0("t=", ed$cohort_time))
    spec$data$event_data <- ed
  }

  spec
}


#' Build diagnostic spec data frames
#' @keywords internal
build_diagnostic_spec <- function(spec, payload) {
  # Convergence trace
  if (!is.null(payload$convergence_vals) && length(payload$convergence_vals) > 0) {
    spec$data$convergence <- data.frame(
      iteration = seq_along(payload$convergence_vals),
      rmse = payload$convergence_vals,
      stringsAsFactors = FALSE
    )
  }

  # Pre-treatment fit
  if (!is.null(payload$treated_pre)) {
    spec$data$pretreatment_fit <- data.frame(
      time = rep(payload$time_pre, 2),
      value = c(payload$treated_pre, payload$synthetic_pre),
      series = rep(c("Treated", "Synthetic Control"), each = length(payload$time_pre)),
      stringsAsFactors = FALSE
    )
    spec$data$residuals <- data.frame(
      time = payload$time_pre,
      residual = payload$residuals_pre,
      stringsAsFactors = FALSE
    )
    spec$annotations$rmse_pre <- payload$rmse_pre
    spec$annotations$time_labels_pre <- payload$time_labels_pre
    time_numeric <- suppressWarnings(as.numeric(payload$time_labels_pre))
    spec$annotations$time_is_numeric <- !any(is.na(time_numeric)) &&
      length(time_numeric) == length(payload$time_pre) &&
      all(time_numeric == payload$time_pre)
  }

  spec$annotations$estimator <- payload$estimator
  spec
}


#' Select a manageable number of axis breaks
#' @param x Numeric vector of sorted x values
#' @param max_breaks Maximum number of breaks
#' @return Numeric vector of break values
#' @keywords internal
select_axis_breaks <- function(x, max_breaks = 10) {
  ux <- sort(unique(x))
  if (length(ux) <= max_breaks) return(ux)
  idx <- round(seq(1, length(ux), length.out = max_breaks))
  ux[unique(idx)]
}


#' Apply categorical-style labels to a numeric time axis when needed
#' @param p A ggplot object
#' @param time_values Numeric x values used in the plot
#' @param time_labels Display labels corresponding to time_values
#' @param time_is_numeric Logical; TRUE if labels are already numeric values
#' @return A ggplot object
#' @keywords internal
apply_time_axis_labels <- function(p, time_values, time_labels, time_is_numeric) {
  if (isTRUE(time_is_numeric)) return(p)

  ux <- sort(unique(time_values))
  labels_chr <- as.character(time_labels)
  if (length(ux) != length(labels_chr)) return(p)

  breaks <- select_axis_breaks(ux)
  break_labels <- labels_chr[match(breaks, ux)]
  p + ggplot2::scale_x_continuous(breaks = breaks, labels = break_labels)
}


# =============================================================================
# LAYER 3: RENDERERS
# =============================================================================

#' Render a trajectory plot
#' @param spec A plot spec from build_plot_spec with type="trajectory"
#' @param show_gap If TRUE, add a gap (treated - synthetic) sub-panel below.
#'   Requires the patchwork package.
#' @return A ggplot2 object
#' @keywords internal
render_trajectory_plot <- function(spec, show_gap = FALSE) {
  traj <- spec$data$trajectories
  onset <- spec$annotations$onset_time

  p <- ggplot2::ggplot()

  # Control trajectories (behind main lines)
  if (!is.null(spec$data$controls)) {
    p <- p + ggplot2::geom_line(
      data = spec$data$controls,
      ggplot2::aes(x = time, y = value, group = unit),
      color = "grey80",
      linewidth = 0.3,
      alpha = 0.4
    )
  }

  # Lambda weights ribbon along bottom (drawn before main lines so it sits behind)
  lambda_df <- spec$data$lambda
  if (any(lambda_df$weight > 0)) {
    # Push ribbon below the data range so it doesn't overlap trajectories
    y_range <- range(traj$value)
    gap <- diff(y_range) * 0.03
    ribbon_height <- diff(y_range) / SYNTHDID_LAMBDA_PLOT_SCALE_DEFAULT
    lambda_df$y_base <- y_range[1] - gap - ribbon_height
    lambda_df$y_top <- lambda_df$y_base +
      lambda_df$weight * ribbon_height / max(lambda_df$weight)

    p <- p + ggplot2::geom_ribbon(
      data = lambda_df[lambda_df$weight > 0, ],
      ggplot2::aes(x = time, ymin = y_base, ymax = y_top),
      fill = "#0072B2",
      alpha = 0.25
    )
  }

  # Main trajectories
  p <- p + ggplot2::geom_line(
    data = traj,
    ggplot2::aes(x = time, y = value, color = series, linetype = series),
    linewidth = SYNTHDID_LINE_WIDTH_DEFAULT * 1.5
  )

  # Treatment onset line
  p <- p + ggplot2::geom_vline(
    xintercept = onset,
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.4
  )

  # Treatment effect annotation
  tau <- spec$annotations$tau
  se <- spec$annotations$se
  tau_fmt <- sprintf("%.1f", tau)
  label_text <- if (!is.na(se)) {
    sprintf("ATT = %s, SE = %s", tau_fmt, sprintf("%.1f", se))
  } else {
    sprintf("ATT = %s", tau_fmt)
  }

  p <- p +
    ggplot2::scale_color_manual(
      values = c("Treated" = "#E69F00", "Synthetic Control" = "#0072B2"),
      name = NULL
    ) +
    ggplot2::scale_linetype_manual(
      values = c("Treated" = "solid", "Synthetic Control" = "22"),
      name = NULL
    ) +
    ggplot2::labs(
      x = "Time",
      y = "Outcome",
      subtitle = label_text
    ) +
    ggplot2::theme_light() +
    ggplot2::theme(
      legend.position = "bottom"
    )

  p <- apply_time_axis_labels(
    p,
    time_values = traj$time,
    time_labels = spec$annotations$time_labels,
    time_is_numeric = spec$annotations$time_is_numeric
  )

  # Gap sub-panel: treated - synthetic over time
  if (isTRUE(show_gap)) {
    treated_vals <- traj$value[traj$series == "Treated"]
    synth_vals <- traj$value[traj$series == "Synthetic Control"]
    gap_time <- traj$time[traj$series == "Treated"]
    gap_df <- data.frame(
      time = gap_time,
      gap = treated_vals - synth_vals,
      post = gap_time >= onset,
      stringsAsFactors = FALSE
    )

    p_gap <- ggplot2::ggplot(gap_df, ggplot2::aes(x = time, y = gap)) +
      ggplot2::geom_hline(yintercept = 0, color = "grey50", linewidth = 0.3) +
      ggplot2::geom_vline(xintercept = onset, linetype = "dashed",
                          color = "grey40", linewidth = 0.4) +
      ggplot2::geom_area(
        data = gap_df[gap_df$post, ],
        fill = "#E69F00", alpha = 0.3
      ) +
      ggplot2::geom_line(color = "#0072B2", linewidth = 0.7) +
      ggplot2::labs(x = "Time", y = "Gap") +
      ggplot2::theme_light() +
      ggplot2::theme(legend.position = "none")

    p_gap <- apply_time_axis_labels(
      p_gap,
      time_values = gap_df$time,
      time_labels = spec$annotations$time_labels,
      time_is_numeric = spec$annotations$time_is_numeric
    )

    if (requireNamespace("patchwork", quietly = TRUE)) {
      return(p / p_gap + patchwork::plot_layout(heights = c(2, 1)))
    }
    message("Install 'patchwork' to see trajectory + gap panel; showing trajectory only.")
  }

  p
}


#' Render a weights plot
#' @param spec A plot spec from build_plot_spec with type="weights"
#' @param subtype "omega" (default), "lambda", "cumulative", or "both".
#'   "both" combines omega lollipop and cumulative curve via patchwork.
#' @return A ggplot2 object
#' @keywords internal
render_weights_plot <- function(spec, subtype = "omega") {
  # Staggered: aggregation weights bar chart
  if (!is.null(spec$data$aggregation_weights)) {
    agg <- spec$data$aggregation_weights
    p <- ggplot2::ggplot(agg, ggplot2::aes(x = cohort_label, y = weight)) +
      ggplot2::geom_col(fill = "#0072B2", alpha = 0.7, width = 0.6) +
      ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%.2f", weight)),
        vjust = -0.5, size = 3.2, color = "grey30"
      ) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
      ggplot2::labs(
        x = "Adoption Cohort",
        y = "Aggregation Weight"
      ) +
      ggplot2::theme_light()
    return(p)
  }

  # Lambda (time) weights
  if (subtype == "lambda") {
    return(.render_lambda_weights(spec))
  }

  # Cumulative weight curve
  if (subtype == "cumulative") {
    return(.render_cumulative_weights(spec))
  }

  # Combined: omega + cumulative side by side
  if (subtype == "both") {
    p_omega <- .render_omega_weights(spec)
    p_cumul <- .render_cumulative_weights(spec)
    if (requireNamespace("patchwork", quietly = TRUE)) {
      return(p_omega | p_cumul)
    }
    message("Install 'patchwork' to see omega + cumulative side by side; showing omega only.")
    return(p_omega)
  }

  # Default: omega lollipop
  .render_omega_weights(spec)
}


#' Render omega (unit) weights as a lollipop chart
#' @keywords internal
.render_omega_weights <- function(spec) {
  weights_df <- spec$data$weights
  is_other <- weights_df$unit == "Other"
  weights_df$color <- ifelse(is_other, "grey50", "#0072B2")

  weights_df$pct <- weights_df$weight * 100

  p <- ggplot2::ggplot(weights_df, ggplot2::aes(x = pct, y = unit)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = pct, y = unit, yend = unit),
      color = weights_df$color,
      linewidth = 0.6
    ) +
    ggplot2::geom_point(
      color = weights_df$color,
      size = 2.5
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.1f%%", pct)),
      hjust = -0.3, size = 2.8, color = "grey30"
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
    ggplot2::labs(
      x = "Weight (%)",
      y = NULL,
      subtitle = sprintf("Top unit weights (effective n = %d of %d)",
                          spec$annotations$n_effective %||% 0,
                          spec$annotations$n_total %||% 0)
    ) +
    ggplot2::theme_light()

  p
}


#' Render cumulative weight curve
#' @keywords internal
.render_cumulative_weights <- function(spec) {
  cumul <- spec$data$cumulative
  n_eff <- spec$annotations$n_effective
  n_total <- spec$annotations$n_total
  n_90 <- spec$annotations$n_90pct

  p <- ggplot2::ggplot(cumul, ggplot2::aes(x = rank, y = cumulative * 100)) +
    ggplot2::geom_area(fill = "#0072B2", alpha = 0.2) +
    ggplot2::geom_line(color = "#0072B2", linewidth = 0.7) +
    ggplot2::geom_hline(yintercept = 90, linetype = "dashed",
                        color = "grey50", linewidth = 0.4)

  if (!is.null(n_90) && n_90 <= nrow(cumul)) {
    p <- p +
      ggplot2::geom_point(
        data = cumul[n_90, ],
        ggplot2::aes(x = rank, y = cumulative * 100),
        color = "#E69F00", size = 3
      ) +
      ggplot2::annotate(
        "text",
        x = n_90 + max(1, nrow(cumul) * 0.03),
        y = cumul$cumulative[n_90] * 100,
        label = sprintf("%d units -> 90%%", n_90),
        color = "grey40", size = 3.2, hjust = 0
      )
  }

  p <- p +
    ggplot2::scale_y_continuous(
      limits = c(0, 100),
      labels = function(x) paste0(x, "%")
    ) +
    ggplot2::labs(
      x = "Ranked Control Units",
      y = "Cumulative Weight (%)",
      subtitle = sprintf("Effective donors = %d of %d",
                          n_eff %||% 0, n_total %||% 0)
    ) +
    ggplot2::theme_light()

  p
}


#' Render lambda (time) weights as a bar chart
#' @keywords internal
.render_lambda_weights <- function(spec) {
  if (is.null(spec$data$lambda_weights)) {
    stop("No lambda (time) weight data available in this spec.")
  }
  lw <- spec$data$lambda_weights
  n_eff_t <- spec$annotations$n_effective_time
  n_pre <- spec$annotations$n_pre
  equal_weight <- if (!is.null(n_pre) && n_pre > 0) 1 / n_pre else NULL

  p <- ggplot2::ggplot(lw, ggplot2::aes(x = time, y = weight)) +
    ggplot2::geom_col(fill = "#E69F00", alpha = 0.7, width = 0.7)

  if (!is.null(equal_weight)) {
    p <- p +
      ggplot2::geom_hline(yintercept = equal_weight, linetype = "dashed",
                          color = "grey50", linewidth = 0.4) +
      ggplot2::annotate(
        "text",
        x = min(lw$time) + diff(range(lw$time)) * 0.02,
        y = equal_weight,
        label = "Equal weighting",
        color = "grey40", size = 3, hjust = 0, vjust = -0.5
      )
  }

  p <- p +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.1))
    ) +
    ggplot2::labs(
      x = "Time (pre-treatment)",
      y = "Weight",
      subtitle = sprintf("Time weights (effective periods = %d of %d)",
                          n_eff_t %||% 0, n_pre %||% 0)
    ) +
    ggplot2::theme_light()

  p <- apply_time_axis_labels(
    p,
    time_values = lw$time,
    time_labels = spec$annotations$lambda_time_labels,
    time_is_numeric = spec$annotations$lambda_time_is_numeric
  )

  p
}


#' Render an effect plot
#' @param spec A plot spec from build_plot_spec with type="effect"
#' @param subtype For staggered objects: "cohort" or "event"
#' @return A ggplot2 object
#' @keywords internal
render_effect_plot <- function(spec, subtype = "cohort") {
  # Simultaneous: effect curve
  if (!is.null(spec$data$effect)) {
    effect_df <- spec$data$effect
    tau <- spec$annotations$tau
    se <- spec$annotations$se

    p <- ggplot2::ggplot(effect_df, ggplot2::aes(x = time, y = effect)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                          color = "grey50", linewidth = 0.4) +
      ggplot2::geom_hline(
        yintercept = tau,
        linetype = "dotted",
        color = "#E69F00",
        linewidth = 0.6
      ) +
      ggplot2::geom_line(color = "#0072B2", linewidth = 0.8) +
      ggplot2::geom_point(color = "#0072B2", size = 2)

    # ATT label + CI error bar at right edge
    x_max <- max(effect_df$time)
    if (!is.na(se)) {
      ci_lo <- tau - SYNTHDID_CI_Z_95 * se
      ci_hi <- tau + SYNTHDID_CI_Z_95 * se
      # Error bar just past the last data point
      x_span <- diff(range(effect_df$time))
      if (!is.finite(x_span) || x_span <= 0) x_span <- 1
      x_nudge <- x_span * 0.04
      p <- p +
        ggplot2::geom_errorbar(
          data = data.frame(x = x_max + x_nudge, ymin = ci_lo, ymax = ci_hi),
          ggplot2::aes(x = x, ymin = ymin, ymax = ymax),
          inherit.aes = FALSE,
          width = x_span * 0.02,
          color = "#E69F00", linewidth = 0.5
        ) +
        ggplot2::geom_point(
          data = data.frame(x = x_max + x_nudge, y = tau),
          ggplot2::aes(x = x, y = y),
          inherit.aes = FALSE,
          color = "#E69F00", size = 2
        )
      att_label <- sprintf("ATT = %.1f [%.1f, %.1f]", tau, ci_lo, ci_hi)
    } else {
      att_label <- sprintf("ATT = %.1f", tau)
    }

    p <- p +
      ggplot2::labs(
        x = "Time",
        y = "Treatment Effect",
        subtitle = att_label
      ) +
      ggplot2::scale_x_continuous(
        expand = ggplot2::expansion(mult = c(0.02, 0.08))
      ) +
      ggplot2::theme_light()

    p <- apply_time_axis_labels(
      p,
      time_values = effect_df$time,
      time_labels = spec$annotations$time_labels,
      time_is_numeric = spec$annotations$time_is_numeric
    )

    return(p)
  }

  # Staggered: cohort bar chart
  if (subtype == "cohort" && !is.null(spec$data$cohort_effects)) {
    ce <- spec$data$cohort_effects
    att <- spec$annotations$att
    se <- spec$annotations$se
    has_se <- !(is.null(se) || length(se) == 0 || is.na(se))
    att_label <- if (has_se) {
      ci_lo <- att - SYNTHDID_CI_Z_95 * se
      ci_hi <- att + SYNTHDID_CI_Z_95 * se
      sprintf("Aggregate ATT = %.3f [%.3f, %.3f]", att, ci_lo, ci_hi)
    } else {
      sprintf("Aggregate ATT = %.3f", att)
    }

    p <- ggplot2::ggplot(ce, ggplot2::aes(x = cohort_label, y = estimate))
    if (has_se) {
      p <- p + ggplot2::annotate(
        "rect",
        xmin = -Inf, xmax = Inf,
        ymin = att - SYNTHDID_CI_Z_95 * se,
        ymax = att + SYNTHDID_CI_Z_95 * se,
        fill = "#E69F00",
        alpha = 0.08
      )
    }

    p <- p +
      ggplot2::geom_col(fill = "#0072B2", alpha = 0.7, width = 0.6) +
      ggplot2::geom_text(
        ggplot2::aes(
          label = sprintf("%.3f", estimate),
          vjust = ifelse(estimate >= 0, -0.5, 1.5)
        ),
        size = 3, color = "grey30"
      ) +
      ggplot2::geom_hline(
        yintercept = att, linetype = "dashed", color = "#E69F00",
        linewidth = 0.6
      ) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.1, 0.12))) +
      ggplot2::labs(
        x = "Adoption Cohort",
        y = "Treatment Effect",
        subtitle = att_label
      ) +
      ggplot2::theme_light()

    return(p)
  }

  # Staggered: event study
  if (subtype == "event" && !is.null(spec$data$event_data)) {
    ed <- spec$data$event_data
    att <- spec$annotations$att
    se <- spec$annotations$se
    has_att <- !(is.null(att) || length(att) == 0 || is.na(att))
    has_se <- !(is.null(se) || length(se) == 0 || is.na(se))

    # Sort cohort labels numerically
    cohort_times <- sort(unique(ed$cohort_time))
    cohort_labels <- paste0("t=", cohort_times)
    ed$cohort_label <- factor(
      paste0("t=", ed$cohort_time),
      levels = cohort_labels
    )

    # Start with Okabe-Ito; smoothly extend for >8 cohorts to avoid dropped data.
    okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                   "#0072B2", "#D55E00", "#CC79A7", "#999999")
    n_cohorts <- length(cohort_times)
    palette_values <- if (n_cohorts <= length(okabe_ito)) {
      okabe_ito[seq_len(n_cohorts)]
    } else {
      grDevices::colorRampPalette(okabe_ito)(n_cohorts)
    }
    palette <- setNames(palette_values, cohort_labels)
    event_subtitle <- if (has_att && has_se) {
      ci_lo <- att - SYNTHDID_CI_Z_95 * se
      ci_hi <- att + SYNTHDID_CI_Z_95 * se
      paste0(
        sprintf("Aggregate ATT = %.3f [%.3f, %.3f]. ", att, ci_lo, ci_hi),
        "Cohort paths show point estimates."
      )
    } else {
      "Cohort paths show point estimates."
    }

    p <- ggplot2::ggplot(
      ed,
      ggplot2::aes(x = relative_time, y = effect, color = cohort_label)
    ) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                          color = "grey50", linewidth = 0.4) +
      ggplot2::geom_line(linewidth = 0.7) +
      ggplot2::geom_point(size = 2) +
      ggplot2::scale_color_manual(values = palette) +
      ggplot2::scale_x_continuous(breaks = function(x) {
        seq(ceiling(x[1]), floor(x[2]))
      }) +
      ggplot2::labs(
        x = "Periods Since Treatment",
        y = "Treatment Effect",
        color = "Cohort",
        subtitle = event_subtitle
      ) +
      ggplot2::theme_light() +
      ggplot2::theme(legend.position = "bottom")

    return(p)
  }

  stop("Cannot render effect plot: no suitable data in spec")
}


#' Render a diagnostic plot
#' @param spec A plot spec from build_plot_spec with type="diagnostic"
#' @param subtype "convergence", "fit", or "both" (default)
#' @return A ggplot2 object
#' @keywords internal
render_diagnostic_plot <- function(spec, subtype = "both") {
  has_convergence <- !is.null(spec$data$convergence)
  has_fit <- !is.null(spec$data$pretreatment_fit)

  if (subtype == "convergence" || (subtype == "both" && has_convergence && !has_fit)) {
    if (!has_convergence) {
      stop("No convergence trace available for this estimate.")
    }
    conv <- spec$data$convergence
    p <- ggplot2::ggplot(conv, ggplot2::aes(x = iteration, y = rmse)) +
      ggplot2::geom_line(color = "#0072B2", linewidth = 0.7) +
      ggplot2::scale_y_log10() +
      ggplot2::labs(
        x = "Iteration",
        y = "RMSE (log scale)",
        subtitle = paste("Convergence trace:", spec$annotations$estimator)
      ) +
      ggplot2::theme_light()
    return(p)
  }

  if (subtype == "fit" || (subtype == "both" && !has_convergence && has_fit)) {
    if (!has_fit) {
      stop("No pre-treatment fit data available.")
    }
    return(.render_pretreatment_fit(spec))
  }

  # subtype == "both" and we have both
  if (!has_convergence || !has_fit) {
    # Render whichever is available
    if (has_convergence) return(render_diagnostic_plot(spec, subtype = "convergence"))
    if (has_fit) return(render_diagnostic_plot(spec, subtype = "fit"))
    stop("No diagnostic data available for this estimate.")
  }

  p1 <- render_diagnostic_plot(spec, subtype = "convergence")
  p2 <- .render_pretreatment_fit(spec)

  if (requireNamespace("patchwork", quietly = TRUE)) {
    return(p1 / p2 + patchwork::plot_layout(heights = c(1, 2)))
  }
  # Without patchwork, return the fit plot (more informative)
  message("Install 'patchwork' to see convergence + fit side by side; showing fit only.")
  p2
}


#' Render the pre-treatment fit panel of a diagnostic plot
#' @param spec A plot spec with pretreatment_fit data
#' @return A ggplot2 object (or patchwork composite if patchwork is available)
#' @keywords internal
.render_pretreatment_fit <- function(spec) {
  fit_df <- spec$data$pretreatment_fit
  resid_df <- spec$data$residuals
  rmse <- spec$annotations$rmse_pre

  # Split trajectories for ribbon
  treated_vals <- fit_df$value[fit_df$series == "Treated"]
  synth_vals <- fit_df$value[fit_df$series == "Synthetic Control"]
  time_vals <- fit_df$time[fit_df$series == "Treated"]
  ribbon_df <- data.frame(
    time = time_vals,
    ymin = pmin(treated_vals, synth_vals),
    ymax = pmax(treated_vals, synth_vals),
    stringsAsFactors = FALSE
  )

  # MAPE
  mape <- mean(abs(resid_df$residual / treated_vals)) * 100
  fit_verdict <- if (mape < 5) "Good fit" else if (mape < 10) "Acceptable" else "Check fit"

  # Main fit panel with ribbon between trajectories
  p_fit <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = ribbon_df,
      ggplot2::aes(x = time, ymin = ymin, ymax = ymax),
      fill = "#56B4E9", alpha = 0.2
    ) +
    ggplot2::geom_line(
      data = fit_df,
      ggplot2::aes(x = time, y = value, color = series, linetype = series),
      linewidth = 0.8
    ) +
    ggplot2::scale_color_manual(
      values = c("Treated" = "#E69F00", "Synthetic Control" = "#0072B2"),
      name = NULL
    ) +
    ggplot2::scale_linetype_manual(
      values = c("Treated" = "solid", "Synthetic Control" = "22"),
      name = NULL
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Outcome",
      subtitle = sprintf("Pre-treatment %s: RMSE = %.3f, MAPE = %.1f%%",
                          fit_verdict, rmse, mape)
    ) +
    ggplot2::theme_light() +
    ggplot2::theme(legend.position = "bottom")

  p_fit <- apply_time_axis_labels(
    p_fit,
    time_values = time_vals,
    time_labels = spec$annotations$time_labels_pre,
    time_is_numeric = spec$annotations$time_is_numeric
  )

  # Residuals panel (percentage of treated)
  resid_df$pct_residual <- 100 * resid_df$residual / treated_vals

  p_resid <- ggplot2::ggplot(resid_df, ggplot2::aes(x = time, y = pct_residual)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey50", linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = c(-5, 5), linetype = "dashed",
                        color = "grey70", linewidth = 0.3) +
    ggplot2::geom_col(fill = "#56B4E9", alpha = 0.5, width = 0.7) +
    ggplot2::labs(x = "Time (pre-treatment)", y = "Residual (%)") +
    ggplot2::theme_light() +
    ggplot2::theme(legend.position = "none")

  p_resid <- apply_time_axis_labels(
    p_resid,
    time_values = resid_df$time,
    time_labels = spec$annotations$time_labels_pre,
    time_is_numeric = spec$annotations$time_is_numeric
  )

  if (requireNamespace("patchwork", quietly = TRUE)) {
    return(p_fit / p_resid + patchwork::plot_layout(heights = c(2, 1)))
  }
  p_fit
}


#' Compose multiple staggered trajectory plots
#' @param plots List of per-cohort ggplot objects
#' @param has_patchwork Logical indicating whether patchwork is available
#' @return A combined plot object (or first plot with warning if unavailable)
#' @keywords internal
compose_staggered_trajectory_plots <- function(
    plots,
    has_patchwork = requireNamespace("patchwork", quietly = TRUE)
) {
  if (length(plots) == 1) return(plots[[1]])
  if (isTRUE(has_patchwork)) {
    return(Reduce(`+`, plots) +
             patchwork::plot_layout(
               ncol = min(length(plots), 2),
               guides = "collect"
             ) &
             ggplot2::theme(legend.position = "bottom"))
  }
  warning(
    "Package 'patchwork' is not installed; returning only the first cohort trajectory. ",
    "Install 'patchwork' to combine all cohorts in one figure."
  )
  plots[[1]]
}


# =============================================================================
# DISPATCHER
# =============================================================================

#' Render a plot for a synthdid object using the new engine
#'
#' Dispatches extraction by class, builds spec, and routes to the appropriate
#' renderer. This is the main entry point for the new plot engine.
#'
#' @param object A synthdid_estimate or synthdid_staggered object
#' @param type Plot type: "trajectory", "weights", or "effect"
#' @param mode "auto" (default), "full", or "top_k"
#' @param top_k Number of top units to show in top_k mode
#' @param subtype For staggered effect plots: "cohort" (default) or "event"
#' @param include_controls Whether to show individual control unit trajectories
#' @param ... Additional options passed through to spec builder
#' @return A ggplot2 object
#' @keywords internal
synthdid_render_plot <- function(object, type = "trajectory",
                                 mode = "auto", top_k = NULL,
                                 subtype = "cohort",
                                 include_controls = FALSE,
                                 show_gap = FALSE, ...) {
  # Diagnostic type uses its own extractor
  if (type == "diagnostic") {
    if (inherits(object, "synthdid_staggered")) {
      stop("Diagnostic plots are not yet supported for staggered estimates.")
    }
    payload <- extract_diagnostic_payload(object)
    spec <- build_plot_spec(payload, type = "diagnostic", mode = mode,
                            top_k = top_k, ...)
    return(render_diagnostic_plot(spec, subtype = subtype))
  }

  # Extract payload based on class
  if (inherits(object, "synthdid_staggered")) {
    payload <- extract_plot_payload.synthdid_staggered(object, type = type)
  } else {
    payload <- extract_plot_payload.synthdid_estimate(
      object,
      include_controls = include_controls
    )
  }

  # For staggered trajectory: render per-cohort faceted plot
  if (inherits(object, "synthdid_staggered") && type == "trajectory") {
    plots <- lapply(payload$cohort_payloads, function(cp) {
      spec <- build_plot_spec(cp, type = "trajectory", mode = mode,
                              top_k = top_k, ...)
      render_trajectory_plot(spec) +
        ggplot2::ggtitle(paste0("Cohort t=", cp$cohort_time))
    })
    return(compose_staggered_trajectory_plots(plots))
  }

  # For staggered weights: build aggregation weight data
  if (inherits(object, "synthdid_staggered") && type == "weights") {
    agg_weights <- payload$aggregation_weights
    agg_df <- data.frame(
      cohort_label = factor(
        paste0("t=", names(agg_weights)),
        levels = paste0("t=", sort(as.integer(names(agg_weights))))
      ),
      weight = as.numeric(agg_weights),
      stringsAsFactors = FALSE
    )
    spec <- list(
      type = "weights",
      mode = mode,
      data = list(aggregation_weights = agg_df),
      annotations = list(),
      options = list(...)
    )
    return(render_weights_plot(spec))
  }

  # Build spec and render
  spec <- build_plot_spec(payload, type = type, mode = mode,
                          top_k = top_k, ...)

  if (type == "trajectory") {
    render_trajectory_plot(spec, show_gap = show_gap)
  } else if (type == "weights") {
    render_weights_plot(spec, subtype = subtype)
  } else if (type == "effect") {
    render_effect_plot(spec, subtype = subtype)
  } else {
    stop("Unknown plot type: ", type)
  }
}
