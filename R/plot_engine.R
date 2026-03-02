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
  event_data <- stats::predict(object, type = "event")
  list(
    type = "effect",
    cohort_effects = object$cohort_effects,
    event_data = event_data,
    att = object$att,
    aggregation_weights = object$aggregation_weights,
    method = object$method
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
  }

  spec
}


#' Build trajectory spec data frames
#' @keywords internal
build_trajectory_spec <- function(spec, payload) {
  time <- payload$time

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

  # Control trajectories (if available)
  if (!is.null(payload$control_trajectories)) {
    ctrl <- payload$control_trajectories
    omega <- payload$omega

    # In top_k mode, select top weighted controls
    if (spec$mode == "top_k" && nrow(ctrl) > spec$top_k) {
      top_idx <- order(omega, decreasing = TRUE)[seq_len(spec$top_k)]
      ctrl <- ctrl[top_idx, , drop = FALSE]
    }

    # Guard against too many points
    n_points <- nrow(ctrl) * ncol(ctrl)
    if (n_points <= SYNTHDID_PLOT_MAX_POINTS) {
      ctrl_names <- rownames(ctrl)
      if (is.null(ctrl_names)) ctrl_names <- paste0("Control_", seq_len(nrow(ctrl)))
      ctrl_df <- data.frame(
        time = rep(time, each = nrow(ctrl)),
        value = as.numeric(t(ctrl)),  # transpose so time varies slower
        unit = rep(ctrl_names, times = length(time)),
        stringsAsFactors = FALSE
      )
      # Fix: time varies faster in column-major, so we reshape correctly
      ctrl_long <- do.call(rbind, lapply(seq_len(nrow(ctrl)), function(i) {
        data.frame(
          time = time,
          value = as.numeric(ctrl[i, ]),
          unit = ctrl_names[i],
          stringsAsFactors = FALSE
        )
      }))
      spec$data$controls <- ctrl_long
    }
  }

  spec
}


#' Build weights spec data frames
#' @keywords internal
build_weights_spec <- function(spec, payload) {
  omega <- payload$omega
  unit_names <- names(omega)
  if (is.null(unit_names)) {
    # Fall back to unit names from payload
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
  spec
}


#' Build effect spec data frames
#' @keywords internal
build_effect_spec <- function(spec, payload) {
  # For simultaneous (non-staggered) estimates
  if (!is.null(payload$effect_curve)) {
    time_post <- payload$time[payload$T0 + seq_len(payload$T1)]
    effect_df <- data.frame(
      time = time_post,
      effect = payload$effect_curve,
      stringsAsFactors = FALSE
    )
    spec$data$effect <- effect_df
    spec$annotations$tau <- payload$tau
    spec$annotations$se <- payload$se
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
  }

  # For staggered: event study
  if (!is.null(payload$event_data)) {
    ed <- payload$event_data
    ed$cohort_label <- factor(paste0("t=", ed$cohort_time))
    spec$data$event_data <- ed
  }

  spec
}


# =============================================================================
# LAYER 3: RENDERERS
# =============================================================================

#' Render a trajectory plot
#' @param spec A plot spec from build_plot_spec with type="trajectory"
#' @return A ggplot2 object
#' @keywords internal
render_trajectory_plot <- function(spec) {
  traj <- spec$data$trajectories
  onset <- spec$annotations$onset_time

  p <- ggplot2::ggplot()

  # Control trajectories (behind main lines)
  if (!is.null(spec$data$controls)) {
    p <- p + ggplot2::geom_line(
      data = spec$data$controls,
      ggplot2::aes(x = time, y = value, group = unit),
      color = "grey70",
      linewidth = SYNTHDID_SPAGHETTI_LINE_WIDTH_DEFAULT,
      alpha = SYNTHDID_SPAGHETTI_LINE_ALPHA_DEFAULT
    )
  }

  # Main trajectories
  p <- p + ggplot2::geom_line(
    data = traj,
    ggplot2::aes(x = time, y = value, color = series),
    linewidth = SYNTHDID_LINE_WIDTH_DEFAULT * 1.5
  )

  # Treatment onset line
  p <- p + ggplot2::geom_vline(
    xintercept = onset,
    linetype = "dashed",
    alpha = SYNTHDID_ONSET_ALPHA_DEFAULT
  )

  # Lambda weights ribbon along bottom
  lambda_df <- spec$data$lambda
  if (any(lambda_df$weight > 0)) {
    # Scale lambda to a fraction of the y-range
    y_range <- range(traj$value)
    ribbon_height <- diff(y_range) / SYNTHDID_LAMBDA_PLOT_SCALE_DEFAULT
    lambda_df$y_base <- y_range[1]
    lambda_df$y_top <- y_range[1] + lambda_df$weight * ribbon_height / max(lambda_df$weight)

    p <- p + ggplot2::geom_ribbon(
      data = lambda_df[lambda_df$weight > 0, ],
      ggplot2::aes(x = time, ymin = y_base, ymax = y_top),
      fill = "steelblue",
      alpha = 0.3
    )
  }

  # Treatment effect annotation
  tau <- spec$annotations$tau
  se <- spec$annotations$se
  label_text <- paste0("Effect = ", format(tau, digits = 3))
  if (!is.na(se)) {
    label_text <- paste0(label_text, " (", format(se, digits = 2), ")")
  }

  p <- p +
    ggplot2::scale_color_manual(
      values = c("Treated" = "#E69F00", "Synthetic Control" = "#0072B2"),
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

  p
}


#' Render a weights plot
#' @param spec A plot spec from build_plot_spec with type="weights"
#' @return A ggplot2 object
#' @keywords internal
render_weights_plot <- function(spec) {
  # Staggered: aggregation weights bar chart
  if (!is.null(spec$data$aggregation_weights)) {
    agg <- spec$data$aggregation_weights
    p <- ggplot2::ggplot(agg, ggplot2::aes(x = cohort_label, y = weight)) +
      ggplot2::geom_col(fill = "steelblue", alpha = 0.8) +
      ggplot2::labs(
        x = "Cohort",
        y = "Aggregation Weight",
        title = "Cohort Aggregation Weights"
      ) +
      ggplot2::theme_light()
    return(p)
  }

  # Simultaneous: lollipop chart of omega weights
  weights_df <- spec$data$weights

  p <- ggplot2::ggplot(weights_df, ggplot2::aes(x = weight, y = unit)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = weight, y = unit, yend = unit),
      color = "steelblue",
      linewidth = 0.7
    ) +
    ggplot2::geom_point(
      color = "steelblue",
      size = 2.5
    ) +
    ggplot2::labs(
      x = "Weight",
      y = NULL,
      title = "Unit Weights (omega)"
    ) +
    ggplot2::theme_light()

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
      ggplot2::geom_line(color = "#0072B2", linewidth = 0.8) +
      ggplot2::geom_point(color = "#0072B2", size = 2) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")

    # ATT reference line
    p <- p + ggplot2::geom_hline(
      yintercept = tau,
      linetype = "dotted",
      color = "#E69F00",
      linewidth = 0.6
    )

    # CI band if SE available
    if (!is.na(se)) {
      p <- p + ggplot2::annotate(
        "rect",
        xmin = -Inf, xmax = Inf,
        ymin = tau - SYNTHDID_CI_Z_95 * se,
        ymax = tau + SYNTHDID_CI_Z_95 * se,
        fill = "#E69F00",
        alpha = 0.1
      )
    }

    p <- p +
      ggplot2::labs(
        x = "Time",
        y = "Treatment Effect",
        title = "Effect Curve"
      ) +
      ggplot2::theme_light()

    return(p)
  }

  # Staggered: cohort bar chart
  if (subtype == "cohort" && !is.null(spec$data$cohort_effects)) {
    ce <- spec$data$cohort_effects
    att <- spec$annotations$att

    p <- ggplot2::ggplot(ce, ggplot2::aes(x = cohort_label, y = estimate)) +
      ggplot2::geom_col(ggplot2::aes(alpha = weight), fill = "steelblue") +
      ggplot2::geom_hline(
        yintercept = att, linetype = 2, color = "red",
        linewidth = 0.7
      ) +
      ggplot2::annotate(
        "text",
        x = Inf, y = att,
        label = paste("ATT =", format(att, digits = 3)),
        hjust = 1.1, vjust = -0.5, color = "red", size = 3.5
      ) +
      ggplot2::scale_alpha_continuous(range = c(0.4, 1), guide = "none") +
      ggplot2::labs(
        x = "Adoption Cohort",
        y = "Treatment Effect",
        title = "Staggered SDID: Cohort-Level Effects"
      ) +
      ggplot2::theme_light()

    return(p)
  }

  # Staggered: event study
  if (subtype == "event" && !is.null(spec$data$event_data)) {
    ed <- spec$data$event_data

    p <- ggplot2::ggplot(
      ed,
      ggplot2::aes(x = relative_time, y = effect, color = cohort_label)
    ) +
      ggplot2::geom_line(linewidth = 0.7) +
      ggplot2::geom_point(size = 2) +
      ggplot2::geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
      ggplot2::labs(
        x = "Periods Since Treatment",
        y = "Treatment Effect",
        color = "Cohort",
        title = "Staggered SDID: Event Study"
      ) +
      ggplot2::theme_light()

    return(p)
  }

  stop("Cannot render effect plot: no suitable data in spec")
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
                                 include_controls = FALSE, ...) {
  # Extract payload based on class
  if (inherits(object, "synthdid_staggered")) {
    payload <- extract_plot_payload.synthdid_staggered(object, type = type)
  } else {
    payload <- extract_plot_payload.synthdid_estimate(
      object,
      include_controls = include_controls || type == "trajectory"
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
    # Return the first plot if only one cohort, otherwise use patchwork if available
    if (length(plots) == 1) return(plots[[1]])
    if (requireNamespace("patchwork", quietly = TRUE)) {
      return(Reduce(`+`, plots) +
               patchwork::plot_layout(ncol = min(length(plots), 2)))
    }
    return(plots[[1]])
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
    render_trajectory_plot(spec)
  } else if (type == "weights") {
    render_weights_plot(spec)
  } else if (type == "effect") {
    render_effect_plot(spec, subtype = subtype)
  } else {
    stop("Unknown plot type: ", type)
  }
}
