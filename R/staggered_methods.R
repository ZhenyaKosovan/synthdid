# S3 Methods for synthdid_staggered Objects

#' Print method for synthdid_staggered
#' @param x A synthdid_staggered object
#' @param ... Additional arguments (currently ignored)
#' @return Invisibly returns x
#' @export
print.synthdid_staggered <- function(x, ...) {
  cat("Staggered Synthetic Difference-in-Differences Estimate\n\n")

  cl <- attr(x, "call")
  if (!is.null(cl)) {
    cat("Call:\n")
    print(cl)
    cat("\n")
  }

  cat("Aggregate ATT: ", format(x$att, digits = 4), "\n")

  se <- attr(x, "se")
  if (!is.null(se) && !is.na(se)) {
    cat("Standard Error: ", format(se, digits = 4),
        " (", attr(x, "se_method"), ")\n", sep = "")
  }

  cat("\nCohorts: ", nrow(x$cohort_effects), "\n")
  cat("Method:  ", x$method, "\n")

  cat("\nCohort-Level Effects:\n")
  ce <- x$cohort_effects
  ce$estimate <- format(round(ce$estimate, 4), nsmall = 4)
  ce$weight <- format(round(ce$weight, 3), nsmall = 3)
  print(ce, row.names = FALSE)

  invisible(x)
}


#' Extract coefficient from synthdid_staggered
#' @param object A synthdid_staggered object
#' @param type Type of coefficient: \code{"aggregate"} returns the scalar ATT,
#'   \code{"cohort"} returns per-cohort effects.
#' @param ... Additional arguments (currently ignored)
#' @return Named numeric vector
#' @export
coef.synthdid_staggered <- function(object, type = c("aggregate", "cohort"), ...) {
  type <- match.arg(type)
  if (type == "aggregate") {
    setNames(object$att, "att")
  } else {
    ce <- object$cohort_effects
    setNames(ce$estimate, paste0("cohort_", ce$cohort_time))
  }
}


#' Confidence intervals for synthdid_staggered
#' @param object A synthdid_staggered object
#' @param parm Ignored (included for S3 generic compatibility)
#' @param level Confidence level (default: 0.95)
#' @param ... Additional arguments (currently ignored)
#' @return A matrix with lower and upper confidence bounds
#' @export
confint.synthdid_staggered <- function(object, parm = NULL, level = 0.95, ...) {
  tau <- object$att
  se <- attr(object, "se")

  if (is.null(se) || is.na(se)) {
    warning("Standard error not available; cannot compute confidence interval. ",
            "Re-estimate with se = TRUE.")
    return(matrix(c(NA, NA), nrow = 1, ncol = 2,
                  dimnames = list("att", c("Lower", "Upper"))))
  }

  z <- qnorm((1 + level) / 2)
  ci <- cbind(
    Lower = tau - z * se,
    Upper = tau + z * se
  )
  rownames(ci) <- "att"
  ci
}


#' Summary for synthdid_staggered
#' @param object A synthdid_staggered object
#' @param ... Additional arguments (currently ignored)
#' @return A summary object of class \code{"summary.synthdid_staggered"}
#' @export
summary.synthdid_staggered <- function(object, ...) {
  se <- attr(object, "se")

  # Compute SE if not cached
  if (is.null(se)) {
    se <- tryCatch(
      sqrt(c(vcov(object, method = "bootstrap"))),
      error = function(e) NA_real_
    )
  }

  summary_obj <- list(
    att = object$att,
    se = se,
    cohort_effects = object$cohort_effects,
    method = object$method,
    n_cohorts = nrow(object$cohort_effects),
    n_units = object$setup$n_units,
    n_periods = object$setup$n_periods,
    control_type = object$setup$control_type,
    call = attr(object, "call")
  )

  class(summary_obj) <- "summary.synthdid_staggered"
  summary_obj
}


#' Print summary of synthdid_staggered
#' @param x A summary.synthdid_staggered object
#' @param ... Additional arguments
#' @return Invisibly returns x
#' @importFrom stats printCoefmat
#' @export
print.summary.synthdid_staggered <- function(x, ...) {
  if (!is.null(x$call)) {
    cat("Call:\n")
    print(x$call)
    cat("\n")
  }

  cat("Staggered Synthetic Difference-in-Differences\n")
  cat("Method: ", x$method, "\n")
  cat("Control set: ", x$control_type, "\n\n")

  # Coefficient table
  tau <- x$att
  se <- x$se
  if (!is.na(se) && se > 0) {
    t_stat <- tau / se
    p_value <- 2 * pnorm(-abs(t_stat))
    coef_table <- cbind(
      Estimate = tau,
      `Std. Error` = se,
      `t value` = t_stat,
      `Pr(>|t|)` = p_value
    )
    rownames(coef_table) <- "ATT"
    cat("Aggregate Treatment Effect:\n")
    printCoefmat(coef_table, digits = 4, signif.stars = TRUE)
  } else {
    cat("Aggregate ATT: ", format(tau, digits = 4), "\n")
    if (is.na(se)) cat("(SE not computed)\n")
  }

  cat("\n")
  cat("Panel:   ", x$n_units, " units, ", x$n_periods, " periods\n", sep = "")
  cat("Cohorts: ", x$n_cohorts, "\n\n")

  cat("Cohort-Level Effects:\n")
  print(x$cohort_effects, row.names = FALSE)
  cat("\n")

  invisible(x)
}


#' Predictions from synthdid_staggered
#' @param object A synthdid_staggered object
#' @param type Type of prediction: \code{"aggregate"} (the scalar ATT),
#'   \code{"cohort"} (per-cohort effects), or \code{"event"} (event-study
#'   style estimates by relative time to treatment).
#' @param ... Additional arguments (currently ignored)
#' @return Depends on \code{type}: a scalar, data frame, or data frame with
#'   event-study columns.
#' @export
predict.synthdid_staggered <- function(object, ...,
                                       type = c("aggregate", "cohort",
                                                 "event")) {
  type <- match.arg(type)

  if (type == "aggregate") {
    return(object$att)
  }

  if (type == "cohort") {
    return(object$cohort_effects)
  }

  # Event study: for each cohort, extract the per-period effect curve
  event_results <- list()
  for (g_char in names(object$subproblem_fits)) {
    fit <- object$subproblem_fits[[g_char]]
    effect_curve <- withCallingHandlers(
      synthdid_effect_curve(fit),
      lifecycle_warning_deprecated = function(cnd) invokeRestart("muffleWarning")
    )

    setup <- attr(fit, "setup")
    T0 <- setup$T0
    T1 <- ncol(setup$Y) - T0

    event_results[[g_char]] <- data.frame(
      cohort_time = as.integer(g_char),
      relative_time = seq_len(T1),
      effect = as.numeric(effect_curve),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, event_results)
}


#' Plot synthdid_staggered
#'
#' Produces plots for staggered synthetic difference-in-differences estimates.
#' The default is an effect plot showing cohort-level treatment effects.
#'
#' For backward compatibility, \code{type = "cohort"} and \code{type = "event"}
#' are translated to \code{type = "effect"} with the corresponding
#' \code{subtype}.
#'
#' @param x A synthdid_staggered object
#' @param type Plot type: \code{"effect"} (default), \code{"trajectory"},
#'   \code{"weights"}, or legacy values \code{"cohort"} / \code{"event"}.
#' @param subtype For effect plots: \code{"cohort"} (default) for cohort
#'   effects bar plot, \code{"event"} for event-study style plot.
#' @param mode Display mode: \code{"auto"} (default), \code{"full"},
#'   or \code{"top_k"}.
#' @param ... Additional arguments passed to the plot engine.
#' @return A ggplot2 object
#' @export
plot.synthdid_staggered <- function(x,
                                    type = c("effect", "trajectory", "weights",
                                             "cohort", "event"),
                                    subtype = c("cohort", "event"),
                                    mode = "auto", ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Plotting requires the 'ggplot2' package.")
  }
  type <- match.arg(type)
  subtype <- match.arg(subtype)

  # Backward compatibility: translate legacy type values
  if (type %in% c("cohort", "event")) {
    lifecycle::deprecate_soft(
      "2.0.0",
      'plot.synthdid_staggered(type = "cohort/event")',
      details = paste0(
        'Use type = "effect", subtype = "', type, '" instead of type = "', type, '".'
      )
    )
    subtype <- type
    type <- "effect"
  }

  synthdid_render_plot(x, type = type, mode = mode, subtype = subtype, ...)
}

