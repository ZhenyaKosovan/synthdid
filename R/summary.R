#' Outputs a table of important synthetic controls and their corresponding weights, sorted by weight.
#' The table is truncated to exclude synthetic controls that do not matter for any estimate ---
#' for each estimate, the truncated controls may have total weight no larger that 1-mass.
#' @param estimates, a list of estimates output by synthdid_estimate. Or a single estimate.
#' @param sort.by, the index of the estimate to sort by. Defaults to 1.
#' @param mass, which controls the length of the table. Defaults to 0.9.
#' @param weight.type, 'omega' for units, 'lambda' for time periods
#' @return A matrix of weights for the top controls/periods, sorted by importance.
#' @examples
#' \donttest{
#' data(california_prop99)
#' setup <- panel.matrices(california_prop99)
#' tau.hat <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
#'
#' # Show top control units
#' synthdid_controls(tau.hat, weight.type = "omega")
#'
#' # Show top time periods
#' synthdid_controls(tau.hat, weight.type = "lambda")
#'
#' # Compare multiple estimates
#' tau.sc <- sc_estimate(setup$Y, setup$N0, setup$T0)
#' synthdid_controls(list(sdid = tau.hat, sc = tau.sc))
#' }
#' @export synthdid_controls
synthdid_controls <- function(estimates, sort.by = 1,
                              mass = SYNTHDID_CONTROLS_MASS_DEFAULT,
                              weight.type = "omega") {
  if (inherits(estimates, "synthdid_estimate")) {
    estimates <- list(estimates)
  }
  if (is.null(names(estimates))) {
    names(estimates) <- sprintf("estimate %d", 1:length(estimates))
  }
  if (!weight.type %in% c("omega", "lambda")) {
    stop('weight.type must be "omega" or "lambda"')
  }
  weights <- do.call(cbind, lapply(estimates, function(est) {
    attr(est, "weights")[[weight.type]]
  }))
  if (is.null(dim(weights))) {
    dim(weights) <- c(length(weights), 1)
  }

  Y <- attr(estimates[[1]], "setup")$Y
  o <- rev(order(weights[, sort.by]))
  tab <- weights[o, , drop = FALSE]
  rownames(tab) <- if (weight.type == "omega") {
    rownames(Y)[o]
  } else {
    colnames(Y)[o]
  }
  colnames(tab) <- "Weight"
  # truncate table to retain a weight sum of at least mass for each unit
  tab.len <- max(apply(tab, 2, function(col) {
    Position(function(x) {
      x >= mass
    }, cumsum(col), nomatch = nrow(tab))
  }))
  tab[1:tab.len, , drop = FALSE]
}

#' Summarize a synthdid object
#' @param object The object to summarize
#' @param weight.digits The number of digits to use when displaying weights (omega, lambda)
#' @param fast Be fast but less accurate, e.g. jackknife instead of bootstrap se estimate
#' @param ... Additional arguments (currently ignored).
#' @return A summary object with estimate, standard error, weights, dimensions, and convergence info.
#' @examples
#' \donttest{
#' data(california_prop99)
#' setup <- panel.matrices(california_prop99)
#' tau.hat <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
#'
#' # Detailed summary with standard error
#' summary(tau.hat)
#'
#' # Fast summary (uses jackknife instead of bootstrap)
#' summary(tau.hat, fast = TRUE)
#' }
#' @method summary synthdid_estimate
#' @export
summary.synthdid_estimate <- function(object,
                                      weight.digits = SYNTHDID_WEIGHT_DIGITS_DEFAULT,
                                      fast = FALSE, ...) {
  N0 <- attr(object, "setup")$N0
  T0 <- attr(object, "setup")$T0
  desired_method <- if (fast) {
    "jackknife"
  } else {
    "bootstrap"
  }
  se_attr <- attr(object, "se")
  se_method <- attr(object, "se_method")
  se_val <- if (!is.null(se_attr) && !is.null(se_method)) {
    se_attr
  } else {
    sqrt(vcov(object, method = desired_method))
  }

  summary_obj <- list(
    estimate = c(object),
    se = se_val,
    controls = round(synthdid_controls(object, weight.type = "omega"), digits = weight.digits),
    periods = round(synthdid_controls(object, weight.type = "lambda"), digits = weight.digits),
    dimensions = c(
      N1 = nrow(Y(object)) - N0, N0 = N0, N0.effective = round(1 / sum(omega(object)^2), weight.digits),
      T1 = ncol(Y(object)) - T0, T0 = T0, T0.effective = round(1 / sum(lambda(object)^2), weight.digits)
    ),
    convergence = attr(object, "convergence")
  )

  # Add formula interface info if available
  if (inherits(object, "synthdid")) {
    summary_obj$call <- attr(object, "call")
    summary_obj$formula <- attr(object, "formula")
  }

  class(summary_obj) <- c("summary.synthdid_estimate", "summary.synthdid")
  summary_obj
}


#' Print summary of synthdid object
#' @param x A summary.synthdid_estimate object
#' @param ... Additional arguments
#' @return Invisibly returns the summary object
#' @examples
#' \donttest{
#' data(california_prop99)
#' setup <- panel.matrices(california_prop99)
#' tau.hat <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
#' s <- summary(tau.hat)
#' print(s)
#' }
#' @importFrom stats printCoefmat
#' @export
print.summary.synthdid_estimate <- function(x, ...) {
  # Print call if available
  if (!is.null(x$call)) {
    cat("Call:\n")
    print(x$call)
    cat("\n")
  }

  # Print coefficient table
  tau <- x$estimate
  se <- x$se
  t_stat <- tau / se
  p_value <- 2 * pnorm(-abs(t_stat))

  cat("Treatment Effect Estimate:\n")
  coef_table <- cbind(
    Estimate = tau,
    `Std. Error` = se,
    `t value` = t_stat,
    `Pr(>|t|)` = p_value
  )
  rownames(coef_table) <- "treated"
  printCoefmat(coef_table, digits = 4, signif.stars = TRUE)
  cat("\n")

  # Print dimensions
  cat("Dimensions:\n")
  dim_df <- data.frame(
    " " = c(
      "Treated units:", "Control units:", "Effective controls:",
      "Post-treatment periods:", "Pre-treatment periods:", "Effective periods:"
    ),
    Value = c(
      x$dimensions["N1"], x$dimensions["N0"], x$dimensions["N0.effective"],
      x$dimensions["T1"], x$dimensions["T0"], x$dimensions["T0.effective"]
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  print(dim_df, row.names = FALSE, right = FALSE)
  cat("\n")

  # Print top weights
  cat("Top Control Units (omega weights):\n")
  n_show <- min(SYNTHDID_SUMMARY_TOP_N_DEFAULT, nrow(x$controls))
  print(head(x$controls, n_show))
  cat("\n")

  cat("Top Time Periods (lambda weights):\n")
  n_show_periods <- min(SYNTHDID_SUMMARY_TOP_N_DEFAULT, nrow(x$periods))
  print(head(x$periods, n_show_periods))
  cat("\n")

  # Print convergence status
  conv <- x$convergence
  if (!is.null(conv)) {
    cat("Convergence Status:\n")
    if (conv$overall_converged) {
      cat("  Overall: CONVERGED\n")
    } else {
      cat("  Overall: NOT CONVERGED\n")
    }

    # Show component status
    if (!is.null(conv$lambda)) {
      status_icon <- if (conv$lambda$converged) "\u2713" else "\u2717"  # ✓ or ✗
      cat(sprintf("  Lambda:  %s (%d/%d iterations, %.1f%% utilization)\n",
                  status_icon,
                  conv$lambda$iterations,
                  conv$lambda$max_iter,
                  100 * conv$lambda$iterations / conv$lambda$max_iter))
    }

    if (!is.null(conv$omega)) {
      status_icon <- if (conv$omega$converged) "\u2713" else "\u2717"
      cat(sprintf("  Omega:   %s (%d/%d iterations, %.1f%% utilization)\n",
                  status_icon,
                  conv$omega$iterations,
                  conv$omega$max_iter,
                  100 * conv$omega$iterations / conv$omega$max_iter))
    }

    if (!is.null(conv$joint)) {
      status_icon <- if (conv$joint$converged) "\u2713" else "\u2717"
      cat(sprintf("  Joint:   %s (%d/%d iterations, %.1f%% utilization)\n",
                  status_icon,
                  conv$joint$iterations,
                  conv$joint$max_iter,
                  100 * conv$joint$iterations / conv$joint$max_iter))
    }

    # Add recommendation if not converged
    if (!conv$overall_converged) {
      cat("\n  Recommendation: Consider increasing max.iter or relaxing min.decrease.\n")
      cat("  Use synthdid_convergence_info() for detailed diagnostics.\n")
    }
  } else {
    cat("Convergence Status: Not available (estimate from older version)\n")
  }

  invisible(x)
}
