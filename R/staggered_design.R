# Staggered Adoption Design Layer
#
# Internal helpers for decomposing staggered-adoption panels into
# cohort-specific 2x2 SDID subproblems, following Porreca (2022).

#' Compute first treatment time for each unit
#'
#' Given a binary treatment matrix W (units x time), returns the first period
#' in which each unit is treated, or \code{Inf} for never-treated units.
#'
#' @param W A binary treatment matrix with units as rows and time periods as
#'   columns. Row and column names are preserved.
#' @return A named numeric vector of length \code{nrow(W)}, where each element
#'   is the column index of the first treatment period (\code{Inf} for
#'   never-treated).
#' @keywords internal
compute_first_treat_time <- function(W) {
  stopifnot(is.matrix(W), all(W %in% c(0, 1)))
  # For each row, find first column where W == 1
  first_treat <- apply(W, 1, function(row) {
    idx <- which(row == 1)
    if (length(idx) == 0) Inf else min(idx)
  })
  names(first_treat) <- rownames(W)
  first_treat
}


#' Validate a staggered panel
#'
#' Checks that a panel's treatment matrix satisfies the assumptions required
#' for staggered SDID estimation:
#' \itemize{
#'   \item Balanced panel (rectangular Y and W)
#'   \item Binary treatment indicator
#'   \item No treated units in the first period
#'   \item Absorbing treatment (once treated, stays treated)
#' }
#'
#' @param Y Outcome matrix (units x time).
#' @param W Binary treatment matrix (units x time), same dimensions as Y.
#' @param absorbing Logical. If \code{TRUE} (default), verify that treatment is
#'   absorbing: once a unit is treated, it remains treated for all subsequent
#'   periods.
#' @return Invisibly returns \code{TRUE} if all checks pass. Throws an
#'   informative error otherwise.
#' @keywords internal
validate_staggered_panel <- function(Y, W, absorbing = TRUE) {
  # Dimension checks
  if (!is.matrix(Y) || !is.matrix(W)) {
    stop("Y and W must be matrices.")
  }
  if (!identical(dim(Y), dim(W))) {
    stop("Y and W must have the same dimensions.")
  }
  if (any(is.na(Y))) {
    stop("Missing values in outcome matrix Y.")
  }
  if (any(is.na(W))) {
    stop("Missing values in treatment matrix W.")
  }

  # Binary treatment
  if (!all(W %in% c(0, 1))) {
    stop("Treatment matrix W must be binary (0/1).")
  }

  # No treatment in first period
  if (any(W[, 1] == 1)) {
    stop("No unit may be treated in the first period. ",
         "At least one pre-treatment period is required for all units.")
  }

  # Must have at least some treated units
  first_treat <- compute_first_treat_time(W)
  if (all(is.infinite(first_treat))) {
    stop("No treated units found in the panel.")
  }

  # Absorbing treatment check

  if (absorbing) {
    for (i in which(is.finite(first_treat))) {
      t_start <- first_treat[i]
      if (t_start <= ncol(W) && !all(W[i, t_start:ncol(W)] == 1)) {
        unit_name <- if (!is.null(rownames(W))) rownames(W)[i] else i
        stop(sprintf(
          "Treatment is not absorbing for unit '%s': treated in period %d but not in all subsequent periods.",
          unit_name, t_start
        ))
      }
    }
  }

  invisible(TRUE)
}


#' Identify adoption cohorts
#'
#' Groups treated units by their first treatment period.
#'
#' @param first_treat A named numeric vector of first treatment times (as
#'   returned by \code{compute_first_treat_time}).
#' @return A named list where each element is named by the adoption time
#'   (as a character string of the column index) and contains the row indices
#'   of units in that cohort. Never-treated units are excluded.
#' @keywords internal
identify_cohorts <- function(first_treat) {
  treated <- first_treat[is.finite(first_treat)]
  if (length(treated) == 0) {
    stop("No treated units found.")
  }
  cohort_times <- sort(unique(treated))
  cohorts <- lapply(cohort_times, function(g) {
    which(first_treat == g)
  })
  names(cohorts) <- as.character(cohort_times)
  cohorts
}


#' Build a single staggered subproblem
#'
#' Constructs the (Y, N0, T0) inputs for one cohort/time-target cell, in the
#' format expected by \code{\link{synthdid_estimate}}.
#'
#' The control set for cohort \code{g} adopting at time \code{g} consists of:
#' \itemize{
#'   \item Never-treated units (always eligible)
#'   \item Not-yet-treated units: units whose first treatment is strictly after
#'     the last post-treatment period we use for cohort \code{g}
#' }
#'
#' @param Y Full outcome matrix (units x time).
#' @param W Full treatment matrix (units x time).
#' @param first_treat Named vector of first treatment times.
#' @param cohort_time Integer. The adoption time (column index) defining this
#'   cohort.
#' @param max_post Integer or \code{NULL}. Maximum number of post-treatment
#'   periods to include. If \code{NULL} (default), uses all available periods
#'   up to the next cohort's adoption or the panel end.
#' @param control_type One of \code{"not_yet_treated"} (default) or
#'   \code{"never_treated"}. Determines which units serve as controls.
#' @return A list with components:
#'   \describe{
#'     \item{Y}{Outcome submatrix with controls first, then treated cohort.}
#'     \item{N0}{Number of control units.}
#'     \item{T0}{Number of pre-treatment periods.}
#'     \item{cohort_units}{Row indices (in original Y) of treated cohort units.}
#'     \item{control_units}{Row indices (in original Y) of control units.}
#'     \item{pre_periods}{Column indices (in original Y) of pre-treatment periods.}
#'     \item{post_periods}{Column indices (in original Y) of post-treatment periods.}
#'   }
#' @keywords internal
build_staggered_subproblem <- function(Y, W, first_treat, cohort_time,
                                       max_post = NULL,
                                       control_type = c("not_yet_treated",
                                                        "never_treated")) {
  control_type <- match.arg(control_type)

  # Pre-treatment periods: all columns before cohort_time
  T0 <- cohort_time - 1
  if (T0 < 2) {
    stop(sprintf(
      "Cohort adopting at period %d has fewer than 2 pre-treatment periods.",
      cohort_time
    ))
  }

  # Determine post-treatment periods
  T_total <- ncol(Y)
  if (is.null(max_post)) {
    # Use all remaining periods
    post_end <- T_total
  } else {
    post_end <- min(cohort_time + max_post - 1, T_total)
  }
  T1 <- post_end - T0
  if (T1 < 1) {
    stop(sprintf(
      "Cohort adopting at period %d has no post-treatment periods.",
      cohort_time
    ))
  }

  pre_periods <- 1:T0
  post_periods <- (T0 + 1):post_end

  # Identify treated units in this cohort
  cohort_units <- which(first_treat == cohort_time)
  if (length(cohort_units) == 0) {
    stop(sprintf("No units found adopting at time %d.", cohort_time))
  }

  # Identify control units
  if (control_type == "never_treated") {
    control_units <- which(is.infinite(first_treat))
  } else {
    # Not-yet-treated: never-treated + those adopting strictly after post_end
    control_units <- which(first_treat > post_end)
  }

  if (length(control_units) < 2) {
    stop(sprintf(
      "Cohort at time %d has fewer than 2 eligible control units (control_type = '%s').",
      cohort_time, control_type
    ))
  }

  # Build subproblem matrices: controls first, then treated
  all_units <- c(control_units, cohort_units)
  all_periods <- c(pre_periods, post_periods)

  Y_sub <- Y[all_units, all_periods, drop = FALSE]
  N0_sub <- length(control_units)

  # Preserve dimnames
  if (!is.null(rownames(Y))) {
    rownames(Y_sub) <- rownames(Y)[all_units]
  }
  if (!is.null(colnames(Y))) {
    colnames(Y_sub) <- colnames(Y)[all_periods]
  }

  list(
    Y = Y_sub,
    N0 = N0_sub,
    T0 = T0,
    cohort_units = cohort_units,
    control_units = control_units,
    pre_periods = pre_periods,
    post_periods = post_periods
  )
}


#' Enumerate all cohort subproblems for a staggered panel
#'
#' Iterates over adoption cohorts and builds the corresponding subproblems.
#' Cohorts that cannot form valid subproblems (too few controls or periods) are
#' skipped with a message.
#'
#' @param Y Outcome matrix (units x time).
#' @param W Treatment matrix (units x time).
#' @param max_post Maximum post-treatment periods per cohort, or \code{NULL}.
#' @param control_type Control set rule: \code{"not_yet_treated"} or
#'   \code{"never_treated"}.
#' @return A named list of subproblem lists (as returned by
#'   \code{build_staggered_subproblem}), keyed by cohort adoption time.
#' @keywords internal
enumerate_staggered_subproblems <- function(Y, W,
                                            max_post = NULL,
                                            control_type = "not_yet_treated") {
  validate_staggered_panel(Y, W)
  first_treat <- compute_first_treat_time(W)
  cohorts <- identify_cohorts(first_treat)

  subproblems <- list()
  for (g_char in names(cohorts)) {
    g <- as.integer(g_char)
    sub <- tryCatch(
      build_staggered_subproblem(Y, W, first_treat, g,
                                 max_post = max_post,
                                 control_type = control_type),
      error = function(e) {
        message(sprintf("Skipping cohort at time %d: %s", g, e$message))
        NULL
      }
    )
    if (!is.null(sub)) {
      subproblems[[g_char]] <- sub
    }
  }

  if (length(subproblems) == 0) {
    stop("No valid cohort subproblems could be constructed. ",
         "Check that the panel has sufficient controls and pre-treatment periods.")
  }

  subproblems
}
