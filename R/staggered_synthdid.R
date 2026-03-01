# Staggered Adoption Synthetic Difference-in-Differences
#
# Implements the Porreca (2022) style staggered SDID estimator by decomposing
# a staggered panel into cohort-specific 2x2 subproblems and aggregating.

#' Staggered Synthetic Difference-in-Differences (Matrix API)
#'
#' Estimates average treatment effects in panels with staggered adoption of
#' treatment. Decomposes the problem into cohort-specific 2x2 SDID subproblems,
#' estimates each using the existing SDID kernel, and aggregates.
#'
#' @param Y Outcome matrix (units x time). Rows are units, columns are time
#'   periods.
#' @param W Binary treatment matrix (units x time). Must be absorbing: once a
#'   unit is treated, it stays treated.
#' @param X Optional 3-D array of time-varying covariates, shaped N x T x C.
#' @param method Estimation method for each subproblem: \code{"synthdid"}
#'   (default), \code{"sc"}, or \code{"did"}.
#' @param max_post Maximum number of post-treatment periods per cohort.
#'   \code{NULL} (default) uses all available.
#' @param control_type Control set rule: \code{"not_yet_treated"} (default)
#'   uses never-treated and not-yet-treated units; \code{"never_treated"} uses
#'   only never-treated units.
#' @param weights Aggregation scheme for cohort effects:
#'   \code{"treated_unit_time"} (default) weights by number of treated
#'   unit-time cells in each cohort, \code{"cohort_size"} weights by number
#'   of treated units, \code{"equal"} gives equal weight to each cohort.
#' @param warm_start Logical. If \code{TRUE} (default), use solved weights from
#'   adjacent cohorts as initialization for the next subproblem.
#' @param ... Additional arguments passed to \code{synthdid_estimate()}.
#'
#' @return An object of class \code{"synthdid_staggered"} with components:
#'   \describe{
#'     \item{att}{Scalar aggregate ATT.}
#'     \item{cohort_effects}{Data frame of cohort-level effects with columns:
#'       cohort_time, estimate, n_treated, n_control, n_pre, n_post, weight.}
#'     \item{aggregation_weights}{Named numeric vector of weights used for aggregation.}
#'     \item{subproblem_fits}{List of per-cohort \code{synthdid_estimate} objects.}
#'     \item{method}{Estimation method used.}
#'     \item{setup}{List with Y, W, and metadata about the full panel.}
#'   }
#'
#' @examples
#' \donttest{
#' # Simulate a staggered panel
#' set.seed(42)
#' N <- 30; TT <- 20
#' Y <- matrix(rnorm(N * TT), N, TT)
#' W <- matrix(0, N, TT)
#' # Cohort 1: units 21-25 treated from period 11
#' W[21:25, 11:TT] <- 1
#' # Cohort 2: units 26-30 treated from period 15
#' W[26:30, 15:TT] <- 1
#' rownames(Y) <- rownames(W) <- paste0("unit", 1:N)
#' colnames(Y) <- colnames(W) <- 1:TT
#'
#' est <- synthdid_staggered_estimate(Y, W)
#' print(est)
#' }
#'
#' @seealso [synthdid_estimate()], [synthdid()]
#' @export
synthdid_staggered_estimate <- function(Y, W,
                                        X = NULL,
                                        method = c("synthdid", "sc", "did"),
                                        max_post = NULL,
                                        control_type = c("not_yet_treated",
                                                         "never_treated"),
                                        weights = c("treated_unit_time",
                                                     "cohort_size",
                                                     "equal"),
                                        warm_start = TRUE,
                                        ...) {
  method <- match.arg(method)
  control_type <- match.arg(control_type)
  weights_scheme <- match.arg(weights)

  # Build subproblems
  subproblems <- enumerate_staggered_subproblems(
    Y, W, max_post = max_post, control_type = control_type
  )

  # Choose the per-subproblem estimator function
  estimator_fn <- switch(method,
    "synthdid" = function(...) {
      withCallingHandlers(
        synthdid_estimate(...),
        lifecycle_warning_deprecated = function(cnd) invokeRestart("muffleWarning")
      )
    },
    "sc" = function(...) {
      withCallingHandlers(
        sc_estimate(...),
        lifecycle_warning_deprecated = function(cnd) invokeRestart("muffleWarning")
      )
    },
    "did" = function(...) {
      withCallingHandlers(
        did_estimate(...),
        lifecycle_warning_deprecated = function(cnd) invokeRestart("muffleWarning")
      )
    }
  )

  # Solve each subproblem
  fits <- list()
  prev_weights <- NULL

  for (g_char in names(subproblems)) {
    sub <- subproblems[[g_char]]

    # Build covariate subarray if X is provided
    X_sub <- if (!is.null(X)) {
      all_units <- c(sub$control_units, sub$cohort_units)
      all_periods <- c(sub$pre_periods, sub$post_periods)
      X[all_units, all_periods, , drop = FALSE]
    } else {
      array(dim = c(dim(sub$Y), 0))
    }

    # Warm-start: use previous cohort's weights as initialization
    extra_args <- list(...)
    if (warm_start && !is.null(prev_weights) && method != "did") {
      # Only warm-start if dimensions are compatible
      prev_omega <- prev_weights$omega
      prev_lambda <- prev_weights$lambda

      if (length(prev_lambda) == sub$T0) {
        extra_args$weights <- list(
          lambda = prev_lambda,
          omega = NULL  # omega dimensions change across cohorts
        )
        extra_args$update.lambda <- TRUE
      }
    }

    fit <- do.call(estimator_fn, c(
      list(
        Y = sub$Y,
        N0 = sub$N0,
        T0 = sub$T0,
        X = X_sub
      ),
      extra_args
    ))

    fits[[g_char]] <- fit
    prev_weights <- attr(fit, "weights")
  }

  # Compute aggregation weights
  agg_weights <- compute_aggregation_weights(subproblems, weights_scheme)

  # Compute aggregate ATT
  cohort_estimates <- vapply(fits, function(f) c(f), numeric(1))
  att <- sum(cohort_estimates * agg_weights)

  # Build cohort effects table
  cohort_effects <- data.frame(
    cohort_time = as.integer(names(subproblems)),
    estimate = cohort_estimates,
    n_treated = vapply(subproblems, function(s) length(s$cohort_units), numeric(1)),
    n_control = vapply(subproblems, function(s) s$N0, numeric(1)),
    n_pre = vapply(subproblems, function(s) s$T0, numeric(1)),
    n_post = vapply(subproblems, function(s) length(s$post_periods), numeric(1)),
    weight = agg_weights,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  # Assemble result
  result <- list(
    att = att,
    cohort_effects = cohort_effects,
    aggregation_weights = agg_weights,
    subproblem_fits = fits,
    method = method,
    setup = list(
      Y = Y,
      W = W,
      X = X,
      n_units = nrow(Y),
      n_periods = ncol(Y),
      control_type = control_type,
      weights_scheme = weights_scheme,
      max_post = max_post
    )
  )

  class(result) <- "synthdid_staggered"
  attr(result, "se") <- NULL
  attr(result, "se_method") <- NULL

  result
}


#' Compute aggregation weights for staggered cohorts
#'
#' @param subproblems List of subproblem specifications (from
#'   \code{enumerate_staggered_subproblems}).
#' @param scheme Weighting scheme: \code{"treated_unit_time"},
#'   \code{"cohort_size"}, or \code{"equal"}.
#' @return Named numeric vector of weights summing to 1.
#' @keywords internal
compute_aggregation_weights <- function(subproblems, scheme) {
  raw <- switch(scheme,
    "treated_unit_time" = vapply(subproblems, function(s) {
      length(s$cohort_units) * length(s$post_periods)
    }, numeric(1)),
    "cohort_size" = vapply(subproblems, function(s) {
      length(s$cohort_units)
    }, numeric(1)),
    "equal" = rep(1, length(subproblems))
  )
  names(raw) <- names(subproblems)
  raw / sum(raw)
}
