# Inference for Staggered SDID Objects
#
# Provides vcov, confint, and standard error computation for
# synthdid_staggered objects.

#' Variance-Covariance for Staggered SDID
#'
#' Computes variance estimates for the aggregate ATT from a staggered SDID
#' object. Supports bootstrap (unit-level resampling across cohorts) and
#' placebo methods.
#'
#' @param object A \code{synthdid_staggered} object.
#' @param method Inference method: \code{"bootstrap"} (default) or
#'   \code{"placebo"}.
#' @param replications Number of bootstrap/placebo replications.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A 1x1 matrix containing the variance estimate.
#'
#' @method vcov synthdid_staggered
#' @export
vcov.synthdid_staggered <- function(object,
                                    method = c("bootstrap", "placebo"),
                                    replications = SYNTHDID_SE_REPLICATIONS_DEFAULT,
                                    ...) {
  method <- match.arg(method)

  # Check for cached SE
  precomputed_se <- attr(object, "se")
  precomputed_method <- attr(object, "se_method")
  if (!is.null(precomputed_se) && !is.null(precomputed_method) &&
      identical(precomputed_method, method)) {
    return(matrix(precomputed_se^2))
  }

  se <- staggered_bootstrap_se(object, replications)
  matrix(se^2)
}


#' Bootstrap SE for staggered SDID
#'
#' Resamples entire units (preserving cohort structure) and re-estimates the
#' staggered ATT for each bootstrap draw.
#'
#' @param object A \code{synthdid_staggered} object.
#' @param replications Number of bootstrap draws.
#' @return Scalar standard error estimate.
#' @keywords internal
staggered_bootstrap_se <- function(object, replications) {
  setup <- object$setup
  Y <- setup$Y
  W <- setup$W
  X <- setup$X
  first_treat <- compute_first_treat_time(W)

  # Identify never-treated and each cohort
  never_treated <- which(is.infinite(first_treat))
  cohorts <- identify_cohorts(first_treat)
  n_never <- length(never_treated)

  boot_estimates <- with_blas_thread_management({
    furrr::future_map_dbl(
      seq_len(replications),
      function(b) {
        # Resample within groups: never-treated and each cohort
        boot_units <- integer(0)
        # Resample never-treated
        if (n_never > 0) {
          boot_units <- c(boot_units, sample(never_treated, replace = TRUE))
        }
        # Resample within each cohort
        for (cohort_units in cohorts) {
          boot_units <- c(boot_units, sample(cohort_units, replace = TRUE))
        }
        boot_units <- sort(boot_units)

        Y_boot <- Y[boot_units, , drop = FALSE]
        W_boot <- W[boot_units, , drop = FALSE]
        X_boot <- if (!is.null(X)) X[boot_units, , , drop = FALSE] else NULL

        # Need unique rownames for panel operations
        rownames(Y_boot) <- paste0("u", seq_len(nrow(Y_boot)))
        rownames(W_boot) <- rownames(Y_boot)

        tryCatch({
          est <- synthdid_staggered_estimate(
            Y_boot, W_boot,
            X = X_boot,
            method = object$method,
            max_post = setup$max_post,
            control_type = setup$control_type,
            weights = setup$weights_scheme,
            warm_start = TRUE
          )
          est$att
        }, error = function(e) NA_real_)
      },
      .options = furrr::furrr_options(seed = TRUE)
    )
  })

  # Drop failed replications
  boot_estimates <- boot_estimates[!is.na(boot_estimates)]
  if (length(boot_estimates) < 10) {
    warning("Fewer than 10 bootstrap replications succeeded. SE may be unreliable.")
  }
  if (length(boot_estimates) == 0) {
    return(NA_real_)
  }

  sqrt((length(boot_estimates) - 1) / length(boot_estimates)) * sd(boot_estimates)
}
