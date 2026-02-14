#' Calculate Variance-Covariance Matrix for a Fitted Model Object
#'
#' Provides variance estimates based on the following three options
#' \itemize{
#'   \item The bootstrap, Algorithm 2 in Arkhangelsky et al.
#'   \item The jackknife, Algorithm 3 in Arkhangelsky et al.
#'   \item Placebo, Algorithm 4 in Arkhangelsky et al.
#' }
#'
#' The jackknife is not recommended for SC, see section 5 in Arkhangelsky et al.
#' "placebo" is the only option that works for only one treated unit.
#'
#' @param object A synthdid model
#' @param method, the CI method. The default is bootstrap (warning: this may be slow on large
#'  data sets, the jackknife option is the fastest, with the caveat that it is not recommended
#'  for SC).
#' @param replications, the number of bootstrap replications
#' @param ... Additional arguments (currently ignored).
#'
#' @return A 1x1 matrix containing the variance estimate.
#' @examples
#' \donttest{
#' # Compute variance using different methods
#' data(california_prop99)
#' setup <- panel.matrices(california_prop99)
#' tau.hat <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
#'
#' # Bootstrap standard error (default, may be slow)
#' se.bootstrap <- sqrt(vcov(tau.hat, method = "bootstrap", replications = 100))
#'
#' # Jackknife standard error (faster)
#' se.jackknife <- sqrt(vcov(tau.hat, method = "jackknife"))
#'
#' # Placebo standard error
#' se.placebo <- sqrt(vcov(tau.hat, method = "placebo", replications = 100))
#'
#' # Display standard errors
#' c(bootstrap = se.bootstrap, jackknife = se.jackknife, placebo = se.placebo)
#'
#' # Or compute SE at estimation time
#' tau.hat.se <- synthdid_estimate(setup$Y, setup$N0, setup$T0,
#'   estimate_se = TRUE, se_method = "jackknife"
#' )
#' # Reuses cached SE without recomputation
#' se.cached <- sqrt(vcov(tau.hat.se, method = "jackknife"))
#' }
#' @references Dmitry Arkhangelsky, Susan Athey, David A. Hirshberg, Guido W. Imbens, and Stefan Wager.
#'  "Synthetic Difference in Differences". arXiv preprint arXiv:1812.09970, 2019.
#'
#' @method vcov synthdid_estimate
#' @export
vcov.synthdid_estimate <- function(
    object,
    method = c("bootstrap", "jackknife", "placebo"),
    replications = SYNTHDID_SE_REPLICATIONS_DEFAULT, ...) {
  method <- match.arg(method)
  precomputed_se <- attr(object, "se")
  precomputed_method <- attr(object, "se_method")
  if (!is.null(precomputed_se) && !is.null(precomputed_method) && identical(precomputed_method, method)) {
    return(matrix(precomputed_se^2))
  }
  if (method == "bootstrap") {
    se <- bootstrap_se(object, replications)
  } else if (method == "jackknife") {
    se <- jackknife_se(object)
  } else if (method == "placebo") {
    se <- placebo_se(object, replications)
  }
  matrix(se^2)
}

#' Calculate the standard error of a synthetic diff in diff estimate. Deprecated. Use vcov.synthdid_estimate.
#' @param ... Any valid arguments for vcov.synthdid_estimate
#' @return A scalar standard error estimate.
#' @examples
#' \donttest{
#' # This function is deprecated. Use sqrt(vcov(...)) instead.
#' data(california_prop99)
#' setup <- panel.matrices(california_prop99)
#' tau.hat <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
#'
#' # Deprecated approach (still works)
#' se.old <- synthdid_se(tau.hat, method = "jackknife")
#'
#' # Preferred approach
#' se.new <- sqrt(vcov(tau.hat, method = "jackknife"))
#'
#' # Both give the same result
#' all.equal(se.old, se.new)
#' }
#' @export synthdid_se
synthdid_se <- function(...) {
  sqrt(vcov(...))
}


#' Bootstrap standard error for synthdid estimates
#'
#' Implements Algorithm 2 from Arkhangelsky et al. to compute a bootstrap
#' standard error.
#'
#' @param estimate A \code{synthdid_estimate} object.
#' @param replications Number of bootstrap replications.
#' @return Scalar bootstrap standard error.
#' @keywords internal
bootstrap_se <- function(estimate, replications) {
  setup <- attr(estimate, "setup")
  opts <- attr(estimate, "opts")
  weights <- attr(estimate, "weights")
  if (setup$N0 == nrow(setup$Y) - 1) {
    return(NA)
  }

  # Use BLAS thread management to prevent oversubscription with parallel workers
  bootstrap_estimates <- with_blas_thread_management({
    furrr::future_map_dbl(
      seq_len(replications),
      \(x) {
        # repeat ensures we discard draws that pick only control or only treated units
        repeat {
          ind <- sample.int(nrow(setup$Y), replace = TRUE)
          if (all(ind <= setup$N0) || all(ind > setup$N0)) {
            next
          }
          sorted_ind <- sort(ind)
          weights.boot <- weights
          weights.boot$omega <- sum_normalize(weights$omega[sorted_ind[sorted_ind <= setup$N0]])

          return(as.vector(synthdid_estimate(
            Y = setup$Y[sorted_ind, ],
            N0 = sum(sorted_ind <= setup$N0),
            T0 = setup$T0,
            X = setup$X[sorted_ind, , ],
            weights = weights.boot,
            zeta.omega = opts$zeta.omega,
            zeta.lambda = opts$zeta.lambda,
            omega.intercept = opts$omega.intercept,
            lambda.intercept = opts$lambda.intercept,
            update.omega = opts$update.omega,
            update.lambda = opts$update.lambda,
            min.decrease = opts$min.decrease,
            max.iter = opts$max.iter,
            suppress_convergence_warning = opts$suppress_convergence_warning
          )))
        }
      },
      .options = furrr::furrr_options(seed = TRUE)
    )
  })

  sqrt((replications - 1) / replications) * sd(bootstrap_estimates)
}


#' Jackknife standard error for synthdid estimates
#'
#' Implements Algorithm 3 from Arkhangelsky et al. Supports both fixed-weight
#' and refit variants depending on whether weights are supplied.
#'
#' @param estimate A \code{synthdid_estimate} object.
#' @param weights Optional weight list; if provided, updates are suppressed to
#'   emulate the fixed-weight jackknife.show_convergence_warning
#' @return Scalar jackknife standard error, or \code{NA} when not applicable.
#' @keywords internal
jackknife_se <- function(estimate, weights = attr(estimate, "weights")) {
  setup <- attr(estimate, "setup")
  opts <- attr(estimate, "opts")
  if (!is.null(weights)) {
    opts$update.omega <- opts$update.lambda <- FALSE
  }
  if (setup$N0 == nrow(setup$Y) - 1 || (!is.null(weights) && sum(weights$omega != 0) == 1)) {
    return(NA)
  }

  # Use BLAS thread management to prevent oversubscription with parallel workers
  jackknife_estimates <- with_blas_thread_management({
    furrr::future_map_dbl(
      seq_len(nrow(setup$Y)),
      \(i) {
        ind <- seq_len(nrow(setup$Y))[-i]
        weights.jk <- weights
        if (!is.null(weights)) {
          weights.jk$omega <- sum_normalize(weights$omega[ind[ind <= setup$N0]])
        }
        as.vector(synthdid_estimate(
          Y = setup$Y[ind, ],
          N0 = sum(ind <= setup$N0),
          T0 = setup$T0,
          X = setup$X[ind, , ],
          weights = weights.jk,
          zeta.omega = opts$zeta.omega,
          zeta.lambda = opts$zeta.lambda,
          omega.intercept = opts$omega.intercept,
          lambda.intercept = opts$lambda.intercept,
          update.omega = opts$update.omega,
          update.lambda = opts$update.lambda,
          min.decrease = opts$min.decrease,
          max.iter = opts$max.iter,
          suppress_convergence_warning = opts$suppress_convergence_warning
        ))
      },
      .options = furrr::furrr_options(seed = TRUE)
    )
  })

  n <- nrow(setup$Y)
  sqrt(((n - 1) / n) * (n - 1) * stats::var(jackknife_estimates))
}


#' Placebo standard error for synthdid estimates
#'
#' Implements Algorithm 4 from Arkhangelsky et al. by repeatedly treating random
#' subsets of control units as pseudo-treated to approximate the estimator's
#' variability.
#'
#' @param estimate A \code{synthdid_estimate} object.
#' @param replications Number of placebo replications.
#'
#' @return Scalar placebo standard error.
#' @keywords internal
placebo_se <- function(estimate, replications) {
  setup <- attr(estimate, "setup")
  opts <- attr(estimate, "opts")
  weights <- attr(estimate, "weights")
  N1 <- nrow(setup$Y) - setup$N0
  if (setup$N0 <= N1) {
    stop("must have more controls than treated units to use the placebo se")
  }

  # Use BLAS thread management to prevent oversubscription with parallel workers
  placebo_estimates <- with_blas_thread_management({
    furrr::future_map_dbl(
      seq_len(replications),
      \(x) {
        ind <- sample.int(setup$N0)
        N0 <- length(ind) - N1
        weights.boot <- weights
        weights.boot$omega <- sum_normalize(weights$omega[ind[1:N0]])
        as.vector(synthdid_estimate(
          Y = setup$Y[ind, ],
          N0 = N0,
          T0 = setup$T0,
          X = setup$X[ind, , ],
          weights = weights.boot,
          zeta.omega = opts$zeta.omega,
          zeta.lambda = opts$zeta.lambda,
          omega.intercept = opts$omega.intercept,
          lambda.intercept = opts$lambda.intercept,
          update.omega = opts$update.omega,
          update.lambda = opts$update.lambda,
          min.decrease = opts$min.decrease,
          max.iter = opts$max.iter,
          suppress_convergence_warning = opts$suppress_convergence_warning
        ))
      },
      .options = furrr::furrr_options(seed = TRUE)
    )
  })

  sqrt((replications - 1) / replications) * sd(placebo_estimates)
}

#' Normalize a weight vector to sum to one
#'
#' Converts an input vector into a simple probability distribution. If all
#' entries are zero, returns a uniform vector of the same length.
#'
#' @param x Numeric vector.
#'
#' @return A numeric vector with entries summing to one.
#' @keywords internal
sum_normalize <- function(x) {
  if (sum(x) != 0) {
    x / sum(x)
  } else {
    # if given a vector of zeros, return uniform weights
    # this is fine when used in bootstrap and placebo standard errors, where it is used only for initialization
    # for jackknife standard errors, where it isn't, we handle the case of a vector of zeros without calling this function.
    rep(1 / length(x), length(x))
  }
}
