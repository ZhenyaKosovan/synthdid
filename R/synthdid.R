#' A function mapping a numeric vector to a (presumably sparser) numeric vector of the same shape to
#' be passed onto synthdid_estimate.
#' @param v a vector
sparsify_function <- function(v) {
  v[v <= max(v) / SYNTHDID_SPARSIFY_DIVISOR] <- 0
  v / sum(v)
}

#' Computes the synthetic diff-in-diff estimate for an average treatment effect on a treated block.
#'
#' See 'Synthetic Difference in Differences' by Arkhangelsky et al. This implements Algorithm 1.
#' @param Y the observation matrix.
#' @param N0 the number of control units (N_co in the paper). Rows 1-N0 of Y correspond to the control units.
#' @param T0 the number of pre-treatment time steps (T_pre in the paper). Columns 1-T0 of Y correspond to pre-treatment time steps.
#' @param X an optional 3-D array of time-varying covariates. Shape should be N X T X C for C covariates.
#' @param noise.level, an estimate of the noise standard deviation sigma. Defaults to the standard deviation of first differences of Y.
#' @param eta.omega  determines the tuning parameter zeta.omega = eta.omega * noise.level. Defaults to the value (N_tr T_post)^(1/4).
#' @param eta.lambda analogous for lambda.  Defaults to an 'infinitesimal' value 1e-6.
#' @param zeta.omega if passed, overrides the default zeta.omega = eta.omega * noise.level. Deprecated.
#' @param zeta.lambda analogous for lambda.
#' @param omega.intercept Binary. Use an intercept when estimating omega.
#' @param lambda.intercept Binary. Use an intercept when estimating lambda.
#' @param weights a list with fields lambda and omega. If non-null weights$lambda is passed,
#'        we use them instead of estimating lambda weights. Same for weights$omega.
#' @param update.omega If true, solve for omega using the passed value of weights$omega only as an initialization.
#'        If false, use it exactly as passed. Defaults to false if a non-null value of weights$omega is passed.
#' @param update.lambda  Analogous.
#' @param min.decrease Tunes a stopping criterion for our weight estimator. Stop after an iteration results in a decrease
#'   in penalized MSE smaller than min.decrease^2.
#' @param max.iter A fallback stopping criterion for our weight estimator. Stop after this number of iterations.
#' @param sparsify A function mapping a numeric vector to a (presumably sparser) numeric vector of the same shape, which must sum to one.
#'   If not null, we try to estimate sparse weights via a second round of Frank-Wolfe optimization initialized at sparsify(the solution to the first round).
#' @param max.iter.pre.sparsify Analogous to max.iter, but for the pre-sparsification first-round of optimization.
#'   Not used if sparsify = NULL.
#' @param estimate_se Logical. If TRUE, attempt to compute a standard error for the estimate (stored as an attribute).
#' @param se_method Standard-error method to use when estimate_se = TRUE; one of "bootstrap", "jackknife", or "placebo".
#' @param se_replications Number of replications when using bootstrap or placebo standard errors.
#' @return An average treatment effect estimate with 'weights' and 'setup' attached as attributes.
#'   'weights' contains the estimated weights lambda and omega and corresponding intercepts, as well as regression coefficients beta if X is passed.
#'   'setup' is a list describing the problem passed in: Y, N0, T0, X.
#'   If estimate_se = TRUE, attributes 'se', 'se_method', and 'se_status' reflect the requested standard error computation.
#' @examples
#' \donttest{
#' # Estimate treatment effect using California Proposition 99 data
#' data(california_prop99)
#' setup <- panel.matrices(california_prop99)
#'
#' # Basic SynthDID estimate
#' tau.hat <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
#' print(tau.hat)
#'
#' # With standard error computation
#' tau.hat.se <- synthdid_estimate(setup$Y, setup$N0, setup$T0,
#'   estimate_se = TRUE, se_method = "jackknife"
#' )
#' sqrt(vcov(tau.hat.se))
#'
#' # Using covariates (if available)
#' # Suppose X is a 3D array of time-varying covariates
#' # tau.hat.cov <- synthdid_estimate(setup$Y, setup$N0, setup$T0, X = X)
#'
#' # Access weights and check convergence
#' omega <- attr(tau.hat, "weights")$omega
#' lambda <- attr(tau.hat, "weights")$lambda
#' synthdid_converged(tau.hat)
#'
#' # Using custom regularization
#' tau.hat.custom <- synthdid_estimate(setup$Y, setup$N0, setup$T0,
#'   eta.omega = 0.5, eta.lambda = 1e-3
#' )
#' }
#' @export synthdid_estimate
synthdid_estimate <- function(Y, N0, T0, X = array(dim = c(dim(Y), 0)),
                              noise.level = NULL,
                              eta.omega = ((nrow(Y) - N0) * (ncol(Y) - T0))^(1 / 4),
                              eta.lambda = SYNTHDID_ETA_LAMBDA_DEFAULT,
                              zeta.omega = NULL, zeta.lambda = NULL,
                              omega.intercept = TRUE, lambda.intercept = TRUE,
                              weights = list(omega = NULL, lambda = NULL),
                              update.omega = is.null(weights$omega), update.lambda = is.null(weights$lambda),
                              min.decrease = NULL, max.iter = SYNTHDID_MAX_ITER_DEFAULT,
                              sparsify = sparsify_function,
                              max.iter.pre.sparsify = SYNTHDID_MAX_ITER_PRE_SPARSIFY,
                              estimate_se = FALSE,
                              se_method = c("bootstrap", "jackknife", "placebo"),
                              se_replications = SYNTHDID_SE_REPLICATIONS_DEFAULT) {
  if (any(is.na(Y))) {
    stop("Missing values in input data.")
  }
  if (any(is.infinite(Y))) {
    stop("Infinite values in input dat.a")
  }
  if (max.iter < 0) {
    stop("max.iter should be positive scalar.")
  }

  stopifnot(
    nrow(Y) > N0, ncol(Y) > T0, length(dim(X)) %in% c(2, 3), dim(X)[1:2] == dim(Y), is.list(weights),
    is.null(weights$lambda) || length(weights$lambda) == T0, is.null(weights$omega) || length(weights$omega) == N0,
    !is.null(weights$lambda) || update.lambda, !is.null(weights$omega) || update.omega
  )
  if (length(dim(X)) == 2) {
    dim(X) <- c(dim(X), 1)
  }
  if (is.null(noise.level)) {
    if (T0 < 2 || N0 < 1) {
      stop("noise.level is undefined with fewer than 2 pre-treatment periods.")
    } else {
      diffs <- Y[1:N0, 2:T0, drop = FALSE] - Y[1:N0, 1:(T0 - 1), drop = FALSE]
      noise.level <- stats::sd(c(diffs))
      if (!is.finite(noise.level)) {
        warning("noise.level could not be computed; using 0.")
        noise.level <- 0
      }
    }
  } else if (!is.finite(noise.level) || length(noise.level) != 1) {
    stop("noise.level must be a finite scalar.")
  } else if (noise.level < 0) {
    stop("noise.level must be a positive scalar.")
  }
  if (is.null(zeta.omega)) {
    zeta.omega <- eta.omega * noise.level
  }
  if (is.null(zeta.lambda)) {
    zeta.lambda <- eta.lambda * noise.level
  }
  if (is.null(min.decrease)) {
    min.decrease <- SYNTHDID_MIN_DECREASE_NOISE_MULTIPLIER * noise.level
  }
  if (is.null(sparsify)) {
    max.iter.pre.sparsify <- max.iter
  }
  N1 <- nrow(Y) - N0
  T1 <- ncol(Y) - T0

  # Initialize convergence tracking structure (zero overhead)
  convergence_info <- list(
    lambda = NULL,
    omega = NULL,
    overall_converged = TRUE
  )

  if (dim(X)[3] == 0) {
    weights$vals <- NULL
    weights$lambda.vals <- NULL
    weights$omega.vals <- NULL

    # PERFORMANCE: Compute collapsed form once if needed by both lambda and omega
    # Avoids duplicate expensive computation when both weights are updated
    Yc <- if (update.lambda || update.omega) {
      collapsed.form(Y, N0, T0)
    } else {
      NULL
    }

    if (update.lambda) {
      lambda.opt <- sc.weight.fw(Yc[1:N0, , drop = FALSE],
        zeta = zeta.lambda, intercept = lambda.intercept, lambda = weights$lambda,
        min.decrease = min.decrease, max.iter = max.iter.pre.sparsify,
        warn_not_converged = FALSE # Suppress warning, we'll aggregate later
      )

      if (!is.null(sparsify)) {
        lambda.opt <- sc.weight.fw(Yc[1:N0, , drop = FALSE],
          zeta = zeta.lambda, intercept = lambda.intercept, lambda = sparsify(lambda.opt$lambda),
          min.decrease = min.decrease, max.iter = max.iter,
          warn_not_converged = FALSE # Suppress warning, we'll aggregate later
        )
      }

      weights$lambda <- lambda.opt$lambda
      weights$lambda.vals <- lambda.opt$vals
      weights$vals <- lambda.opt$vals

      # Store convergence info (lightweight)
      convergence_info$lambda <- list(
        converged = lambda.opt$converged,
        iterations = lambda.opt$iterations,
        max_iter = max.iter
      )
      if (!lambda.opt$converged) {
        convergence_info$overall_converged <- FALSE
      }
    }
    if (update.omega) {
      omega.opt <- sc.weight.fw(t(Yc[, 1:T0, drop = FALSE]),
        zeta = zeta.omega, intercept = omega.intercept, lambda = weights$omega,
        min.decrease = min.decrease, max.iter = max.iter.pre.sparsify,
        warn_not_converged = FALSE # Suppress warning, we'll aggregate later
      )

      if (!is.null(sparsify)) {
        omega.opt <- sc.weight.fw(t(Yc[, 1:T0, drop = FALSE]),
          zeta = zeta.omega, intercept = omega.intercept, lambda = sparsify(omega.opt$lambda),
          min.decrease = min.decrease, max.iter = max.iter,
          warn_not_converged = FALSE # Suppress warning, we'll aggregate later
        )
      }
      weights$omega <- omega.opt$lambda
      weights$omega.vals <- omega.opt$vals
      if (is.null(weights$vals)) {
        weights$vals <- omega.opt$vals
      } else {
        weights$vals <- pairwise.sum.decreasing(weights$vals, omega.opt$vals)
      }

      # Store convergence info (lightweight)
      convergence_info$omega <- list(
        converged = omega.opt$converged,
        iterations = omega.opt$iterations,
        max_iter = max.iter
      )
      if (!omega.opt$converged) {
        convergence_info$overall_converged <- FALSE
      }
    }
  } else {
    Yc <- collapsed.form(Y, N0, T0)
    Xc <- apply(X, 3, function(Xi) {
      collapsed.form(Xi, N0, T0)
    })
    dim(Xc) <- c(dim(Yc), dim(X)[3])
    weights <- sc.weight.fw.covariates(Yc, Xc,
      zeta.lambda = zeta.lambda, zeta.omega = zeta.omega,
      lambda.intercept = lambda.intercept, omega.intercept = omega.intercept,
      min.decrease = min.decrease, max.iter = max.iter,
      lambda = weights$lambda, omega = weights$omega,
      update.lambda = update.lambda, update.omega = update.omega,
      warn_not_converged = FALSE # Suppress warning, we'll aggregate later
    )

    # Store convergence info from covariate solver (lightweight)
    convergence_info$joint <- list(
      converged = weights$converged,
      iterations = weights$iterations,
      max_iter = max.iter
    )
    convergence_info$overall_converged <- weights$converged
  }

  X.beta <- contract3(X, weights$beta)
  estimate <- t(c(-weights$omega, rep(1 / N1, N1))) %*% (Y - X.beta) %*% c(-weights$lambda, rep(1 / T1, T1))

  class(estimate) <- "synthdid_estimate"
  attr(estimate, "estimator") <- "synthdid_estimate"
  attr(estimate, "weights") <- weights
  # PERFORMANCE: Cache noise.level in setup for reuse in bootstrap/jackknife
  # Avoids recomputing expensive sd() in each replication
  attr(estimate, "setup") <- list(Y = Y, X = X, N0 = N0, T0 = T0, noise.level = noise.level)
  attr(estimate, "opts") <- list(
    zeta.omega = zeta.omega, zeta.lambda = zeta.lambda,
    omega.intercept = omega.intercept, lambda.intercept = lambda.intercept,
    update.omega = update.omega, update.lambda = update.lambda,
    min.decrease = min.decrease, max.iter = max.iter
  )

  # Attach convergence information (zero overhead, just metadata)
  attr(estimate, "convergence") <- convergence_info

  # Issue warning if optimization did not converge (lightweight check)
  if (!convergence_info$overall_converged) {
    warning_parts <- character()

    if (!is.null(convergence_info$lambda) && !convergence_info$lambda$converged) {
      warning_parts <- c(warning_parts, sprintf(
        "lambda weights (%d/%d iterations)",
        convergence_info$lambda$iterations,
        convergence_info$lambda$max_iter
      ))
    }
    if (!is.null(convergence_info$omega) && !convergence_info$omega$converged) {
      warning_parts <- c(warning_parts, sprintf(
        "omega weights (%d/%d iterations)",
        convergence_info$omega$iterations,
        convergence_info$omega$max_iter
      ))
    }
    if (!is.null(convergence_info$joint) && !convergence_info$joint$converged) {
      warning_parts <- c(warning_parts, sprintf(
        "joint optimization (%d/%d iterations)",
        convergence_info$joint$iterations,
        convergence_info$joint$max_iter
      ))
    }

    if (length(warning_parts) > 0) {
      warning(sprintf(
        "synthdid optimization did not converge: %s\nConsider increasing max.iter or relaxing min.decrease threshold.",
        paste(warning_parts, collapse = ", ")
      ), call. = FALSE)
    }
  }
  se_method <- match.arg(se_method)
  se_value <- NULL
  se_status <- "not_requested"
  if (estimate_se) {
    se_status <- "failed"
    se_warning_emitted <- FALSE
    se_value <- tryCatch(
      {
        if (se_method == "placebo" && N0 <= N1) {
          stop("placebo standard errors require more controls than treated units.")
        }
        if (se_method == "bootstrap" && N1 == 1) {
          stop("bootstrap standard errors require more than one treated unit.")
        }
        if (se_method == "jackknife" && (N1 == 1 ||
          (!is.null(weights$omega) && sum(weights$omega != 0) == 1 && !update.omega))) {
          stop("jackknife standard errors require more than one treated unit and at least two controls with weight.")
        }
        if (se_method == "bootstrap") {
          bootstrap_se(estimate, se_replications)
        } else if (se_method == "jackknife") {
          jackknife_se(estimate)
        } else {
          placebo_se(estimate, se_replications)
        }
      },
      error = function(e) {
        warning(e$message)
        se_warning_emitted <<- TRUE
        NA_real_
      }
    )
    if (!is.na(se_value)) {
      se_status <- "computed"
    } else if (!se_warning_emitted) {
      warning(sprintf("%s standard errors could not be computed; returning NA.", se_method))
    }
  }
  attr(estimate, "se") <- se_value
  attr(estimate, "se_method") <- if (estimate_se) se_method else NULL
  attr(estimate, "se_status") <- se_status
  return(estimate)
}

#' synthdid_estimate for synthetic control estimates.
#' Takes all the same parameters, but by default, passes options to use the synthetic control estimator
#' By default, this uses only 'infinitesimal' ridge regularization when estimating the weights.
#' @param Y the observation matrix.
#' @param N0 the number of control units. Rows 1-N0 of Y correspond to the control units.
#' @param T0 the number of pre-treatment time steps. Columns 1-T0 of Y correspond to pre-treatment time steps.
#' @param eta.omega determines the level of ridge regularization, zeta.omega = eta.omega * noise.level, as in synthdid_estimate.
#' @param ... additional options for synthdid_estimate
#' @return an object like that returned by synthdid_estimate
#' @examples
#' \donttest{
#' # Estimate treatment effect using synthetic control method
#' data(california_prop99)
#' setup <- panel.matrices(california_prop99)
#'
#' # Synthetic control estimate
#' tau.sc <- sc_estimate(setup$Y, setup$N0, setup$T0)
#' print(tau.sc)
#'
#' # Compare with SynthDID
#' tau.sdid <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
#' c(sc = tau.sc, sdid = tau.sdid)
#'
#' # With standard error
#' tau.sc.se <- sc_estimate(setup$Y, setup$N0, setup$T0,
#'   estimate_se = TRUE, se_method = "placebo"
#' )
#' sqrt(vcov(tau.sc.se))
#' }
#' @export sc_estimate
sc_estimate <- function(Y, N0, T0, eta.omega = SYNTHDID_ETA_OMEGA_SC_DEFAULT, ...) {
  estimate <- synthdid_estimate(Y, N0, T0,
    eta.omega = eta.omega,
    weights = list(lambda = rep(0, T0)), omega.intercept = FALSE, ...
  )
  attr(estimate, "estimator") <- "sc_estimate"
  estimate
}

#' synthdid_estimate for diff-in-diff estimates.
#' Takes all the same parameters, but by default, passes options to use the diff-in-diff estimator
#' @param Y the observation matrix.
#' @param N0 the number of control units. Rows 1-N0 of Y correspond to the control units.
#' @param T0 the number of pre-treatment time steps. Columns 1-T0 of Y correspond to pre-treatment time steps.
#' @param ... additional  options for synthdid_estimate
#' @return an object like that returned by synthdid_estimate
#' @examples
#' \donttest{
#' # Estimate treatment effect using difference-in-differences
#' data(california_prop99)
#' setup <- panel.matrices(california_prop99)
#'
#' # DID estimate
#' tau.did <- did_estimate(setup$Y, setup$N0, setup$T0)
#' print(tau.did)
#'
#' # Compare all three estimators
#' tau.sc <- sc_estimate(setup$Y, setup$N0, setup$T0)
#' tau.sdid <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
#' estimates <- list(did = tau.did, sc = tau.sc, sdid = tau.sdid)
#' sapply(estimates, function(x) x)
#'
#' # Visualize the differences
#' # synthdid_plot(estimates)
#' }
#' @export did_estimate
did_estimate <- function(Y, N0, T0, ...) {
  estimate <- synthdid_estimate(Y, N0, T0, weights = list(lambda = rep(1 / T0, T0), omega = rep(1 / N0, N0)), ...)
  attr(estimate, "estimator") <- "did_estimate"
  estimate
}

#' Computes a placebo variant of our estimator using pre-treatment data only
#' @param estimate, as output by synthdid_estimate
#' @param treated.fraction, the fraction of pre-treatment data to use as a placebo treatment period
#'        Defaults to NULL, which indicates that it should be the fraction of post-treatment to pre-treatment data
#' @return A placebo estimate using only pre-treatment data
#' @examples
#' \donttest{
#' data(california_prop99)
#' setup <- panel.matrices(california_prop99)
#' tau.hat <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
#'
#' # Compute placebo estimate
#' tau.placebo <- synthdid_placebo(tau.hat)
#'
#' # Should be close to zero if parallel trends holds
#' c(estimate = tau.hat, placebo = tau.placebo)
#' }
#' @export synthdid_placebo
synthdid_placebo <- function(estimate, treated.fraction = NULL) {
  setup <- attr(estimate, "setup")
  opts <- attr(estimate, "opts")
  weights <- attr(estimate, "weights")
  X.beta <- contract3(setup$X, weights$beta)
  estimator <- attr(estimate, "estimator")

  if (is.null(treated.fraction)) {
    treated.fraction <- 1 - setup$T0 / ncol(setup$Y)
  }
  placebo.T0 <- floor(setup$T0 * (1 - treated.fraction))

  do.call(estimator, c(list(Y = setup$Y[, 1:setup$T0], N0 = setup$N0, T0 = placebo.T0, X = setup$X[, 1:setup$T0, ]), opts))
}

#' Outputs the effect curve that was averaged to produce our estimate
#' @param estimate, as output by synthdid_estimate
#' @return A vector of treatment effects for each post-treatment period
#' @examples
#' \donttest{
#' data(california_prop99)
#' setup <- panel.matrices(california_prop99)
#' tau.hat <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
#'
#' # Get effect curve over time
#' effect_curve <- synthdid_effect_curve(tau.hat)
#' plot(effect_curve, type = "l", xlab = "Post-treatment period", ylab = "Effect")
#' }
#' @export synthdid_effect_curve
synthdid_effect_curve <- function(estimate) {
  setup <- attr(estimate, "setup")
  weights <- attr(estimate, "weights")
  X.beta <- contract3(setup$X, weights$beta)
  N1 <- nrow(setup$Y) - setup$N0
  T1 <- ncol(setup$Y) - setup$T0

  tau.sc <- t(c(-weights$omega, rep(1 / N1, N1))) %*% (setup$Y - X.beta)
  tau.curve <- tau.sc[setup$T0 + (1:T1)] - c(tau.sc[1:setup$T0] %*% weights$lambda)
  tau.curve
}
