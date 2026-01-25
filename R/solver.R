#' Contract a 3-D array along its covariate dimension
#'
#' Multiplies each slice \code{X[,,k]} by the corresponding coefficient in
#' \code{v} and sums the results. When \code{v} is empty, returns a zero matrix
#' of the appropriate shape.
#'
#' @param X A three-dimensional array of covariates shaped \code{N x T x C}.
#' @param v Numeric vector of length \code{C} containing coefficients.
#'
#' @return A numeric matrix with dimension \code{dim(X)[1:2]}.
#' @keywords internal
contract3 <- function(X, v) {
  stopifnot(length(dim(X)) == 3, dim(X)[3] == length(v))
  if (length(v) == 0) {
    return(array(0, dim = dim(X)[1:2]))
  }
  contract3_cpp(X, v)
}

#' Frank-Wolfe solver for synthetic control time weights
#'
#' Estimates \code{lambda} weights using exact line search under simplex
#' constraints with optional ridge regularization.
#'
#' @param Y Collapsed outcome matrix (controls first, pre-treatment first).
#' @param zeta Ridge penalty multiplier.
#' @param intercept Logical; include an intercept adjustment when fitting
#'   weights.
#' @param lambda Optional initialization for \code{lambda}; defaults to uniform.
#' @param min.decrease Minimum squared-improvement threshold to continue
#'   iterating.
#' @param max.iter Maximum number of Frank-Wolfe iterations.
#'
#' @return A list with fields \code{lambda} (weights), \code{vals} (objective
#'   trace), \code{converged} (boolean), and \code{iterations} (count).
#' @keywords internal
sc.weight.fw <- function(Y,
                         zeta,
                         intercept = TRUE,
                         lambda = NULL,
                         min.decrease = SYNTHDID_MIN_DECREASE_DEFAULT,
                         max.iter = SYNTHDID_MAX_ITER_DEFAULT,
                         warn_not_converged = TRUE) {
  T0 <- ncol(Y) - 1
  if (is.null(lambda)) {
    lambda <- rep(1 / T0, T0)
  }
  result <- sc_weight_fw_cpp(Y, zeta, intercept, lambda, min.decrease, max.iter)

  # Issue warning if optimization did not converge (lightweight check)
  if (warn_not_converged && !result$converged) {
    # Compute final decrease for diagnostic message (uses existing vals)
    final_idx <- result$iterations
    final_decrease <- if (final_idx >= 2 && !is.na(result$vals[final_idx]) && !is.na(result$vals[final_idx - 1])) {
      result$vals[final_idx - 1] - result$vals[final_idx]
    } else {
      NA_real_
    }

    warning(sprintf(
      "Frank-Wolfe optimization did not converge within %d iterations (final decrease: %.2e, threshold: %.2e)",
      max.iter, final_decrease, min.decrease^2
    ), call. = FALSE)
  }

  result
}

#' Joint Frank-Wolfe and gradient solver with covariates
#'
#' Alternates Frank-Wolfe updates for \code{lambda} and \code{omega} with
#' gradient steps for covariate coefficients \code{beta}. Supports holding
#' either weight vector fixed by setting \code{update.lambda} or
#' \code{update.omega} to \code{FALSE}.
#'
#' @param Y Collapsed outcome matrix (controls first, pre-treatment first).
#' @param X Array of covariates aligned with \code{Y}, shaped \code{N x T x C}.
#' @param zeta.lambda,zeta.omega Ridge penalties for \code{lambda} and
#'   \code{omega}.
#' @param lambda.intercept,omega.intercept Logical flags for whether to demean
#'   when fitting each set of weights.
#' @param min.decrease Minimum squared-improvement threshold to continue
#'   iterating.
#' @param max.iter Maximum number of iterations.
#' @param lambda,omega,beta Optional initial values; defaults are uniform
#'   weights and zero coefficients.
#' @param update.lambda,update.omega Logical flags indicating which weights to
#'   update.
#'
#' @return A list with estimated \code{lambda}, \code{omega}, \code{beta},
#'   the vector of objective values \code{vals}, convergence flag
#'   \code{converged}, and iteration count \code{iterations}.
#' @keywords internal
sc.weight.fw.covariates <- function(Y,
                                    X = array(0, dim = c(dim(Y), 0)),
                                    zeta.lambda = 0,
                                    zeta.omega = 0,
                                    lambda.intercept = TRUE,
                                    omega.intercept = TRUE,
                                    min.decrease = SYNTHDID_MIN_DECREASE_DEFAULT,
                                    max.iter = SYNTHDID_MAX_ITER_DEFAULT,
                                    lambda = NULL,
                                    omega = NULL,
                                    beta = NULL,
                                    update.lambda = TRUE,
                                    update.omega = TRUE,
                                    warn_not_converged = TRUE) {
  stopifnot(
    length(dim(Y)) == 2,
    length(dim(X)) == 3,
    all(dim(Y) == dim(X)[1:2]),
    all(is.finite(Y)),
    all(is.finite(X))
  )
  T0 <- ncol(Y) - 1
  N0 <- nrow(Y) - 1
  if (length(dim(X)) == 2) {
    dim(X) <- c(dim(X), 1)
  }
  if (is.null(lambda)) {
    lambda <- rep(1 / T0, T0)
  }
  if (is.null(omega)) {
    omega <- rep(1 / N0, N0)
  }
  if (is.null(beta)) {
    beta <- rep(0, dim(X)[3])
  }



  vals <- rep(NA, max.iter)
  t <- 0
  converged <- FALSE
  Y.beta <- Y - contract3(X, beta)
  weights <- update.weights(
    Y.beta,
    lambda, omega,
    lambda.intercept,
    omega.intercept,
    N0,
    T0,
    update.lambda,
    update.omega,
    zeta.lambda,
    zeta.omega
  )
  # state is kept in weights$lambda, weights$omega, beta
  while (t < max.iter && (t < 2 || abs(vals[t - 1] - vals[t]) > min.decrease^2)) {
    t <- t + 1
    grad.beta <- -if (dim(X)[3] == 0) {
      c()
    } else {
      grad_beta_cpp(
        X,
        weights$lambda,
        weights$omega,
        weights$err.lambda,
        weights$err.omega,
        N0,
        T0
      )
    }

    alpha <- 1 / t
    beta <- beta - alpha * grad.beta
    Y.beta <- Y - contract3(X, beta)
    weights <- update.weights(
      Y.beta, weights$lambda, weights$omega,
      lambda.intercept,
      omega.intercept,
      N0,
      T0,
      update.lambda,
      update.omega,
      zeta.lambda,
      zeta.omega
    )
    vals[t] <- weights$val

    # Check if converged (hit stopping criterion before max_iter)
    if (t >= 2 && abs(vals[t - 1] - vals[t]) <= min.decrease^2) {
      converged <- TRUE
    }
  }

  # Issue warning if not converged (lightweight check)
  if (warn_not_converged && !converged && t >= max.iter) {
    final_decrease <- if (t >= 2 && !is.na(vals[t]) && !is.na(vals[t - 1])) {
      abs(vals[t - 1] - vals[t])
    } else {
      NA_real_
    }

    warning(sprintf(
      "Covariate optimization did not converge within %d iterations (final decrease: %.2e, threshold: %.2e)",
      max.iter, final_decrease, min.decrease^2
    ), call. = FALSE)
  }

  list(
    lambda = weights$lambda,
    omega = weights$omega,
    beta = beta,
    vals = vals,
    converged = converged,
    iterations = t
  )
}


#' Update weight vectors for synthdid optimization
#'
#' Performs a single round of Frank-Wolfe updates for \code{lambda} and
#' \code{omega}, returning residuals and objective values used to monitor
#' convergence.
#'
#' @param Y Collapsed outcome matrix (controls first, pre-treatment first).
#' @param lambda,omega Current weight vectors.
#' @param lambda.intercept,omega.intercept Logical flags for whether to demean
#'   when fitting each set of weights.
#' @param N0,T0 Numbers of control units and pre-treatment periods.
#' @param update.lambda,update.omega Logical flags indicating which weights to
#'   update.
#' @param zeta.lambda,zeta.omega Ridge penalties for \code{lambda} and
#'   \code{omega}.
#'
#' @return A list with updated \code{lambda}, \code{omega}, residual vectors
#'   \code{err.lambda}, \code{err.omega}, and the objective value \code{val}.
#' @keywords internal
update.weights <- function(Y,
                           lambda,
                           omega,
                           lambda.intercept,
                           omega.intercept,
                           N0,
                           T0,
                           update.lambda,
  update.omega,
  zeta.lambda,
  zeta.omega) {
  Y.lambda <- if (lambda.intercept) {
    sweep(Y[1:N0, , drop = FALSE], 2, colMeans(Y[1:N0, , drop = FALSE]))
  } else {
    Y[1:N0, , drop = FALSE]
  }
  if (update.lambda) {
    lambda <- fw_step_cpp(
      Y.lambda[, 1:T0, drop = FALSE],
      lambda,
      Y.lambda[, T0 + 1],
      N0 * Re(zeta.lambda^2)
    )
  }
  err.lambda <- Y.lambda %*% c(lambda, -1)

  Y.omega <- if (omega.intercept) {
    sweep(t(Y[, 1:T0, drop = FALSE]), 2, colMeans(t(Y[, 1:T0, drop = FALSE])))
  } else {
    t(Y[, 1:T0, drop = FALSE])
  }
  if (update.omega) {
    omega <- fw_step_cpp(Y.omega[, 1:N0, drop = FALSE], omega, Y.omega[, N0 + 1], T0 * Re(zeta.omega^2))
  }
  err.omega <- Y.omega %*% c(omega, -1)

  val <- Re(zeta.omega^2) * sum(omega^2) + Re(zeta.lambda^2) * sum(lambda^2) + sum(err.omega^2) / T0 + sum(err.lambda^2) / N0
  list(
    val = val,
    lambda = lambda,
    omega = omega,
    err.lambda = err.lambda,
    err.omega = err.omega
  )
}
