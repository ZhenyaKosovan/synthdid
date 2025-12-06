contract3 <- function(X, v) {
  stopifnot(length(dim(X)) == 3, dim(X)[3] == length(v))
  if (length(v) == 0) {
    return(array(0, dim = dim(X)[1:2]))
  }
  contract3_cpp(X, v)
}

# a Frank-Wolfe solver for synthetic control weights using exact line search
sc.weight.fw <- function(Y,
                         zeta,
                         intercept = TRUE,
                         lambda = NULL,
                         min.decrease = 1e-3,
                         max.iter = 1000) {
  T0 <- ncol(Y) - 1
  if (is.null(lambda)) {
    lambda <- rep(1 / T0, T0)
  }
  sc_weight_fw_cpp(Y, zeta, intercept, lambda, min.decrease, max.iter)
}

# A Frank-Wolfe + Gradient solver for lambda, omega, and beta when there are covariates
# Uses the exact line search Frank-Wolfe steps for lambda, omega and (1/t)*gradient steps for beta
# pass update.lambda=FALSE/update.omega=FALSE to fix those weights at initial values, defaulting to uniform 1/T0 and 1/N0
sc.weight.fw.covariates <- function(Y,
                                    X = array(0, dim = c(dim(Y), 0)),
                                    zeta.lambda = 0,
                                    zeta.omega = 0,
                                    lambda.intercept = TRUE,
                                    omega.intercept = TRUE,
                                    min.decrease = 1e-3,
                                    max.iter = 1000,
                                    lambda = NULL,
                                    omega = NULL,
                                    beta = NULL,
                                    update.lambda = TRUE,
                                    update.omega = TRUE) {
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
  }
  list(lambda = weights$lambda, omega = weights$omega, beta = beta, vals = vals)
}


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
    sweep(Y[1:N0, ], 2, colMeans(Y[1:N0, ]))
  } else {
    Y[1:N0, ]
  }
  if (update.lambda) {
    lambda <- fw_step_cpp(
      Y.lambda[, 1:T0],
      lambda,
      Y.lambda[, T0 + 1],
      N0 * Re(zeta.lambda^2)
    )
  }
  err.lambda <- Y.lambda %*% c(lambda, -1)

  Y.omega <- if (omega.intercept) {
    sweep(t(Y[, 1:T0]), 2, colMeans(t(Y[, 1:T0])))
  } else {
    t(Y[, 1:T0])
  }
  if (update.omega) {
    omega <- fw_step_cpp(Y.omega[, 1:N0], omega, Y.omega[, N0 + 1], T0 * Re(zeta.omega^2))
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
