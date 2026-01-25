#' Collapse an outcome matrix to synthetic control form
#'
#' Aggregates post-treatment rows and columns into averages so the final row and
#' column represent treated units and treated periods. This is a convenience for
#' algorithms that expect a \code{(N0 + 1) x (T0 + 1)} matrix with the treated
#' block compressed to a single row and column.
#'
#' @param Y Outcome matrix with control units first and pre-treatment periods first.
#' @param N0 Number of control units (rows).
#' @param T0 Number of pre-treatment periods (columns).
#'
#' @return A numeric matrix of dimension \code{N0 + 1} by \code{T0 + 1} where
#'   the final row and column contain treated averages.
#' @keywords internal
collapsed.form <- function(Y, N0, T0) {
  N <- nrow(Y)
  T <- ncol(Y)
  rbind(
    cbind(Y[1:N0, 1:T0, drop = FALSE], rowMeans(Y[1:N0, (T0 + 1):T, drop = FALSE])),
    cbind(t(colMeans(Y[(N0 + 1):N, 1:T0, drop = FALSE])), mean(Y[(N0 + 1):N, (T0 + 1):T, drop = FALSE]))
  )
}

#' Combine decreasing sequences while respecting terminal NA markers
#'
#' Treats \code{NA} as a sentinel indicating that a series has stopped
#' decreasing and should hold its last non-\code{NA} value. Returns the
#' element-wise sum with terminal \code{NA} preserved where both series ended.
#'
#' @param x,y Numeric vectors, typically monotone decreasing with trailing
#'   \code{NA} entries.
#'
#' @return Numeric vector of pairwise sums with shared terminal \code{NA}
#'   retained.
#' @keywords internal
pairwise.sum.decreasing <- function(x, y) {
  na.x <- is.na(x)
  na.y <- is.na(y)
  x[is.na(x)] <- min(x[!na.x])
  y[is.na(y)] <- min(y[!na.y])
  pairwise.sum <- x + y
  pairwise.sum[na.x & na.y] <- NA
  return(pairwise.sum)
}

#' Convert a long (balanced) panel to a wide matrix
#'
#' Converts a data set in panel form to matrix format required by synthdid estimators.
#' A typical long panel date set looks like \[unit, time, outcome, treatment\]. Synthdid
#' requires a balanced panel with simultaneous adoption of treatment: each unit must be observed
#' at all times, and all treated units must begin treatment simultaneosly. This function
#' creates num.units x num.time.periods matrices Y and W of outcomes and treatment indicators.
#' In these matrices, columns are sorted by time, and by default (when treated.last=TRUE),
#' rows for control units appear before those of treated units.
#'
#' @param panel A data.frame with columns consisting of units, time, outcome, and treatment indicator.
#' @param unit The column number/name corresponding to the unit identifier. Default is 1.
#' @param time The column number/name corresponding to the time identifier. Default is 2.
#' @param outcome The column number/name corresponding to the outcome identifier. Default is 3.
#' @param treatment The column number/name corresponding to the treatment status. Default is 4.
#' @param treated.last Should we sort the rows of Y and W so treated units are last. If FALSE, sort by unit number/name. Default is TRUE.
#' @return A list with entries `Y`: the data matrix, `N0`: the number of control units, `T0`:
#'  the number of time periods before treatment, `W`: the matrix of treatment indicators.
#'
#' @examples
#' \donttest{
#' # Load tobacco sales in long panel format.
#' data("california_prop99")
#' # Transform to N*T matrix format required for synthdid,
#' # where N is the number of units and T the time periods.
#' setup <- panel.matrices(california_prop99, unit = 1, time = 2, outcome = 3, treatment = 4)
#'
#' # Compute synthdid estimate
#' synthdid_estimate(setup$Y, setup$N0, setup$T0)
#' }
#'
#' @export
panel.matrices <- function(panel, unit = 1, time = 2, outcome = 3, treatment = 4, treated.last = TRUE) {
  # TODO: add support for covariates X, i.e. could keep all other columns
  keep <- c(unit, time, outcome, treatment)
  if (!all(keep %in% 1:ncol(panel) | keep %in% colnames(panel))) {
    stop("Column identifiers should be either integer or column names in `panel`.")
  }
  index.to.name <- function(x) {
    if (x %in% 1:ncol(panel)) {
      colnames(panel)[x]
    } else {
      x
    }
  }
  unit <- index.to.name(unit)
  time <- index.to.name(time)
  outcome <- index.to.name(outcome)
  treatment <- index.to.name(treatment)
  keep <- c(unit, time, outcome, treatment)

  panel <- panel[keep]
  if (!is.data.frame(panel)) {
    stop("Unsupported input type `panel.`")
  }
  if (anyNA(panel)) {
    stop("Missing values in `panel`.")
  }
  if (length(unique(panel[, treatment])) == 1) {
    stop("There is no variation in treatment status.")
  }
  if (!all(panel[, treatment] %in% c(0, 1))) {
    stop("The treatment status should be in 0 or 1.")
  }
  # Convert potential factor/date columns to character
  panel <- data.frame(
    lapply(panel, function(col) {
      if (is.factor(col) || inherits(col, "Date")) as.character(col) else col
    }),
    stringsAsFactors = FALSE
  )
  val <- as.vector(table(panel[, unit], panel[, time]))
  if (!all(val == 1)) {
    stop("Input `panel` must be a balanced panel: it must have an observation for every unit at every time.")
  }

  panel <- panel[order(panel[, unit], panel[, time]), ]
  num.years <- length(unique(panel[, time]))
  num.units <- length(unique(panel[, unit]))
  Y <- matrix(panel[, outcome], num.units, num.years,
    byrow = TRUE,
    dimnames = list(unique(panel[, unit]), unique(panel[, time]))
  )
  W <- matrix(panel[, treatment], num.units, num.years,
    byrow = TRUE,
    dimnames = list(unique(panel[, unit]), unique(panel[, time]))
  )
  W_binary <- W != 0
  w <- rowSums(W_binary) > 0 # indicator for units that are treated at any time
  T0 <- unname(which(colSums(W_binary) > 0)[1] - 1) # last period nobody is treated
  N0 <- sum(!w)

  if (T0 < 1) {
    stop("Treatment starts in the first period; at least one pre-treatment period is required.")
  }

  if (!(all(W_binary[!w, ] == 0) && all(W_binary[, 1:T0] == 0) && all(W_binary[w, (T0 + 1):ncol(Y)] == 1))) {
    stop("The package cannot use this data. Treatment adoption is not simultaneous.")
  }

  unit.order <- if (treated.last) {
    order(W[, T0 + 1], rownames(Y))
  } else {
    1:nrow(Y)
  }
  return(list(Y = Y[unit.order, ], N0 = N0, T0 = T0, W = W[unit.order, ]))
}

#' Get timesteps from panel matrix Y
#'
#' timesteps are stored as colnames(Y), but column names cannot be Date objects.
#' Instead, we use strings. If they are strings convertible to dates, return that
#'
#' @param Y a matrix
#' @return its column names interpreted as Dates if possible
#' @examples
#' \donttest{
#' data(california_prop99)
#' setup <- panel.matrices(california_prop99)
#' timesteps(setup$Y)
#' }
#' @export
timesteps <- function(Y) {
  labels <- colnames(Y)
  if (is.null(labels)) {
    return(labels)
  }
  parsed <- tryCatch(
    suppressWarnings(as.Date(labels)),
    error = function(e) rep(NA, length(labels))
  )
  if (length(parsed) == 0 || any(is.na(parsed))) {
    return(labels)
  }
  return(parsed)
}


## define some convenient accessors
setOldClass("synthdid_estimate")
setOldClass("synthdid")

#' Create a slot accessor function
#'
#' Builds a small closure used to define S4 generics that delegate to list-style
#' elements stored on synthdid objects.
#'
#' @param name Name of the element to retrieve.
#'
#' @return A function that extracts the named element from its input.
#' @keywords internal
get_slot <- function(name) {
  function(object) {
    object[[name]]
  }
}

setGeneric("weights")
setGeneric("Y", get_slot("Y"))
setGeneric("lambda", get_slot("lambda"))
setGeneric("omega", get_slot("omega"))
setMethod(weights,
  signature = "synthdid_estimate",
  definition = function(object) {
    attr(object, "weights")
  }
)
setMethod(weights,
  signature = "synthdid",
  definition = function(object) {
    attr(object, "weights")
  }
)
setMethod(Y,
  signature = "synthdid_estimate",
  definition = function(object) {
    attr(object, "setup")$Y
  }
)
setMethod(Y,
  signature = "synthdid",
  definition = function(object) {
    attr(object, "setup")$Y
  }
)
setMethod(lambda,
  signature = "synthdid_estimate",
  definition = function(object) {
    lambda(weights(object))
  }
)
setMethod(lambda,
  signature = "synthdid",
  definition = function(object) {
    lambda(weights(object))
  }
)
setMethod(omega,
  signature = "synthdid_estimate",
  definition = function(object) {
    omega(weights(object))
  }
)
setMethod(omega,
  signature = "synthdid",
  definition = function(object) {
    omega(weights(object))
  }
)


#' Generate a synthetic low-rank panel for testing
#'
#' Creates a reproducible, low-rank outcome matrix with a block treatment
#' indicator. Useful for unit tests or demonstrations that need structured data
#' with controlled noise and treatment effects.
#'
#' @return A list containing the outcome matrix \code{Y}, latent matrix
#'   \code{L}, and integers \code{N0} and \code{T0} describing the number of
#'   control units and pre-treatment periods.
#' @examples
#' \donttest{
#' # Generate synthetic data for testing
#' data <- random.low.rank()
#'
#' # Run estimation
#' tau.hat <- synthdid_estimate(data$Y, data$N0, data$T0)
#' print(tau.hat)
#'
#' # True treatment effect is 1
#' tau.hat # Should be close to 1
#' }
#' @export
random.low.rank <- function() {
  n_0 <- 100
  n_1 <- 10
  T_0 <- 120
  T_1 <- 20
  n <- n_0 + n_1
  T <- T_0 + T_1
  tau <- 1
  sigma <- 0.5
  rank <- 2
  rho <- 0.7
  var <- outer(1:T, 1:T, FUN = function(x, y) rho^(abs(x - y)))

  W <- (1:n > n_0) %*% t(1:T > T_0)
  U <- matrix(rpois(rank * n, sqrt(sample(1:n)) / sqrt(n)), n, rank)
  V <- matrix(rpois(rank * T, sqrt(1:T) / sqrt(T)), T, rank)
  alpha <- outer(10 * sample(1:n) / n, rep(1, T))
  beta <- outer(rep(1, n), 10 * (1:T) / T)
  mu <- U %*% t(V) + alpha + beta
  error <- mvtnorm::rmvnorm(n, sigma = var, method = "chol")
  Y <- mu + tau * W + sigma * error
  rownames(Y) <- 1:n
  colnames(Y) <- 1:T
  return(list(Y = Y, L = mu, N0 = n_0, T0 = T_0))
}

#' Check convergence status of synthdid estimate
#'
#' Lightweight function to check if the optimization converged properly.
#' Uses information already computed during estimation (zero overhead).
#'
#' @param estimate A synthdid_estimate object
#' @return Logical indicating whether optimization converged
#' @examples
#' \donttest{
#' data(california_prop99)
#' setup <- panel.matrices(california_prop99)
#' tau.hat <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
#'
#' # Quick convergence check
#' synthdid_converged(tau.hat)
#'
#' # If FALSE, get more details
#' if (!synthdid_converged(tau.hat)) {
#'   synthdid_convergence_info(tau.hat)
#' }
#' }
#' @export
synthdid_converged <- function(estimate) {
  conv <- attr(estimate, "convergence")
  if (is.null(conv)) {
    # Old estimate without convergence info
    return(NA)
  }
  conv$overall_converged
}

#' Get convergence diagnostics for synthdid estimate
#'
#' Returns detailed convergence information including iteration counts
#' and which components (lambda, omega, joint) converged.
#'
#' @param estimate A synthdid_estimate object
#' @return List with convergence diagnostics
#' @examples
#' \donttest{
#' data(california_prop99)
#' setup <- panel.matrices(california_prop99)
#' tau.hat <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
#'
#' # Get detailed convergence diagnostics
#' conv_info <- synthdid_convergence_info(tau.hat)
#' print(conv_info)
#'
#' # Access components
#' conv_info$lambda$iterations
#' conv_info$omega$iterations
#' conv_info$overall_converged
#' }
#' @export
synthdid_convergence_info <- function(estimate) {
  conv <- attr(estimate, "convergence")
  if (is.null(conv)) {
    message("This estimate does not have convergence information.")
    message("Re-run with a recent version of synthdid to get convergence diagnostics.")
    return(NULL)
  }

  # Format output for readability
  info <- list(
    overall_converged = conv$overall_converged
  )

  if (!is.null(conv$lambda)) {
    info$lambda <- data.frame(
      converged = conv$lambda$converged,
      iterations = conv$lambda$iterations,
      max_iter = conv$lambda$max_iter,
      utilization = sprintf("%.1f%%", 100 * conv$lambda$iterations / conv$lambda$max_iter)
    )
  }

  if (!is.null(conv$omega)) {
    info$omega <- data.frame(
      converged = conv$omega$converged,
      iterations = conv$omega$iterations,
      max_iter = conv$omega$max_iter,
      utilization = sprintf("%.1f%%", 100 * conv$omega$iterations / conv$omega$max_iter)
    )
  }

  if (!is.null(conv$joint)) {
    info$joint <- data.frame(
      converged = conv$joint$converged,
      iterations = conv$joint$iterations,
      max_iter = conv$joint$max_iter,
      utilization = sprintf("%.1f%%", 100 * conv$joint$iterations / conv$joint$max_iter)
    )
  }

  class(info) <- c("synthdid_convergence", "list")
  info
}

#' Print method for synthdid convergence info
#' @param x A synthdid_convergence object
#' @param ... Additional arguments (ignored)
#' @return Invisibly returns the original object
#' @examples
#' \donttest{
#' data(california_prop99)
#' setup <- panel.matrices(california_prop99)
#' tau.hat <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
#' conv_info <- synthdid_convergence_info(tau.hat)
#' print(conv_info)
#' }
#' @export
print.synthdid_convergence <- function(x, ...) {
  cat("Synthdid Convergence Diagnostics\n")
  cat("=================================\n\n")

  cat("Overall Status:", if (x$overall_converged) "CONVERGED" else "NOT CONVERGED", "\n\n")

  if (!is.null(x$lambda)) {
    cat("Lambda weights optimization:\n")
    print(x$lambda, row.names = FALSE)
    cat("\n")
  }

  if (!is.null(x$omega)) {
    cat("Omega weights optimization:\n")
    print(x$omega, row.names = FALSE)
    cat("\n")
  }

  if (!is.null(x$joint)) {
    cat("Joint optimization (with covariates):\n")
    print(x$joint, row.names = FALSE)
    cat("\n")
  }

  if (!x$overall_converged) {
    cat("Recommendation: Consider increasing max.iter or relaxing min.decrease threshold.\n")
  }

  invisible(x)
}


#' Estimate memory requirements for synthdid computation
#'
#' Provides rough estimates of peak memory usage for different synthdid operations.
#' Useful for assessing feasibility of large-scale computations before running them.
#'
#' @param N Total number of units (control + treated)
#' @param T Total number of time periods (pre + post)
#' @param K Number of covariates (default: 0)
#' @param replications Number of bootstrap/placebo replications for SE (default: 200)
#' @param include_se Logical. Include memory for SE computation (default: TRUE)
#'
#' @return List with memory estimates in GB for different components
#' @export
#'
#' @examples
#' \dontrun{
#' # Estimate memory for California Prop 99 scale problem
#' synthdid_memory_estimate(N = 39, T = 31, K = 0, replications = 1000)
#'
#' # Large problem with covariates
#' synthdid_memory_estimate(N = 1000, T = 100, K = 50, replications = 1000)
#' }
synthdid_memory_estimate <- function(N, T, K = 0,
                                     replications = SYNTHDID_SE_REPLICATIONS_DEFAULT,
                                     include_se = TRUE) {
  if (N <= 0 || is.na(N)) {
    stop("N should be positive scalar.")
  }
  if (T <= 0 || is.na(T)) {
    stop("T should be positive scalar.")
  }
  if (K <= 0 || is.na(K)) {
    stop("K should be positive scalar.")
  }
  if (replications <= 0) {
    stop("replications should be positive scalar.")
  }
  # Use package constants
  bytes_per_double <- SYNTHDID_BYTES_PER_DOUBLE
  gb_divisor <- SYNTHDID_GB_DIVISOR

  # Core data structures
  Y_size <- N * T * bytes_per_double
  X_size <- if (K > 0) N * T * K * bytes_per_double else 0

  # Weight vectors
  weights_size <- (N + T) * bytes_per_double

  # Working memory for optimization (gradient, objective values, etc.)
  # Frank-Wolfe needs: gradient (N or T), residuals (N*T), etc.
  working_memory <- max(N, T) * bytes_per_double + N * T * bytes_per_double

  # Core estimation memory
  core_memory_gb <- (Y_size + X_size + weights_size + working_memory) / gb_divisor

  result <- list(
    N = N,
    T = T,
    K = K,
    core_estimation_gb = core_memory_gb
  )

  if (include_se) {
    # Bootstrap/jackknife: needs to store all estimates
    # With parallel processing, may need multiple copies
    n_parallel <- min(parallel::detectCores(), SYNTHDID_MAX_PARALLEL_CORES)
    se_memory_gb <- (replications * bytes_per_double +
      n_parallel * (Y_size + X_size)) / gb_divisor

    result$se_estimation_gb <- se_memory_gb
    result$total_gb <- core_memory_gb + se_memory_gb
    result$replications <- replications
  } else {
    result$total_gb <- core_memory_gb
  }

  class(result) <- "synthdid_memory_estimate"
  result
}


#' Print method for synthdid memory estimate
#' @param x A synthdid_memory_estimate object
#' @param ... Additional arguments (ignored)
#' @return Invisibly returns the original object
#' @examples
#' \donttest{
#' # Estimate memory for a large problem
#' mem <- synthdid_memory_estimate(N = 1000, T = 200, K = 10, replications = 500)
#' print(mem)
#' }
#' @export
print.synthdid_memory_estimate <- function(x, ...) {
  cat("Synthdid Memory Estimate\n")
  cat("========================\n\n")
  cat(sprintf("Problem size: N=%d units, T=%d periods", x$N, x$T))
  if (x$K > 0) cat(sprintf(", K=%d covariates", x$K))
  cat("\n\n")

  cat(sprintf("Core estimation:     %.2f GB\n", x$core_estimation_gb))

  if (!is.null(x$se_estimation_gb)) {
    cat(sprintf(
      "SE computation:      %.2f GB (%d replications)\n",
      x$se_estimation_gb, x$replications
    ))
    cat(sprintf("Total (peak):        %.2f GB\n", x$total_gb))
  } else {
    cat(sprintf("Total:               %.2f GB\n", x$total_gb))
  }

  cat("\n")
  if (x$total_gb > SYNTHDID_MEMORY_WARNING_HIGH_GB) {
    cat(sprintf(
      "WARNING: Large memory requirement (>%dGB). Consider:\n",
      SYNTHDID_MEMORY_WARNING_HIGH_GB
    ))
    cat("  - Using a high-memory machine\n")
    cat("  - Reducing number of SE replications\n")
    cat("  - Using jackknife SE (faster, less memory than bootstrap)\n")
    cat("  - Setting estimate_se=FALSE and computing SE separately\n")
  } else if (x$total_gb > SYNTHDID_MEMORY_WARNING_MODERATE_GB) {
    cat(sprintf(
      "NOTE: Moderate memory requirement (>%dGB).\n",
      SYNTHDID_MEMORY_WARNING_MODERATE_GB
    ))
    cat("  - Ensure sufficient RAM is available\n")
    cat("  - Consider closing other applications\n")
  } else {
    cat("Memory requirement is reasonable for most systems.\n")
  }

  invisible(x)
}
