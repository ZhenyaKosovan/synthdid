#' Synthetic Difference-in-Differences Estimation with Formula Interface
#'
#' This function provides a formula-based interface to synthdid estimation,
#' similar to \code{lm()}, \code{plm()}, and \code{glm()}. It automatically
#' handles panel data conversion and returns a rich model object compatible
#' with standard R modeling tools.
#'
#' @importFrom stats coef fitted pnorm qnorm setNames terms
#' @importFrom utils head
#'
#' @param formula A formula of the form \code{outcome ~ treatment} or
#'   \code{outcome ~ treatment | covariates}. The treatment variable should
#'   be a binary indicator (0/1).
#' @param data A data.frame containing the panel data in long format.
#' @param index A character vector of length 2 specifying the names of the
#'   unit and time variables, e.g., \code{c("state", "year")}. If NULL, assumes
#'   the first two columns are unit and time.
#' @param method Estimation method: "synthdid" (default), "sc" (synthetic control),
#'   or "did" (difference-in-differences).
#' @param adoption Treatment adoption pattern: \code{"auto"} (default) detects
#'   whether adoption is simultaneous or staggered, \code{"simultaneous"} forces
#'   the classic single-cohort estimator, \code{"staggered"} forces the
#'   staggered-adoption estimator from Porreca (2022).
#' @param se Logical. If TRUE, compute standard errors. Default is FALSE.
#' @param se_method Standard error method: "bootstrap" (default), "jackknife", or "placebo".
#' @param se_replications Number of replications for bootstrap/placebo standard errors.
#' @param ... Additional arguments passed to \code{synthdid_estimate()} (for
#'   simultaneous adoption) or \code{synthdid_staggered_estimate()} (for
#'   staggered adoption).
#'
#' @return For simultaneous adoption, an object of class
#'   \code{c("synthdid", "synthdid_estimate")}. For staggered adoption, an
#'   object of class \code{"synthdid_staggered"}. Both support standard methods
#'   (\code{print}, \code{summary}, \code{coef}, \code{confint}, \code{plot}).
#'
#'   Simultaneous adoption objects have components:
#'   \item{coefficients}{Treatment effect estimate}
#'   \item{call}{The matched call}
#'   \item{formula}{The formula used}
#'   \item{terms}{The terms object from the formula}
#'   \item{model}{The model frame (if requested)}
#'   \item{Y}{The outcome matrix}
#'   \item{N0}{Number of control units}
#'   \item{T0}{Number of pre-treatment periods}
#'   \item{weights}{List with lambda, omega, and beta weights}
#'   \item{setup}{List describing the problem}
#'   \item{estimator}{Name of estimator used}
#'   \item{index}{Panel index variables}
#'   \item{data_info}{Information about the original data}
#'
#' @examples
#' \donttest{
#' data(california_prop99)
#'
#' # Basic usage
#' result <- synthdid(PacksPerCapita ~ treated,
#'   data = california_prop99,
#'   index = c("State", "Year")
#' )
#'
#' # With standard errors
#' result <- synthdid(PacksPerCapita ~ treated,
#'   data = california_prop99,
#'   index = c("State", "Year"),
#'   se = TRUE,
#'   se_method = "bootstrap"
#' )
#'
#' # Standard R methods work
#' print(result)
#' summary(result)
#' coef(result)
#' confint(result)
#' plot(result)
#'
#' # Compare methods
#' did_result <- synthdid(PacksPerCapita ~ treated,
#'   data = california_prop99,
#'   index = c("State", "Year"),
#'   method = "did"
#' )
#' }
#'
#' @seealso \code{\link{synthdid_estimate}}, \code{\link{panel.matrices}}
#' @export
synthdid <- function(formula,
                     data,
                     index = NULL,
                     method = c("synthdid", "sc", "did"),
                     adoption = c("auto", "simultaneous", "staggered"),
                     se = FALSE,
                     se_method = c("bootstrap", "jackknife", "placebo"),
                     se_replications = SYNTHDID_SE_REPLICATIONS_DEFAULT,
                     ...) {
  # Capture the call
  cl <- match.call()

  # Match arguments
  method <- match.arg(method)
  adoption <- match.arg(adoption)
  se_method <- match.arg(se_method)

  # Validate inputs
  if (!inherits(data, "data.frame")) {
    stop("'data' must be a data.frame")
  }
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula object")
  }

  # Parse formula
  formula_parts <- parse_synthdid_formula(formula)
  outcome_var <- formula_parts$outcome
  treatment_var <- formula_parts$treatment
  covariate_vars <- formula_parts$covariates

  # Determine index variables
  if (is.null(index)) {
    if (ncol(data) < 2) {
      stop("'index' must be specified or data must have at least 2 columns for unit and time")
    }
    index <- colnames(data)[1:2]
    message("Using first two columns as panel index: ", paste(index, collapse = ", "))
  }
  if (length(index) != 2) {
    stop("'index' must be a character vector of length 2: c('unit', 'time')")
  }

  unit_var <- index[1]
  time_var <- index[2]

  # Validate that all required variables exist
  required_vars <- c(unit_var, time_var, outcome_var, treatment_var, covariate_vars)
  missing_vars <- setdiff(required_vars, colnames(data))
  if (length(missing_vars) > 0) {
    stop("Variables not found in data: ", paste(missing_vars, collapse = ", "))
  }

  # Build panel matrices (Y and W) from long-format data
  panel_mats <- build_panel_matrices(data, unit_var, time_var, outcome_var,
                                     treatment_var)
  Y <- panel_mats$Y
  W <- panel_mats$W

  # Detect adoption pattern
  is_staggered <- detect_staggered(W)
  if (adoption == "auto") {
    adoption <- if (is_staggered) "staggered" else "simultaneous"
  } else if (adoption == "simultaneous" && is_staggered) {
    stop("adoption = 'simultaneous' but treatment adoption is staggered. ",
         "Use adoption = 'staggered' or adoption = 'auto'.")
  }

  # Handle covariates if present
  X <- NULL
  if (length(covariate_vars) > 0) {
    X_list <- lapply(covariate_vars, function(cov) {
      cov_mats <- build_panel_matrices(data, unit_var, time_var, cov,
                                       treatment_var)
      cov_mats$Y
    })
    X <- array(unlist(X_list), dim = c(dim(Y), length(covariate_vars)))
  }

  # Route to appropriate estimator

  if (adoption == "staggered") {
    # Staggered adoption: use Porreca-style decomposition
    X_arg <- if (!is.null(X)) X else NULL
    estimate <- synthdid_staggered_estimate(
      Y, W, X = X_arg, method = method, ...
    )

    # Attach formula interface metadata
    attr(estimate, "call") <- cl
    attr(estimate, "formula") <- formula
    attr(estimate, "terms") <- terms(formula)
    attr(estimate, "index") <- setNames(index, c("unit", "time"))
    attr(estimate, "data_info") <- list(
      outcome = outcome_var,
      treatment = treatment_var,
      covariates = covariate_vars,
      n_units = nrow(Y),
      n_periods = ncol(Y),
      adoption = "staggered"
    )

    # Handle SE for staggered
    if (se) {
      se_val <- vcov.synthdid_staggered(estimate, method = se_method,
                                         replications = se_replications)
      attr(estimate, "se") <- sqrt(c(se_val))
      attr(estimate, "se_method") <- se_method
    }

    return(estimate)
  }

  # Simultaneous adoption: existing path
  # Determine N0 and T0 from the simultaneous panel
  W_binary <- W != 0
  w <- rowSums(W_binary) > 0
  T0 <- unname(which(colSums(W_binary) > 0)[1] - 1)
  N0 <- sum(!w)

  # Sort: controls first
  unit_order <- order(W[, T0 + 1], rownames(Y))
  Y <- Y[unit_order, ]
  W <- W[unit_order, ]

  X_sim <- array(dim = c(dim(Y), 0))
  if (!is.null(X)) {
    X_sim <- X[unit_order, , , drop = FALSE]
  }

  # Call appropriate estimator (suppress deprecation warning for internal use)
  estimate <- withCallingHandlers(
    switch(method,
      "synthdid" = synthdid_estimate(Y, N0, T0,
        X = X_sim,
        estimate_se = se,
        se_method = se_method,
        se_replications = se_replications,
        ...
      ),
      "sc" = sc_estimate(Y, N0, T0,
        X = X_sim,
        estimate_se = se,
        se_method = se_method,
        se_replications = se_replications,
        ...
      ),
      "did" = did_estimate(Y, N0, T0,
        X = X_sim,
        estimate_se = se,
        se_method = se_method,
        se_replications = se_replications,
        ...
      )
    ),
    lifecycle_warning_deprecated = function(cnd) invokeRestart("muffleWarning")
  )

  # Enrich the object with formula interface attributes
  attr(estimate, "call") <- cl
  attr(estimate, "formula") <- formula
  attr(estimate, "terms") <- terms(formula)
  attr(estimate, "index") <- setNames(index, c("unit", "time"))
  attr(estimate, "data_info") <- list(
    outcome = outcome_var,
    treatment = treatment_var,
    covariates = covariate_vars,
    n_units = nrow(Y),
    n_periods = ncol(Y),
    adoption = "simultaneous"
  )

  # Update class
  class(estimate) <- c("synthdid", class(estimate))
  return(estimate)
}


#' Build panel matrices from long-format data
#'
#' Converts long-format panel data into Y (outcome) and W (treatment) matrices
#' without requiring simultaneous adoption. This is the staggered-friendly
#' replacement for the internal use of \code{panel.matrices()}.
#'
#' @param data A data.frame in long panel format.
#' @param unit Name of the unit column.
#' @param time Name of the time column.
#' @param outcome Name of the outcome column.
#' @param treatment Name of the treatment column.
#' @return A list with \code{Y} (outcome matrix) and \code{W} (treatment matrix),
#'   both with units as rows and time as columns.
#' @keywords internal
build_panel_matrices <- function(data, unit, time, outcome, treatment) {
  panel <- as.data.frame(data[, c(unit, time, outcome, treatment), drop = FALSE])

  if (anyNA(panel)) {
    stop("Missing values in panel data.")
  }
  if (!all(panel[[treatment]] %in% c(0, 1))) {
    stop("Treatment variable must be binary (0/1).")
  }

  # Convert factors/dates to character for table operations
  panel <- data.frame(
    lapply(panel, function(col) {
      if (is.factor(col) || inherits(col, "Date")) as.character(col) else col
    }),
    stringsAsFactors = FALSE
  )

  # Check balanced panel
  val <- as.vector(table(panel[[unit]], panel[[time]]))
  if (!all(val == 1)) {
    stop("Panel must be balanced: every unit must be observed at every time period.")
  }

  panel <- panel[order(panel[[unit]], panel[[time]]), ]
  units <- unique(panel[[unit]])
  times <- unique(panel[[time]])
  num_units <- length(units)
  num_times <- length(times)

  Y <- matrix(panel[[outcome]], num_units, num_times, byrow = TRUE,
              dimnames = list(units, times))
  W <- matrix(panel[[treatment]], num_units, num_times, byrow = TRUE,
              dimnames = list(units, times))

  list(Y = Y, W = W)
}


#' Detect whether treatment adoption is staggered
#'
#' Returns \code{TRUE} if units adopt treatment at different times,
#' \code{FALSE} if all treated units adopt simultaneously.
#'
#' @param W Binary treatment matrix (units x time).
#' @return Logical scalar.
#' @keywords internal
detect_staggered <- function(W) {
  first_treat <- compute_first_treat_time(W)
  treated_times <- first_treat[is.finite(first_treat)]
  if (length(treated_times) == 0) {
    stop("No treated units found in the panel.")
  }
  length(unique(treated_times)) > 1
}


#' Parse synthdid formula
#'
#' Internal function to parse formula of the form outcome ~ treatment | covariates
#'
#' @param formula A formula object
#' @return A list with outcome, treatment, and covariates
#' @keywords internal
parse_synthdid_formula <- function(formula) {
  if (length(formula) != 3) {
    stop("Formula must have both left and right hand sides")
  }

  # Extract outcome (left-hand side)
  outcome <- as.character(formula[[2]])

  # Extract right-hand side
  rhs <- formula[[3]]

  # Check for covariate separator |
  if (length(rhs) > 1 && as.character(rhs[[1]]) == "|") {
    # Formula has covariates: outcome ~ treatment | cov1 + cov2
    treatment_part <- rhs[[2]]
    covariate_part <- rhs[[3]]

    treatment <- as.character(treatment_part)
    covariates <- extract_vars_from_expr(covariate_part)
  } else {
    # Simple formula: outcome ~ treatment
    treatment <- as.character(rhs)
    covariates <- character(0)
  }

  # Validate: treatment must be a single variable name
  if (length(treatment) != 1) {
    stop("Unsupported formula syntax. Only 'outcome ~ treatment' or ",
         "'outcome ~ treatment | covariates' is supported.")
  }

  list(
    outcome = outcome,
    treatment = treatment,
    covariates = covariates
  )
}


#' Extract variable names from formula expression
#'
#' @param expr An R expression from formula parsing
#' @return Character vector of variable names
#' @keywords internal
extract_vars_from_expr <- function(expr) {
  if (is.name(expr)) {
    return(as.character(expr))
  } else if (is.call(expr)) {
    if (as.character(expr[[1]]) == "+") {
      # Handle addition: recurse on both sides
      return(c(
        extract_vars_from_expr(expr[[2]]),
        extract_vars_from_expr(expr[[3]])
      ))
    } else {
      # For other functions like I(), log(), etc., extract the variable
      # This is a simplified version - could be enhanced
      vars <- unlist(lapply(as.list(expr)[-1], extract_vars_from_expr))
      return(vars[vars != ""])
    }
  }
  character(0)
}


#' Print method for synthdid objects
#' @param x A synthdid object
#' @param ... Additional arguments (currently ignored)
#' @return Invisibly returns the original object
#' @examples
#' \donttest{
#' data(california_prop99)
#' result <- synthdid(PacksPerCapita ~ treated,
#'   data = california_prop99,
#'   index = c("State", "Year")
#' )
#' print(result)
#' }
#' @export
print.synthdid <- function(x, ...) {
  cat("Synthetic Difference-in-Differences Estimate\n\n")

  # Show call
  cat("Call:\n")
  print(attr(x, "call"))
  cat("\n")

  # Show treatment effect
  tau <- coef(x)
  cat("Treatment Effect: ", format(tau, digits = 4), "\n")

  # Show standard error if available
  se_attr <- attr(x, "se")
  if (!is.null(se_attr) && !is.na(se_attr)) {
    cat("Standard Error:   ", format(se_attr, digits = 4),
      " (", attr(x, "se_method"), ")\n",
      sep = ""
    )
  }

  # Show dimensions
  setup <- attr(x, "setup")
  N0 <- setup$N0
  N1 <- nrow(setup$Y) - N0
  T0 <- setup$T0
  T1 <- ncol(setup$Y) - T0

  cat("\n")
  cat("Units:        ", N0, " control, ", N1, " treated\n", sep = "")
  cat("Time Periods: ", T0, " pre-treatment, ", T1, " post-treatment\n", sep = "")

  # Show convergence status
  conv <- attr(x, "convergence")
  if (!is.null(conv)) {
    cat("\nConvergence:  ")
    if (conv$overall_converged) {
      cat("CONVERGED\n")
    } else {
      cat("NOT CONVERGED - ")
      # Show which components failed
      failed_parts <- character()
      if (!is.null(conv$lambda) && !conv$lambda$converged) {
        failed_parts <- c(failed_parts, sprintf(
          "lambda (%d/%d iters)",
          conv$lambda$iterations,
          conv$lambda$max_iter
        ))
      }
      if (!is.null(conv$omega) && !conv$omega$converged) {
        failed_parts <- c(failed_parts, sprintf(
          "omega (%d/%d iters)",
          conv$omega$iterations,
          conv$omega$max_iter
        ))
      }
      if (!is.null(conv$joint) && !conv$joint$converged) {
        failed_parts <- c(failed_parts, sprintf(
          "joint (%d/%d iters)",
          conv$joint$iterations,
          conv$joint$max_iter
        ))
      }
      cat(paste(failed_parts, collapse = ", "), "\n")
      cat("              Consider increasing max.iter. Use summary() for details.\n")
    }
  }

  invisible(x)
}


#' Extract coefficients from synthdid object
#' @param object A synthdid object
#' @param ... Additional arguments (currently ignored)
#' @return Named numeric vector with the treatment effect estimate
#' @examples
#' \donttest{
#' data(california_prop99)
#' result <- synthdid(PacksPerCapita ~ treated,
#'   data = california_prop99,
#'   index = c("State", "Year")
#' )
#' coef(result)
#' }
#' @export
coef.synthdid <- function(object, ...) {
  setNames(c(object), "treated")
}


#' Compute confidence intervals for synthdid object
#' @param object A synthdid object
#' @param parm Ignored (included for S3 generic compatibility)
#' @param level Confidence level (default: 0.95)
#' @param ... Additional arguments (currently ignored)
#' @return A matrix with lower and upper confidence bounds
#' @examples
#' \donttest{
#' data(california_prop99)
#' result <- synthdid(PacksPerCapita ~ treated,
#'   data = california_prop99,
#'   index = c("State", "Year"),
#'   se = TRUE,
#'   se_method = "jackknife"
#' )
#' confint(result)
#' confint(result, level = 0.90)
#' }
#' @export
confint.synthdid <- function(object, parm, level = 0.95, ...) {
  tau <- coef(object)
  se <- attr(object, "se")
  if (is.null(se) || is.na(se)) {
    warning("Standard error not available; cannot compute confidence interval")
    return(matrix(c(NA, NA),
      nrow = 1, ncol = 2,
      dimnames = list("att", c("Lower", "Upper"))
    ))
  }

  z <- qnorm((1 + level) / 2)

  ci <- cbind(
    Lower = tau - z * se,
    Upper = tau + z * se
  )
  rownames(ci) <- "att"

  ci
}


#' Extract fitted values from synthdid object
#' @param object A synthdid object
#' @param ... Additional arguments (currently ignored)
#' @return A matrix of fitted values for the control units
#' @examples
#' \donttest{
#' data(california_prop99)
#' result <- synthdid(PacksPerCapita ~ treated,
#'   data = california_prop99,
#'   index = c("State", "Year")
#' )
#' fitted_vals <- fitted(result)
#' head(fitted_vals)
#' }
#' @export
fitted.synthdid <- function(object, ...) {
  setup <- attr(object, "setup")
  weights <- attr(object, "weights")

  # Compute fitted values for control units
  Y <- setup$Y
  X.beta <- contract3(setup$X, weights$beta)

  # Create synthetic control predictions
  N0 <- setup$N0
  N1 <- nrow(Y) - N0
  T0 <- setup$T0
  T1 <- ncol(Y) - T0

  # Fitted values for treated units based on synthetic control
  fitted_matrix <- matrix(NA, nrow = nrow(Y), ncol = ncol(Y))

  # For control units, fitted = actual (they define the synthetic control)
  fitted_matrix[1:N0, ] <- Y[1:N0, ] - X.beta[1:N0, ]

  # For treated units, use synthetic control
  synthetic_control <- t(weights$omega) %*% (Y[1:N0, ] - X.beta[1:N0, ])
  for (i in (N0 + 1):nrow(Y)) {
    fitted_matrix[i, ] <- synthetic_control + X.beta[i, ]
  }

  fitted_matrix
}


#' Extract residuals from synthdid object
#' @param object A synthdid object
#' @param type Type of residuals: "control" (default), "pretreatment", or "all"
#' @param ... Additional arguments (currently ignored)
#' @return A matrix of residuals
#' @examples
#' \donttest{
#' data(california_prop99)
#' result <- synthdid(PacksPerCapita ~ treated,
#'   data = california_prop99,
#'   index = c("State", "Year")
#' )
#' resid_control <- residuals(result, type = "control")
#' resid_pretreat <- residuals(result, type = "pretreatment")
#' head(resid_control)
#' }
#' @export
residuals.synthdid <- function(object, type = c("control", "pretreatment", "all"), ...) {
  type <- match.arg(type)
  setup <- attr(object, "setup")
  weights <- attr(object, "weights")

  Y <- setup$Y
  X.beta <- contract3(setup$X, weights$beta)
  N0 <- setup$N0
  T0 <- setup$T0

  if (type == "control") {
    # Residuals for control units
    Y_control <- (Y - X.beta)[1:N0, , drop = FALSE]
    synthetic <- matrix(t(weights$omega) %*% Y_control,
      nrow = N0,
      ncol = ncol(Y),
      byrow = TRUE
    )
    omega_resid <- Y_control - synthetic
    return(omega_resid)
  } else if (type == "pretreatment") {
    # Residuals for pre-treatment periods
    Y_pre <- (Y - X.beta)[, 1:T0, drop = FALSE]
    synthetic <- matrix(Y_pre %*% weights$lambda,
      nrow = nrow(Y),
      ncol = T0,
      byrow = FALSE
    )
    lambda_resid <- Y_pre - synthetic
    return(lambda_resid)
  } else {
    # All residuals
    fitted_vals <- fitted(object)
    return(Y - fitted_vals)
  }
}


#' Predictions from synthdid object
#' @param object A synthdid object
#' @param newdata Currently not supported (reserved for future use)
#' @param type Type of prediction: "counterfactual" (default), "treated", or "effect"
#' @param ... Additional arguments (currently ignored)
#' @return A matrix or vector of predictions
#' @examples
#' \donttest{
#' data(california_prop99)
#' result <- synthdid(PacksPerCapita ~ treated,
#'   data = california_prop99,
#'   index = c("State", "Year")
#' )
#' # Counterfactual (what would have happened without treatment)
#' pred_counterfactual <- predict(result, type = "counterfactual")
#' # Observed treated outcomes
#' pred_treated <- predict(result, type = "treated")
#' # Treatment effect over time
#' pred_effect <- predict(result, type = "effect")
#' }
#' @export
predict.synthdid <- function(object,
                             newdata = NULL,
                             type = c("counterfactual", "treated", "effect"), ...) {
  type <- match.arg(type)
  setup <- attr(object, "setup")
  weights <- attr(object, "weights")

  Y <- setup$Y
  X.beta <- contract3(setup$X, weights$beta)
  N0 <- setup$N0
  N1 <- nrow(Y) - N0
  T0 <- setup$T0
  T1 <- ncol(Y) - T0

  if (type == "counterfactual") {
    # What would have happened to treated units without treatment?
    synthetic_control <- as.numeric(t(weights$omega) %*% (Y[1:N0, ] - X.beta[1:N0, ]))
    counterfactual <- matrix(rep(synthetic_control, each = N1),
      nrow = N1,
      ncol = ncol(Y)
    ) + X.beta[(N0 + 1):nrow(Y), , drop = FALSE]
    rownames(counterfactual) <- rownames(Y)[(N0 + 1):nrow(Y)]
    colnames(counterfactual) <- colnames(Y)
    return(counterfactual)
  } else if (type == "treated") {
    # Actual treated outcomes
    treated_outcomes <- Y[(N0 + 1):nrow(Y), , drop = FALSE]
    return(treated_outcomes)
  } else {
    # Treatment effect by period (effect curve)
    return(withCallingHandlers(
      synthdid_effect_curve(object),
      lifecycle_warning_deprecated = function(cnd) invokeRestart("muffleWarning")
    ))
  }
}


#' Update a synthdid model
#'
#' Update and re-fit a synthdid model with different options
#'
#' @param object A synthdid object
#' @param formula. An optional new formula
#' @param method An optional new method
#' @param ... Additional arguments
#' @return An updated synthdid object
#' @examples
#' \donttest{
#' data(california_prop99)
#' result <- synthdid(PacksPerCapita ~ treated,
#'   data = california_prop99,
#'   index = c("State", "Year")
#' )
#' # Update to use different method
#' result_did <- update(result, method = "did")
#' c(sdid = coef(result), did = coef(result_did))
#' }
#' @export
update.synthdid <- function(object, formula. = NULL, method = NULL, ...) {
  call <- attr(object, "call")

  if (!is.null(formula.)) {
    call$formula <- formula.
  }
  if (!is.null(method)) {
    call$method <- method
  }

  # Update additional arguments
  extras <- list(...)
  if (length(extras) > 0) {
    for (nm in names(extras)) {
      call[[nm]] <- extras[[nm]]
    }
  }

  eval(call, parent.frame())
}


#' Model frame for synthdid
#'
#' Extract model frame from synthdid object
#' @param formula A synthdid object
#' @param ... Additional arguments (currently ignored)
#' @return A list with call and terms information
#' @examples
#' \donttest{
#' data(california_prop99)
#' result <- synthdid(PacksPerCapita ~ treated,
#'   data = california_prop99,
#'   index = c("State", "Year")
#' )
#' model.frame(result)
#' }
#' @export
model.frame.synthdid <- function(formula, ...) {
  # For now, return NULL or minimal info
  # Could enhance to return actual model frame
  list(
    call = attr(formula, "call"),
    terms = attr(formula, "terms")
  )
}
