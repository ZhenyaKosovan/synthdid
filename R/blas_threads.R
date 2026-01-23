#' BLAS Thread Management for Parallel Processing
#'
#' These functions manage BLAS/LAPACK thread counts to prevent thread
#' oversubscription when using parallel workers with multi-threaded BLAS
#' libraries (OpenBLAS, MKL, Apple Accelerate, etc.).
#'
#' @importFrom future plan
#' @name blas_threads
#' @keywords internal
NULL

#' Detect if future plan is using parallel workers
#'
#' @return Logical indicating if parallel processing is active
#' @keywords internal
is_future_parallel <- function() {
  # Check if future is loaded
  if (!requireNamespace("future", quietly = TRUE)) {
    return(FALSE)
  }

  # Get current plan
  plan_class <- class(future::plan())[1]

  # Parallel plans include: multicore, multisession, cluster
  # Sequential plans: sequential, transparent
  parallel_plans <- c("multicore", "multisession", "cluster")

  plan_class %in% parallel_plans
}

#' Get current BLAS thread count
#'
#' Attempts to detect the number of threads used by the BLAS library.
#' Supports OpenBLAS, MKL, and Apple Accelerate.
#'
#' @return Integer number of BLAS threads, or NULL if cannot be determined
#' @keywords internal
get_blas_threads <- function() {
  # Try OpenBLAS
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    return(RhpcBLASctl::blas_get_num_procs())
  }

  # Try environment variables
  for (var in c("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")) {
    val <- Sys.getenv(var)
    if (val != "") {
      return(as.integer(val))
    }
  }

  # Cannot determine
  return(NULL)
}

#' Set BLAS thread count
#'
#' Sets the number of threads for BLAS operations. This is critical when
#' using parallel workers to avoid thread oversubscription.
#'
#' @param n Number of threads to use (typically 1 when using parallel workers)
#' @return Previous thread count (for restoration), or NULL if unsupported
#' @keywords internal
set_blas_threads <- function(n) {
  # Try OpenBLAS/MKL via RhpcBLASctl (most reliable)
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    old <- RhpcBLASctl::blas_get_num_procs()
    RhpcBLASctl::blas_set_num_threads(n)
    return(old)
  }

  # Try environment variables (less reliable, may require restart)
  old_vals <- list()
  for (var in c("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")) {
    old_val <- Sys.getenv(var)
    if (old_val != "") {
      old_vals[[var]] <- old_val
    }
    Sys.setenv(setNames(as.character(n), var))
  }

  if (length(old_vals) > 0) {
    return(old_vals)
  }

  # Cannot control BLAS threads
  return(NULL)
}

#' Restore BLAS thread count
#'
#' Restores BLAS threads to a previously saved state.
#'
#' @param old Previous state from set_blas_threads()
#' @keywords internal
restore_blas_threads <- function(old) {
  if (is.null(old)) {
    return(invisible(NULL))
  }

  if (is.numeric(old)) {
    # RhpcBLASctl case
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      RhpcBLASctl::blas_set_num_threads(old)
    }
  } else if (is.list(old)) {
    # Environment variable case
    for (var in names(old)) {
      Sys.setenv(setNames(old[[var]], var))
    }
  }

  invisible(NULL)
}

#' Execute code with BLAS thread management
#'
#' Wrapper function that automatically manages BLAS threads when using
#' parallel processing. Sets BLAS to single-threaded mode if parallel
#' workers are detected, then restores original thread count.
#'
#' @param expr Expression to evaluate
#' @return Result of evaluating expr
#' @keywords internal
#'
#' @details
#' Thread oversubscription occurs when:
#' - Future workers: e.g., 4 parallel R processes
#' - BLAS threads per worker: e.g., 8 threads each
#' - Total threads: 4 * 8 = 32 threads on a 4-core machine!
#'
#' This causes severe performance degradation. The solution is to set
#' BLAS to single-threaded (n=1) when using parallel workers, so:
#' - Future workers: 4 parallel R processes
#' - BLAS threads per worker: 1 thread each
#' - Total threads: 4 * 1 = 4 threads (optimal for 4 cores)
#'
#' @examples
#' \dontrun{
#' # Without management (bad - 32 threads on 4 cores)
#' future::plan(multisession, workers = 4)
#' result <- furrr::future_map(1:100, heavy_blas_function)
#'
#' # With management (good - 4 threads on 4 cores)
#' future::plan(multisession, workers = 4)
#' result <- with_blas_thread_management({
#'   furrr::future_map(1:100, heavy_blas_function)
#' })
#' }
with_blas_thread_management <- function(expr) {
  # Check if we're using parallel processing
  if (!is_future_parallel()) {
    # Sequential execution - no thread management needed
    return(expr)
  }

  # Get current BLAS thread count
  old_threads <- get_blas_threads()

  # Set BLAS to single-threaded for parallel workers
  old_state <- set_blas_threads(1)

  # Execute with proper cleanup
  on.exit({
    restore_blas_threads(old_state)
  }, add = TRUE)

  # Inform user if thread management is active
  if (!is.null(old_state)) {
    message(sprintf(
      "Parallel processing detected: Setting BLAS to single-threaded mode (was %d threads)",
      if (is.numeric(old_state)) old_state else old_threads %||% "unknown"
    ))
  }

  expr
}

# Helper operator for NULL coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x
