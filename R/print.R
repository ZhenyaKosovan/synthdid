#' Print a synthdid object
#' @param x The object to print
#' @param ... Additional arguments (currently ignored).
#' @return Invisibly returns the original object
#' @examples
#' \donttest{
#' data(california_prop99)
#' setup <- panel.matrices(california_prop99)
#' tau.hat <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
#' print(tau.hat)
#' }
#' @method print synthdid_estimate
#' @export
print.synthdid_estimate <- function(x, ...) {
  cat(format(x, ...), "\n")
  invisible(x)
}
