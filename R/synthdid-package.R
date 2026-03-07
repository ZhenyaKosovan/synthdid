#' @description
#' This package implements the synthetic difference in difference estimator (SDID) for the average treatment effect in panel data,
#' as proposed in Arkhangelsky et al (2019). We observe matrices of outcomes Y and binary treatment indicators W
#' that we think of as satisfying Y\[i,j\] = L\[i,j\] + tau\[i,j\] W\[i,j\] + noise\[i,j\].
#' Here tau\[i,j\] is the effect of treatment on the unit i at time j, and we estimate the average effect of
#' treatment when and where it happened: the average of tau\[i,j\] over the observations with W\[i,j\]=1.
#'
#' The package supports two treatment adoption patterns:
#' * **Simultaneous adoption**: all treated units begin treatment at the same time (classic SDID).
#' * **Staggered adoption**: treated units adopt at different times. The staggered estimator
#'   decomposes the panel into cohort-specific 2x2 SDID subproblems and aggregates, following
#'   the approach of Porreca (2022).
#'
#' The formula interface [synthdid()] auto-detects the adoption pattern and routes to the
#' appropriate estimator. Use `adoption = "staggered"` or `adoption = "simultaneous"` to
#' override auto-detection.
#'
#' Some helpful links for getting started:
#'
#' * The [R package documentation](https://synth-inference.github.io/synthdid/) contains usage examples and method reference.
#' * The [online vignettes](https://synth-inference.github.io/synthdid/articles/more-plotting.html) contains a gallery of plot examples.
#' * For community questions and answers around usage, see [Github issues page](https://github.com/synth-inference/synthdid/issues).
#'
#' @examples
#' \donttest{
#' # Estimate the effect of California Proposition 99 on cigarette consumption
#' data("california_prop99")
#' fit <- synthdid(PacksPerCapita ~ treated,
#'   data = california_prop99,
#'   index = c("State", "Year")
#' )
#' summary(fit)
#' plot(fit)
#' }
#'
#' @keywords internal
#' @useDynLib synthdid, .registration = TRUE
## usethis namespace: start
#' @importFrom lifecycle deprecated
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
"_PACKAGE"
