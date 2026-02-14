#' @description
#' This package implements the synthetic difference in difference estimator (SDID) for the average treatment effect in panel data,
#' as proposed in Arkhangelsky et al (2019). We observe matrices of outcomes Y and binary treatment indicators W
#' that we think of as satisfying Y\[i,j\] = L\[i,j\] + tau\[i,j\] W\[i,j\] + noise\[i,j\].
#' Here tau\[i,j\] is the effect of treatment on the unit i at time j, and we estimate the average effect of
#' treatment when and where it happened: the average of tau\[i,j\] over the observations with W\[i,j\]=1.
#' All treated units must begin treatment simultaneously, so W is a block matrix: W\[i,j\] = 1 for i > N0 and j > T0
#' and zero otherwise, with N0 denoting the number of control units and T0 the number of observation times
#' before onset of treatment. This applies, in particular, to the case of a single treated unit or treated period.
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
