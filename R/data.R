#' California proposition 99
#'
#' A dataset containing per-capita cigarette consumption (in packs).
#' In year 1989 California imposed a Tobacco tax. The column `treated` is 1 from then on for California.
#'
#' @docType data
#' @name california_prop99
#'
#' @format A data frame with 1209 rows and 4 variables:
#' \describe{
#'   \item{State}{US state name, character string}
#'   \item{Year}{Year, integer}
#'   \item{PacksPerCapita}{per-capita cigarette consumption, numeric}
#'   \item{treated}{the treatmed indicator 0: control, 1: treated, numeric}
#' }
#' @source Abadie, Alberto, Alexis Diamond, and Jens Hainmueller.
#'  "Synthetic control methods for comparative case studies: Estimating the effect of California’s tobacco control program."
#'   Journal of the American statistical Association 105, no. 490 (2010): 493-505.
#'
#' @usage data(california_prop99)
#'
#' @examples
#' \donttest{
#' # Load tobacco sales in long panel format.
#' data("california_prop99")
#' # Transform to N*T matrix format required for synthdid,
#' # where N is the number of units and T the time periods.
#' setup <- panel.matrices(california_prop99)
#' }
#'
NULL

#' PENN
#'
#' @docType data
#' @name PENN
#'
#' @format A data frame with 3219 rows and 5 variables.
#' \describe{
#'   \item{country}{country}
#'   \item{year}{year}
#'   \item{log_gdp}{log_gdp}
#'   \item{dem}{dem}
#'   \item{educ}{educ}
#' }
#'
#' @usage data(PENN)
#'
NULL

#' CPS
#'
#' @docType data
#' @name CPS
#'
#' @format A data frame with 2000 rows and 8 variables.
#' \describe{
#'   \item{state}{state}
#'   \item{year}{year}
#'   \item{log_wage}{log_wage}
#'   \item{hours}{hours}
#'   \item{urate}{urate}
#'   \item{min_wage}{min_wage}
#'   \item{open_carry}{open_carry}
#'   \item{abort_ban}{abort_ban}
#' }
#'
#' @usage data(CPS)
#'
NULL


#' Castle Doctrine / Stand Your Ground Laws
#'
#' State-level panel data on homicide rates and the adoption of Castle Doctrine
#' (Stand Your Ground) laws across US states, 2000--2010. Used by Cheng and
#' Hoekstra (2013) to study whether strengthening self-defense law deters crime
#' or escalates violence. This is a classic example of staggered treatment
#' adoption: 21 states adopted castle doctrine laws at different times during
#' the panel.
#'
#' @docType data
#' @name castle_doctrine
#'
#' @format A data frame with 550 rows (50 states x 11 years) and 5 variables:
#' \describe{
#'   \item{state}{Numeric state identifier (1--50).}
#'   \item{year}{Year (2000--2010).}
#'   \item{l_homicide}{Log homicide rate per 100,000 population.}
#'   \item{post}{Binary treatment indicator: 1 if the state's castle doctrine
#'     law is in effect (i.e., \code{year > effyear}), 0 otherwise.}
#'   \item{effyear}{Year in which the state's castle doctrine law became
#'     effective. \code{NA} for the 29 never-treated states.}
#' }
#'
#' @source Cheng, Cheng and Mark Hoekstra. "Does Strengthening Self-Defense Law
#'   Deter Crime or Escalate Violence? Evidence from Expansions to Castle Doctrine."
#'   \emph{Journal of Human Resources} 48, no. 3 (2013): 821--854.
#'
#'   Data obtained from the \code{did2s} R package by Kyle Butts.
#'
#' @usage data(castle_doctrine)
#'
#' @examples
#' \donttest{
#' data(castle_doctrine)
#'
#' # Staggered SDID estimate
#' est <- synthdid(l_homicide ~ post,
#'   data = castle_doctrine,
#'   index = c("state", "year")
#' )
#' print(est)
#' plot(est)
#' }
NULL


#' Sentencing Enhancement Laws (Abrams 2012)
#'
#' State-level panel data on armed robbery rates and the adoption of sentencing
#' enhancement ("add-on") laws for firearm offenses across US states,
#' 1965--2002. Originally studied by Abrams (2012) and used by Porreca (2022)
#' to demonstrate the staggered synthetic difference-in-differences estimator.
#'
#' This is a balanced panel of 40 states over 38 years. Five states with
#' incomplete outcome data and five states that adopted enhancements before
#' 1970 (treated from the first panel period) are excluded from the original
#' 50-state dataset. Treatment is coded as absorbing: once a state enacts
#' an add-on gun law, it remains treated for all subsequent periods.
#'
#' Twenty-three states adopted sentencing enhancements at 14 different times
#' between 1972 and 1996, making this a rich example of staggered adoption
#' with many cohorts. Seventeen states never adopted.
#'
#' @docType data
#' @name abrams_sentencing
#'
#' @format A data frame with 1520 rows (40 states x 38 years) and 5 variables:
#' \describe{
#'   \item{state}{Numeric FIPS state code.}
#'   \item{year}{Year (1965--2002).}
#'   \item{pcrrobgun}{Per-capita armed (gun) robbery rate per 100,000
#'     population.}
#'   \item{post}{Binary treatment indicator: 1 if the state's sentencing
#'     enhancement is in effect (\code{year > addonyr}), 0 otherwise.}
#'   \item{addonyr}{Year the sentencing enhancement was enacted.
#'     \code{NA} for the 17 never-treated states.}
#' }
#'
#' @source Abrams, David S. "Estimating the Deterrent Effect of Incarceration
#'   using Sentencing Enhancements." \emph{American Economic Journal: Applied
#'   Economics} 4, no. 4 (2012): 32--56.
#'
#'   Replication data from openICPSR project 113838.
#'
#' @references
#' Porreca, Zachary. "Synthetic Difference in Differences Estimation with
#'   Staggered Treatment Timing." \emph{Economics Letters} 220 (2022): 110873.
#'
#' @usage data(abrams_sentencing)
#'
#' @examples
#' \donttest{
#' data(abrams_sentencing)
#'
#' # Replicate Porreca (2022) Table 1: SynthDiD = -16.697
#' est <- synthdid(pcrrobgun ~ post,
#'   data = abrams_sentencing,
#'   index = c("state", "year"),
#'   control_type = "never_treated",
#'   weights = "treated_unit_time"
#' )
#' print(est)
#' coef(est)  # -16.697
#' }
NULL
