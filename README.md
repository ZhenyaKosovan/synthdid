
# synthdid: Synthetic Difference in Differences Estimation

<!-- badges: start -->

[![R-CMD-check](https://github.com/ZhenyaKosovan/synthdid/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ZhenyaKosovan/synthdid/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/ZhenyaKosovan/synthdid/graph/badge.svg)](https://app.codecov.io/gh/ZhenyaKosovan/synthdid)

<!-- badges: end -->

This package implements the synthetic difference in difference estimator
(SDID) for the average treatment effect in panel data, as proposed in
Arkhangelsky et al. (2019). We observe matrices of outcomes `Y` and
binary treatment indicators `W` that satisfy
$Y_{ij} = L_{ij} + \tau_{ij} W_{ij} + \varepsilon_{ij}$. Here
$\tau_{ij}$ is the effect of treatment on unit $i$ at time $j$, and we
estimate the average effect of treatment when and where it happened (the
average of $\tau_{ij}$ over the observations with $W_{ij} = 1$). All
treated units must begin treatment simultaneously, so $W$ is a block
matrix: $W_{ij} = 1$ for $i > N_0$ and $j > T_0$ and zero otherwise,
with $N_0$ denoting the number of control units and $T_0$ the number of
observation times before onset of treatment. This applies, in
particular, to the case of a single treated unit or treated period.

This package is currently in beta and the functionality and interface
are subject to change.

Some helpful links for getting started:

- The [R package
  documentation](https://zhenyakosovan.github.io/synthdid/) contains
  usage examples and method reference.
- The [online
  vignettes](https://zhenyakosovan.github.io/synthdid/articles/more-plotting.html)
  contain a gallery of plot examples.
- For community questions and answers around usage, see the [GitHub
  issues page](https://github.com/ZhenyaKosovan/synthdid/issues).

## Installation

The current development version can be installed from source using
devtools.

``` r
devtools::install_github("ZhenyaKosovan/synthdid")
```

## Example

``` r
library(synthdid)

# Estimate the effect of California Proposition 99 on cigarette consumption
data("california_prop99")
setup <- panel.matrices(california_prop99)
tau_hat <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
# Note: SE estimation requires re-estimation of the model N-replications times. It can be time consuming!
se <- sqrt(vcov(tau_hat, method = "placebo"))
sprintf("point estimate: %1.2f", tau_hat)
sprintf("95%% CI (%1.2f, %1.2f)", tau_hat - 1.96 * se, tau_hat + 1.96 * se)
plot(tau_hat)
```

## Speeding up standard error computation

Bootstrap and placebo standard errors use `furrr` under the hood. You
can enable parallel execution by setting a `future` plan. The snippet
below uses multiple cores via `multisession` when supported and resets
to sequential afterward.

``` r
library(future)
library(synthdid)

data("california_prop99")
setup <- panel.matrices(california_prop99)

# cache current plan (by default it's `sequential`)
old_plan <- future::plan()
on.exit(future::plan(old_plan), add = TRUE)

# allow R to spawn 4 parallel processes, please adjust it if you have a machine with small number of cores.
future::plan(future::multisession, workers = 4)

tau_hat <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
se <- sqrt(vcov(tau_hat, method = "placebo"))
sprintf("point estimate: %1.2f", tau_hat)
sprintf("95%% CI (%1.2f, %1.2f)", tau_hat - 1.96 * se, tau_hat + 1.96 * se)
plot(tau_hat)
```

Note: on some platforms (e.g., CRAN macOS/Windows builders) multisession
may be restricted; in that case `future::plan()` will fall back to
sequential execution.

### References

Dmitry Arkhangelsky, Susan Athey, David A. Hirshberg, Guido W. Imbens,
and Stefan Wager. **Synthetic Difference in Differences**, 2019.
[arXiv](https://arxiv.org/abs/1812.09970)
