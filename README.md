# synthdid: Synthetic Difference in Differences Estimation

<!-- badges: start -->

[![R-CMD-check](https://github.com/ZhenyaKosovan/synthdid/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ZhenyaKosovan/synthdid/actions/workflows/R-CMD-check.yaml) [![Codecov test coverage](https://codecov.io/gh/ZhenyaKosovan/synthdid/graph/badge.svg)](https://app.codecov.io/gh/ZhenyaKosovan/synthdid)

<!-- badges: end -->

This package implements the synthetic difference in difference estimator (SDID) for the average treatment effect in panel data, as proposed in Arkhangelsky et al. (2019). We observe matrices of outcomes `Y` and binary treatment indicators `W` that satisfy $Y_{ij} = L_{ij} + \tau_{ij} W_{ij} + \varepsilon_{ij}$. Here $\tau_{ij}$ is the effect of treatment on unit $i$ at time $j$, and we estimate the average effect of treatment when and where it happened (the average of $\tau_{ij}$ over the observations with $W_{ij} = 1$). All treated units must begin treatment simultaneously, so $W$ is a block matrix: $W_{ij} = 1$ for $i > N_0$ and $j > T_0$ and zero otherwise, with $N_0$ denoting the number of control units and $T_0$ the number of observation times before onset of treatment. This applies, in particular, to the case of a single treated unit or treated period.

## What's New

**Version 2.0** introduces a modern formula-based interface similar to `lm()` and `plm()`, while maintaining 100% backward compatibility:

-   **Formula interface**: `synthdid(outcome ~ treatment, data, index = c("unit", "time"))`
-   **Standard R methods**: `summary()`, `coef()`, `confint()`, `predict()`, `residuals()`, `fitted()`
-   **Easy method comparison**: Switch between estimators with `update(result, method = "did")`
-   **Enhanced performance**: RcppArmadillo with AVX vectorization for 3-8x speedup

See the [new vignette](https://zhenyakosovan.github.io/synthdid/articles/formula-interface.html) for details.

## Helpful Links

-   The [R package documentation](https://zhenyakosovan.github.io/synthdid/) contains usage examples and method reference.
-   The [formula interface vignette](https://zhenyakosovan.github.io/synthdid/articles/formula-interface.html) demonstrates the new modern interface.
-   The [plotting vignette](https://zhenyakosovan.github.io/synthdid/articles/more-plotting.html) contains a gallery of plot examples.
-   For community questions and answers around usage, see the [GitHub issues page](https://github.com/ZhenyaKosovan/synthdid/issues).

## Installation

The current development version can be installed from source using devtools.

``` r
devtools::install_github("ZhenyaKosovan/synthdid")
```

## Quick Start

### New Formula Interface (Recommended)

``` r
set.seed(12345)
library(synthdid)

# Estimate the effect of California Proposition 99 on cigarette consumption
data("california_prop99")

# Modern formula interface
result <- synthdid(PacksPerCapita ~ treated,
                   data = california_prop99,
                   index = c("State", "Year"),
                   se = TRUE,
                   se_method = "placebo")

# Use standard R methods
print(result)
summary(result)
coef(result)
confint(result)
plot(result)

# Compare estimators easily
did_result <- update(result, method = "did")
sc_result <- update(result, method = "sc")

# Predictions and diagnostics
effect_curve <- predict(result, type = "effect")
counterfactual <- predict(result, type = "counterfactual")
residuals(result, type = "control")
```

### Classic Interface (Still Supported)

The original matrix-based interface continues to work exactly as before:

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

## Key Features

### Multiple Estimators

Switch between different panel data estimators:

``` r
# Synthetic Difference-in-Differences (default)
synthdid_est <- synthdid(outcome ~ treatment, data, index, method = "synthdid")

# Pure Difference-in-Differences
did_est <- synthdid(outcome ~ treatment, data, index, method = "did")

# Synthetic Control
sc_est <- synthdid(outcome ~ treatment, data, index, method = "sc")
```

### Covariate Adjustment

Include time-varying covariates:

``` r
result <- synthdid(PacksPerCapita ~ treated | log_income + unemployment,
                   data = california_prop99,
                   index = c("State", "Year"))
```

### Standard Error Methods

Choose from multiple SE estimation methods:

``` r
# Bootstrap (default, most reliable)
result <- synthdid(outcome ~ treatment, data, index,
                   se = TRUE, se_method = "bootstrap", se_replications = 200)

# Jackknife (faster)
result <- synthdid(outcome ~ treatment, data, index,
                   se = TRUE, se_method = "jackknife")

# Placebo (for single treated unit)
result <- synthdid(outcome ~ treatment, data, index,
                   se = TRUE, se_method = "placebo", se_replications = 100)
```

## Speeding Up Standard Error Computation

Bootstrap and placebo standard errors use `furrr` under the hood. You can enable parallel execution by setting a `future` plan:

``` r
library(future)
library(synthdid)

data("california_prop99")

# Cache current plan (by default it's `sequential`)
old_plan <- future::plan()
on.exit(future::plan(old_plan), add = TRUE)

# Allow R to spawn 4 parallel processes
future::plan(future::multisession, workers = 4)

# Estimate with parallel SE computation
result <- synthdid(PacksPerCapita ~ treated,
                   data = california_prop99,
                   index = c("State", "Year"),
                   se = TRUE,
                   se_method = "bootstrap",
                   se_replications = 200)

print(result)
summary(result)
```

Note: on some platforms (e.g., CRAN macOS/Windows builders) multisession may be restricted; in that case `future::plan()` will fall back to sequential execution.

## Performance Optimizations

This package uses **RcppArmadillo** with AVX vectorization for significant performance improvements:

-   **3-8x faster** matrix-vector operations using optimized BLAS
-   **2-4x faster** tensor operations for covariate adjustment
-   Automatic use of AVX128/256 SIMD instructions on modern CPUs
-   Zero overhead for formula interface

## Package Interface Comparison

| Task | Old Interface | New Interface |
|------------------|---------------------------|---------------------------|
| Basic estimation | `setup <- panel.matrices(data)`<br>`synthdid_estimate(setup$Y, setup$N0, setup$T0)` | `synthdid(outcome ~ treatment, data, index)` |
| Get coefficient | `c(result)` | `coef(result)` |
| Summary | Custom function | `summary(result)` |
| Confidence intervals | Manual calculation | `confint(result)` |
| Predictions | Custom calculation | `predict(result, type = "counterfactual")` |
| Method comparison | Re-run with different function | `update(result, method = "did")` |

## Documentation

For detailed examples and use cases, see:

-   `vignette("formula-interface")` - Modern formula-based interface
-   `vignette("synthdid")` - Original introduction
-   `vignette("more-plotting")` - Plotting gallery
-   `?synthdid` - Function documentation

## References

Dmitry Arkhangelsky, Susan Athey, David A. Hirshberg, Guido W. Imbens, and Stefan Wager. **Synthetic Difference in Differences**, 2019. [arXiv](https://arxiv.org/abs/1812.09970)

## Contributing

Contributions are welcome! Please see the [GitHub repository](https://github.com/ZhenyaKosovan/synthdid) for more information.
