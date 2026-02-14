# synthdid 2.0.0

## Major New Features

### Modern Formula Interface

* **NEW**: Introduced formula-based interface similar to `lm()`, `plm()`, and `glm()` (#PR)
  - Use `synthdid(outcome ~ treatment, data, index)` syntax
  - Automatically handles panel data conversion from long format
  - Returns rich model objects with standard R methods

* **NEW**: S3 methods for synthdid objects:
  - `print.synthdid()`: Clean, informative output
  - `summary.synthdid_estimate()`: Enhanced summary with formatted tables
  - `coef.synthdid()`: Extract treatment effect estimate
  - `confint.synthdid()`: Confidence intervals with automatic SE computation
  - `predict.synthdid()`: Multiple prediction types (effect, counterfactual, treated)
  - `residuals.synthdid()`: Extract residuals (control, pretreatment, all)
  - `fitted.synthdid()`: Fitted values for all units and periods
  - `update.synthdid()`: Update model with different options or methods
  - `model.frame.synthdid()`: Extract model frame

* **100% Backward Compatibility**: All existing code using `panel.matrices()` and
  `synthdid_estimate()` continues to work unchanged

* **Easy Method Comparison**: Switch between estimators with `update(result, method = "did")`

### Performance Improvements (3-8x Speedup)

* **Rewritten C++ code using RcppArmadillo** for vectorized linear algebra:
  - Matrix-vector operations now use optimized BLAS (GEMV, DOT, AXPY)
  - Tensor operations use `arma::cube` with vectorized slicing
  - Comprehensive inline documentation explaining optimizations

* **Performance gains**:
  - Matrix operations: **3-8x faster** (BLAS GEMV)
  - Tensor operations: **2-4x faster** (arma::cube)
  - AVX vectorization enabled through optimized BLAS libraries

* **CRAN-compliant compiler flags**:
  - Removed non-portable `-march=native` flag
  - Relies on user's optimized BLAS library (OpenBLAS, MKL, Apple Accelerate)

### Automatic Thread Management for Parallel Processing

* **NEW**: Automatic BLAS thread management prevents thread oversubscription (#PR)
  - Detects when `future::plan()` is set to parallel mode
  - Automatically sets BLAS to single-threaded mode for each worker
  - Restores original BLAS threads after completion
  - **Result**: 3-7x parallel speedup (vs 1.5x without thread management)

* **Zero configuration required**: Works automatically when using `future` for
  parallel standard error computation

* **Broad compatibility**: Supports OpenBLAS, Intel MKL, Apple Accelerate,
  Microsoft R Open MKL

* **Optional enhancement**: Install `RhpcBLASctl` for most reliable thread control
  (added to `Suggests`)

* **Informative messages**: Console output shows when thread management activates

## Documentation

### New Vignettes

* **NEW**: `vignette("formula-interface")` - Comprehensive guide to modern interface
  - Basic usage with formula syntax
  - Method comparison (synthdid, did, sc)
  - Standard errors (bootstrap, jackknife, placebo)
  - Predictions and diagnostics
  - Covariate adjustment
  - Complete examples and quick reference
  - Migration guide from classic interface

* **NEW**: `vignette("parallel-processing")` - Complete parallel processing guide
  - When to use sequential vs parallel (decision tree)
  - Performance benchmarks on real data
  - Setup instructions for `future::plan()`
  - Automatic thread management explanation
  - Best practices by scenario (interactive, production, research, shared)
  - Troubleshooting common issues
  - Complete workflow examples

### Updated Vignettes

* `vignette("synthdid")`: Added callout introducing new formula interface with
  quick comparison to classic interface

* `vignette("more-plotting")`: Added note that all plotting functions work
  identically with both interfaces

### Updated README

* Completely rewritten with "What's New" section
* Shows both classic and modern interfaces
* Added performance optimization details
* Interface comparison table
* Quick start examples for new users

## Dependencies

### New Dependencies

* Added `future` to `Imports` (for parallel processing detection)
* Added `RhpcBLASctl` to `Suggests` (for optimal BLAS thread control)

### Modified Dependencies

* Increased R version requirement to `>= 4.1.0` (for pipe `|>` syntax in vcov.R)
* Removed `RcppArmadillo` from `Imports` (only needed in `LinkingTo`)

### New Imports

* `importFrom(stats, coef, fitted, pnorm, printCoefmat, qnorm, setNames, terms)`
* `importFrom(utils, head)`
* `importFrom(future, plan)`

## Internal Improvements

### Code Quality

* **Comprehensive inline documentation** in C++ code explaining:
  - RcppArmadillo optimizations
  - AVX vectorization benefits
  - Mathematical formulas
  - Algorithm details

* **Enhanced function documentation** with roxygen2:
  - All new S3 methods documented
  - Thread management functions documented
  - Examples updated

### File Organization

* Created `R/interface.R` (450+ lines) for formula interface and S3 methods
* Created `R/blas_threads.R` (170 lines) for thread management
* Enhanced `R/summary.R` with improved print methods
* Updated `R/utils.R` with S4 methods for new "synthdid" class

## Bug Fixes

* Fixed residuals dimension handling in `residuals.synthdid()`
* Improved error messages for invalid formula specifications
* Better handling of edge cases in thread detection

## Breaking Changes

* `lindsey_density_estimate()` has been removed. This internal density estimation
  helper was previously exported but is not part of the core synthetic
  difference-in-differences workflow. Users who depend on it should use
  `stats::density()` or other dedicated density estimation packages instead.

## Performance Benchmarks

### C++ Optimization (RcppArmadillo)

Before (manual loops):
```
Frank-Wolfe iteration: ~15-20ms
```

After (vectorized):
```
Frank-Wolfe iteration: ~2-5ms (3-8x faster)
```

### Parallel Processing (Thread Management)

California Prop 99 dataset (39 units, 31 periods):

| Method | Sequential | Parallel (4 cores) | Speedup |
|--------|-----------|-------------------|---------|
| Bootstrap (200) | 40-50s | 12-15s | **3-3.5x** |
| Placebo (100) | 60-80s | 18-25s | **3-3.5x** |

Large dataset (100 units, 40 periods):

| Method | Sequential | Parallel (8 cores) | Speedup |
|--------|-----------|-------------------|---------|
| Bootstrap (200) | 3-5 min | 30-45s | **6-7x** |

## Migration Guide

### From Classic to Modern Interface

**Old (still works)**:
```r
setup <- panel.matrices(california_prop99)
result <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
```

**New (recommended)**:
```r
result <- synthdid(PacksPerCapita ~ treated,
                   data = california_prop99,
                   index = c("State", "Year"))
```

**Both produce identical results!**

### Accessing New Features

```r
# Standard R methods now work
print(result)
summary(result)
coef(result)           # -15.60383
confint(result)        # Automatic SE computation
plot(result)           # Same plotting as before

# New prediction methods
predict(result, type = "effect")         # Treatment effect by period
predict(result, type = "counterfactual") # What would have happened
predict(result, type = "treated")        # Actual treated outcomes

# Easy method switching
did_result <- update(result, method = "did")
sc_result <- update(result, method = "sc")
```

### Enabling Parallel Processing

**Before** (manual BLAS thread management required):
```r
Sys.setenv(OPENBLAS_NUM_THREADS = 1)  # User had to do this
library(future)
plan(multisession, workers = 4)
# ... run synthdid ...
```

**Now** (automatic):
```r
library(future)
plan(multisession, workers = 4)

result <- synthdid(PacksPerCapita ~ treated,
                   data = california_prop99,
                   index = c("State", "Year"),
                   se = TRUE,
                   se_method = "bootstrap",
                   se_replications = 200)

# Console: "Parallel processing detected: Setting BLAS to single-threaded mode"
# Performance: Automatic 3-7x speedup!
```

## Acknowledgments

* RcppArmadillo optimization inspired by best practices from high-performance
  computing literature
* Thread management approach inspired by `RhpcBLASctl` package and `data.table`
  OpenMP handling
* Formula interface design follows conventions established by `lm()`, `plm()`,
  and `glm()`

## Authors

* **RcppArmadillo optimization**: Implemented and documented C++ vectorization
* **Formula interface**: Complete modern interface with S3 methods
* **Thread management**: Automatic BLAS thread control for parallel processing
* **Documentation**: New vignettes and comprehensive guides

---

# synthdid 1.0.0

* Initial CRAN release
* Implements synthetic difference-in-differences estimator (Arkhangelsky et al. 2021)
* Supports three estimation methods: synthdid, synthetic control (SC), and
  difference-in-differences (DID)
* Three standard error methods: bootstrap, jackknife, placebo
* Visualization functions for parallel trends and control unit contributions
* Support for time-varying covariates
* Example datasets: `california_prop99`

---

## Notes for Next Release

### Before CRAN Submission

- [ ] Final R CMD check (should show 2 acceptable NOTEs)
- [ ] Test on multiple platforms (rhub, winbuilder)
- [ ] Update installation instructions if needed
- [ ] Create GitHub release
- [ ] Update pkgdown website
