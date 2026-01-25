# CLAUDE.md

This file provides guidance for Claude Code when working with the synthdid R package.

## Build & Development Commands

```r
devtools::check()        # Full package check (build, test, docs)
devtools::test()         # Run all tests
testthat::test_file("tests/testthat/test_synthdid.R")  # Run single test file
devtools::load_all()     # Load package for interactive development
devtools::document()     # Regenerate roxygen2 documentation
Rcpp::compileAttributes()  # Regenerate C++ exports after modifying src/
pkgdown::build_site()    # Build documentation website locally
```

## Architecture Overview

### Key Components
- **R/synthdid.R** - Core `synthdid_estimate()` function and related estimators (sc_estimate, did_estimate)
- **R/solver.R** - Frank-Wolfe optimization algorithm for weight computation
- **R/vcov.R** - Standard error methods (bootstrap, jackknife, placebo)
- **R/plot.R** - ggplot2-based visualization functions
- **src/*.cpp** - Rcpp C++ implementations of performance-critical computations

### Data Flow
1. `panel.matrices()` converts long-format panel data to matrix form (Y, N0, T0)
2. `synthdid_estimate()` computes lambda (time weights) and omega (unit weights) via Frank-Wolfe
3. Treatment effect computed as weighted difference-in-differences
4. `vcov()` methods compute standard errors via resampling

## Key Conventions
- Uses `furrr` for parallelization (requires `future::plan()` configuration)
- roxygen2 with markdown for documentation
- testthat edition 3 for testing
- C++11 via Rcpp
