# Package-Level Constants for synthdid
#
# This file defines named constants used throughout the synthdid package.
# Using named constants improves code readability and maintainability by:
# - Making the purpose of numeric values clear
# - Centralizing threshold values for easy tuning
# - Avoiding "magic numbers" scattered throughout the codebase
#
# These constants are internal to the package and not exported.
# They are grouped by functional area:
#
# - Optimization Parameters: Frank-Wolfe solver behavior and convergence criteria
# - Regularization Defaults: Default regularization parameters for weight estimation
# - Standard Error Computation: Bootstrap, jackknife, and placebo SE settings
# - Weight Processing: Sparsification and truncation thresholds
# - Numerical Stability: Tolerances and safeguards
# - Visualization: Default plot parameters
# - Memory Estimation: Computational resource estimation constants


# =============================================================================
# OPTIMIZATION PARAMETERS
# =============================================================================

# Default maximum iterations for Frank-Wolfe optimization
SYNTHDID_MAX_ITER_DEFAULT <- 1e4

# Maximum iterations for pre-sparsification optimization round
SYNTHDID_MAX_ITER_PRE_SPARSIFY <- 100

# Default minimum decrease threshold for solver functions
# Used as stopping criterion: stop if objective decrease < min.decrease^2
SYNTHDID_MIN_DECREASE_DEFAULT <- 1e-3

# Multiplier for noise.level to compute min.decrease in synthdid_estimate
# min.decrease = SYNTHDID_MIN_DECREASE_NOISE_MULTIPLIER * noise.level
SYNTHDID_MIN_DECREASE_NOISE_MULTIPLIER <- 1e-5


# =============================================================================
# REGULARIZATION DEFAULTS
# =============================================================================

# Default eta.lambda (infinitesimal ridge regularization for time weights)
SYNTHDID_ETA_LAMBDA_DEFAULT <- 1e-6

# Default eta.omega for synthetic control estimator
# Uses infinitesimal regularization like eta.lambda
SYNTHDID_ETA_OMEGA_SC_DEFAULT <- 1e-6


# =============================================================================
# STANDARD ERROR COMPUTATION
# =============================================================================

# Default number of replications for bootstrap/placebo standard errors
SYNTHDID_SE_REPLICATIONS_DEFAULT <- 200

# Critical value for 95% confidence intervals (standard normal quantile)
# Corresponds to qnorm(0.975) for two-sided 95% CI
SYNTHDID_CI_Z_95 <- 1.96


# =============================================================================
# WEIGHT PROCESSING
# =============================================================================

# Sparsify threshold divisor
# Weights <= max(weight) / SYNTHDID_SPARSIFY_DIVISOR are set to zero
SYNTHDID_SPARSIFY_DIVISOR <- 4

# Default mass for synthdid_controls table truncation
# Table shows controls accounting for at least this fraction of total weight
SYNTHDID_CONTROLS_MASS_DEFAULT <- 0.9


# =============================================================================
# NUMERICAL STABILITY
# =============================================================================

# Minimum noise level to prevent division by zero
# Used when noise.level is effectively zero
SYNTHDID_MIN_NOISE_LEVEL <- 1e-8

# Generic small epsilon for numerical comparisons
SYNTHDID_EPSILON <- 1e-6


# =============================================================================
# VISUALIZATION PARAMETERS
# =============================================================================

# Default lambda plot scale divisor
# Controls the height of the lambda weight ribbon in plots
SYNTHDID_LAMBDA_PLOT_SCALE_DEFAULT <- 3

# Default curvature for treatment effect arrows in plots
SYNTHDID_EFFECT_CURVATURE_DEFAULT <- 0.3

# Default line width for plot elements
SYNTHDID_LINE_WIDTH_DEFAULT <- 0.5

# Default point size for plot elements
SYNTHDID_POINT_SIZE_DEFAULT <- 1

# Default alpha for trajectory lines
SYNTHDID_TRAJECTORY_ALPHA_DEFAULT <- 0.5

# Default alpha for diff-in-diff diagram elements
SYNTHDID_DIAGRAM_ALPHA_DEFAULT <- 0.95

# Default alpha for treatment effect arrows
SYNTHDID_EFFECT_ALPHA_DEFAULT <- 0.95

# Default alpha for onset vertical lines
SYNTHDID_ONSET_ALPHA_DEFAULT <- 0.3

# Default alpha for confidence interval arrows
SYNTHDID_CI_ALPHA_DEFAULT <- 0.3

# Weight threshold below which units are shown as negligible
# Units with weight <= this value shown as transparent xs in units_plot
SYNTHDID_NEGLIGIBLE_WEIGHT_THRESHOLD <- 0.001

# Alpha for negligible weight units in plots
SYNTHDID_NEGLIGIBLE_ALPHA_DEFAULT <- 0.3

# Spaghetti line width for individual unit trajectories
SYNTHDID_SPAGHETTI_LINE_WIDTH_DEFAULT <- 0.2

# Spaghetti label size for unit names
SYNTHDID_SPAGHETTI_LABEL_SIZE_DEFAULT <- 2

# Spaghetti line alpha
SYNTHDID_SPAGHETTI_LINE_ALPHA_DEFAULT <- 0.3

# Spaghetti label alpha
SYNTHDID_SPAGHETTI_LABEL_ALPHA_DEFAULT <- 0.5


# =============================================================================
# MEMORY ESTIMATION
# =============================================================================

# Bytes per double-precision floating point number
SYNTHDID_BYTES_PER_DOUBLE <- 8

# GB divisor (bytes to gigabytes)
SYNTHDID_GB_DIVISOR <- 1024^3

# Memory warning threshold (GB) - high
# Problems exceeding this may require special handling
SYNTHDID_MEMORY_WARNING_HIGH_GB <- 16

# Memory warning threshold (GB) - moderate
# Problems exceeding this warrant a note to the user
SYNTHDID_MEMORY_WARNING_MODERATE_GB <- 8

# Maximum assumed parallel cores for memory estimation
SYNTHDID_MAX_PARALLEL_CORES <- 8


# =============================================================================
# SUMMARY AND DISPLAY
# =============================================================================

# Default number of digits for weight display in summary output
SYNTHDID_WEIGHT_DIGITS_DEFAULT <- 3

# Default number of top controls/periods to display in summary
SYNTHDID_SUMMARY_TOP_N_DEFAULT <- 5
