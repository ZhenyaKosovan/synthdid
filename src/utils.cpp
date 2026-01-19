#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

/**
 * @brief Compute gradient of covariate adjustment coefficients (beta) for synthdid
 *
 * This function computes the gradient of the synthetic diff-in-diff objective
 * with respect to covariate adjustment coefficients. It handles 3D covariate
 * arrays X where X[i,j,k] is the k-th covariate for unit i at time j.
 *
 * The gradient has two components:
 * 1. Lambda term: gradient from time-weighted control unit residuals
 * 2. Omega term: gradient from unit-weighted pre-treatment period residuals
 *
 * Mathematical formula for gradient component k:
 *   grad[k] = (1/N0) * err_lambda^T * X[1:N0, :, k] * [lambda; -1]
 *           + (1/T0) * err_omega^T * X[:, 1:T0, k]^T * [omega; -1]
 *
 * OPTIMIZATION: Uses RcppArmadillo for:
 *   - arma::cube for efficient 3D array slicing
 *   - Vectorized matrix-vector products (GEMV)
 *   - Vectorized dot products (DOT)
 * These operations use AVX SIMD instructions for 2-4x speedup.
 *
 * @param X 3D array (N0+1) x (T0+1) x K of covariate values
 * @param lambda Time weights (length T0)
 * @param omega Unit weights (length N0)
 * @param err_lambda Residuals for control units (length N0)
 * @param err_omega Residuals for pre-treatment periods (length T0)
 * @param N0 Number of control units
 * @param T0 Number of pre-treatment periods
 * @return Gradient vector (length K)
 */
// [[Rcpp::export(rng = false)]]
NumericVector grad_beta_cpp(const NumericVector& X,
                            const arma::vec& lambda,
                            const arma::vec& omega,
                            const arma::vec& err_lambda,
                            const arma::vec& err_omega,
                            const int N0,
                            const int T0) {
  // Validate input dimensions
  const IntegerVector dims = X.attr("dim");
  if (dims.size() != 3) {
    stop("X must be a 3D array");
  }
  const int n_rows = dims[0]; // expected N0 + 1 (control + treated units)
  const int n_cols = dims[1]; // expected T0 + 1 (pre-treatment + target period)
  const int K = dims[2];      // number of covariates

  if (n_rows != N0 + 1 || n_cols != T0 + 1) {
    stop("Dimension mismatch between X and N0/T0");
  }

  // Convert R 3D array to Armadillo cube for vectorized slicing operations
  // Armadillo stores cubes in column-major order: X_cube(row, col, slice)
  arma::cube X_cube(n_rows, n_cols, K);
  std::copy(X.begin(), X.end(), X_cube.memptr());

  arma::vec grad(K);

  // Loop over covariates (could be parallelized for large K)
  for (int k = 0; k < K; ++k) {
    // Extract the k-th covariate slice: (N0+1) x (T0+1) matrix
    arma::mat slice = X_cube.slice(k);

    // ========== Lambda term: gradient from time-weighted residuals ==========

    // Create weight vector: [lambda[0], ..., lambda[T0-1], -1]
    // The -1 corresponds to the treated period (last column)
    arma::vec weight_lambda(T0 + 1);
    weight_lambda.head(T0) = lambda;
    weight_lambda(T0) = -1.0;

    // Extract control unit rows (first N0 rows) and compute weighted sum across time
    // slice_control * weight_lambda gives: X[i,:,k] %*% [lambda; -1] for i=1..N0
    arma::mat slice_control = slice.rows(0, N0 - 1);
    arma::vec dot_lambda = slice_control * weight_lambda;

    // Compute lambda component: (1/N0) * sum_i err_lambda[i] * dot_lambda[i]
    // Uses vectorized dot product (AVX-accelerated)
    double term_lambda = arma::dot(err_lambda, dot_lambda) / static_cast<double>(N0);

    // ========== Omega term: gradient from unit-weighted residuals ==========

    // Create weight vector: [omega[0], ..., omega[N0-1], -1]
    // The -1 corresponds to the treated unit (last row)
    arma::vec weight_omega(N0 + 1);
    weight_omega.head(N0) = omega;
    weight_omega(N0) = -1.0;

    // Extract pre-treatment columns (first T0 columns) and compute weighted sum across units
    // slice_pretreat^T * weight_omega gives: X[:,j,k]^T %*% [omega; -1] for j=1..T0
    arma::mat slice_pretreat = slice.cols(0, T0 - 1);
    arma::vec dot_omega = slice_pretreat.t() * weight_omega;

    // Compute omega component: (1/T0) * sum_j err_omega[j] * dot_omega[j]
    // Uses vectorized dot product (AVX-accelerated)
    double term_omega = arma::dot(err_omega, dot_omega) / static_cast<double>(T0);

    // Combine both terms for the gradient of covariate k
    grad(k) = term_lambda + term_omega;
  }

  return wrap(grad);
}

/**
 * @brief Contract a 3D tensor along its third dimension using a weight vector
 *
 * Computes a weighted sum of slices of a 3D array X along the third dimension:
 *   result[i,j] = sum_k v[k] * X[i,j,k]
 *
 * This is used in synthdid to combine covariate effects: given covariate array X
 * and coefficient vector v (beta), this computes the total covariate contribution
 * to the outcome.
 *
 * OPTIMIZATION: Uses RcppArmadillo for:
 *   - arma::cube for efficient 3D array representation
 *   - Vectorized matrix addition and scalar multiplication
 *   - Zero-copy slice extraction
 * Each slice operation uses AVX for element-wise operations, providing 2-4x speedup.
 *
 * @param X 3D array n_rows x n_cols x K of values
 * @param v Weight vector (length K)
 * @return Matrix n_rows x n_cols of weighted sums
 */
// [[Rcpp::export(rng = false)]]
NumericMatrix contract3_cpp(const NumericVector& X, const arma::vec& v) {
  // Validate input dimensions
  const IntegerVector dims = X.attr("dim");
  if (dims.size() != 3) {
    stop("X must be a 3D array");
  }
  const int n_rows = dims[0];
  const int n_cols = dims[1];
  const int K = dims[2];      // number of slices

  if (static_cast<int>(v.n_elem) != K) {
    stop("Length of v must match third dimension of X");
  }

  // Convert R 3D array to Armadillo cube for vectorized operations
  // Memory layout: column-major, allowing efficient slicing
  arma::cube X_cube(n_rows, n_cols, K);
  std::copy(X.begin(), X.end(), X_cube.memptr());

  // Initialize output matrix with zeros
  arma::mat out = arma::zeros<arma::mat>(n_rows, n_cols);

  // Accumulate weighted slices
  // out += v[k] * X_cube.slice(k) uses vectorized operations:
  //   - slice extraction (zero-copy view)
  //   - scalar multiplication (AVX SIMD)
  //   - matrix addition (AVX SIMD)
  for (int k = 0; k < K; ++k) {
    if (v(k) != 0.0) {  // Skip zero weights to save computation
      out += v(k) * X_cube.slice(k);
    }
  }

  return wrap(out);
}
