#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export(rng = false)]]
NumericVector grad_beta_cpp(const NumericVector& X,
                            const NumericVector& lambda,
                            const NumericVector& omega,
                            const NumericVector& err_lambda,
                            const NumericVector& err_omega,
                            const int N0,
                            const int T0) {
  const IntegerVector dims = X.attr("dim");
  if (dims.size() != 3) {
    stop("X must be a 3D array");
  }
  const int n_rows = dims[0]; // expected N0 + 1
  const int n_cols = dims[1]; // expected T0 + 1
  const int K = dims[2];

  if (n_rows != N0 + 1 || n_cols != T0 + 1) {
    stop("Dimension mismatch between X and N0/T0");
  }

  NumericVector grad(K);
  const int slice_size = n_rows * n_cols;

  for (int k = 0; k < K; ++k) {
    const int offset = k * slice_size;

    double term_lambda = 0.0;
    for (int i = 0; i < N0; ++i) {
      double dot = 0.0;
      for (int j = 0; j <= T0; ++j) {
        const double weight = (j < T0) ? lambda[j] : -1.0;
        dot += X[offset + i + n_rows * j] * weight;
      }
      term_lambda += err_lambda[i] * dot;
    }
    term_lambda /= static_cast<double>(N0);

    double term_omega = 0.0;
    for (int j = 0; j < T0; ++j) {
      double dot = 0.0;
      for (int i = 0; i <= N0; ++i) {
        const double weight = (i < N0) ? omega[i] : -1.0;
        dot += X[offset + i + n_rows * j] * weight;
      }
      term_omega += err_omega[j] * dot;
    }
    term_omega /= static_cast<double>(T0);

    grad[k] = term_lambda + term_omega;
  }

  return grad;
}

// [[Rcpp::export(rng = false)]]
NumericMatrix contract3_cpp(const NumericVector& X, const NumericVector& v) {
  const IntegerVector dims = X.attr("dim");
  if (dims.size() != 3) {
    stop("X must be a 3D array");
  }
  const int n_rows = dims[0];
  const int n_cols = dims[1];
  const int K = dims[2];

  if (v.size() != K) {
    stop("Length of v must match third dimension of X");
  }

  NumericMatrix out(n_rows, n_cols);
  if (K == 0) {
    return out;
  }

  const int slice_size = n_rows * n_cols;
  for (int k = 0; k < K; ++k) {
    const double weight = v[k];
    if (weight == 0.0) continue;
    const int offset = k * slice_size;
    for (int j = 0; j < n_cols; ++j) {
      for (int i = 0; i < n_rows; ++i) {
        out(i, j) += weight * X[offset + i + n_rows * j];
      }
    }
  }

  return out;
}
