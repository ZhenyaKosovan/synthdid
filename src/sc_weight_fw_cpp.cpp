#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]



// [[Rcpp::export(rng = false)]]
NumericVector fw_step_cpp(const NumericMatrix& A,
                          const NumericVector& x,
                          const NumericVector& b,
                          const double eta,
                          const Rcpp::Nullable<double> alpha = R_NilValue) {
  const int n_rows = A.nrow();
  const int n_cols = A.ncol();

  NumericVector Ax(n_rows, 0.0);
  for (int col = 0; col < n_cols; ++col) {
    const double x_col = x[col];
    if (x_col == 0.0) continue;
    for (int row = 0; row < n_rows; ++row) {
      Ax[row] += A(row, col) * x_col;
    }
  }

  NumericVector half_grad(n_cols, 0.0);
  for (int col = 0; col < n_cols; ++col) {
    double acc = 0.0;
    for (int row = 0; row < n_rows; ++row) {
      acc += (Ax[row] - b[row]) * A(row, col);
    }
    half_grad[col] = acc + eta * x[col];
  }

  // index of minimum entry in half_grad (ties choose first)
  int i_min = 0;
  double min_val = half_grad[0];
  for (int col = 1; col < n_cols; ++col) {
    const double val = half_grad[col];
    if (val < min_val) {
      min_val = val;
      i_min = col;
    }
  }

  if (alpha.isNotNull()) {
    const double a = Rcpp::as<double>(alpha);
    NumericVector out = clone(x);
    for (int idx = 0; idx < n_cols; ++idx) {
      out[idx] *= (1.0 - a);
    }
    out[i_min] += a;
    return out;
  }

  NumericVector dx = clone(x);
  for (int idx = 0; idx < n_cols; ++idx) {
    dx[idx] = -dx[idx];
  }
  dx[i_min] = 1.0 - x[i_min];

  bool all_zero = true;
  for (int idx = 0; idx < n_cols; ++idx) {
    if (dx[idx] != 0.0) {
      all_zero = false;
      break;
    }
  }
  if (all_zero) {
    return x;
  }

  NumericVector d_err(n_rows);
  for (int row = 0; row < n_rows; ++row) {
    d_err[row] = A(row, i_min) - Ax[row];
  }

  double sum_err_sq = 0.0;
  for (int row = 0; row < n_rows; ++row) {
    sum_err_sq += d_err[row] * d_err[row];
  }

  double sum_dx_sq = 0.0;
  for (int idx = 0; idx < n_cols; ++idx) {
    sum_dx_sq += dx[idx] * dx[idx];
  }

  double num = 0.0;
  for (int idx = 0; idx < n_cols; ++idx) {
    num += half_grad[idx] * dx[idx];
  }
  double step = -num / (sum_err_sq + eta * sum_dx_sq);
  double constrained_step = std::min(1.0, std::max(0.0, step));

  NumericVector out = clone(x);
  for (int idx = 0; idx < n_cols; ++idx) {
    out[idx] += constrained_step * dx[idx];
  }
  return out;
}

// [[Rcpp::export(rng = false)]]
List sc_weight_fw_cpp(NumericMatrix Y,
                      const double zeta,
                      const bool intercept,
                      NumericVector lambda,
                      const double min_decrease,
                      const int max_iter) {
  const int N0 = Y.nrow();
  const int T0 = Y.ncol() - 1;
  const double min_dec_sq = min_decrease * min_decrease;

  if (lambda.size() == 0) {
    lambda = NumericVector(T0, 1.0 / static_cast<double>(T0));
  }

  if (intercept) {
    for (int col = 0; col < Y.ncol(); ++col) {
      double mean_col = 0.0;
      for (int row = 0; row < N0; ++row) {
        mean_col += Y(row, col);
      }
      mean_col /= static_cast<double>(N0);
      for (int row = 0; row < N0; ++row) {
        Y(row, col) -= mean_col;
      }
    }
  }

  NumericMatrix A(N0, T0);
  NumericVector b(N0);
  for (int row = 0; row < N0; ++row) {
    for (int col = 0; col < T0; ++col) {
      A(row, col) = Y(row, col);
    }
    b[row] = Y(row, T0);
  }

  const double eta = N0 * (zeta * zeta);
  NumericVector vals(max_iter, NumericVector::get_na());

  NumericVector lambda_work = clone(lambda);
  NumericVector err(N0);

  for (int t = 0; t < max_iter; ++t) {
    lambda_work = fw_step_cpp(A, lambda_work, b, eta, R_NilValue);

    std::fill(err.begin(), err.end(), 0.0);
    for (int col = 0; col < T0; ++col) {
      const double lambda_col = lambda_work[col];
      if (lambda_col == 0.0) continue;
      for (int row = 0; row < N0; ++row) {
        err[row] += A(row, col) * lambda_col;
      }
    }
    for (int row = 0; row < N0; ++row) {
      err[row] -= b[row];
    }

    double sum_err_sq = 0.0;
    for (int row = 0; row < N0; ++row) {
      sum_err_sq += err[row] * err[row];
    }

    double sum_lambda_sq = 0.0;
    for (int col = 0; col < T0; ++col) {
      sum_lambda_sq += lambda_work[col] * lambda_work[col];
    }

    vals[t] = (zeta * zeta) * sum_lambda_sq + sum_err_sq / static_cast<double>(N0);

    if (t >= 1) {
      double prev = vals[t - 1];
      double curr = vals[t];
      if (!NumericVector::is_na(prev) && !NumericVector::is_na(curr)) {
        if (!((prev - curr) > min_dec_sq)) {
          break;
        }
      }
    }
  }

  return List::create(
    _["lambda"] = lambda_work,
    _["vals"] = vals
  );
}
