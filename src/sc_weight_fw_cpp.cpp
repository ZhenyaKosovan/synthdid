#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

/**
 * @brief Internal C++ Frank-Wolfe step (no R/Rcpp marshalling overhead)
 *
 * Pure Armadillo implementation called directly from sc_weight_fw_cpp.
 * Performs one iteration of the Frank-Wolfe algorithm to solve:
 *   min_(x >= 0, sum(x) = 1) 0.5 * ||A*x - b||^2 + 0.5 * eta * ||x||^2
 */
static arma::vec fw_step_internal(const arma::mat& A,
                                  const arma::vec& x,
                                  const arma::vec& b,
                                  const double eta) {
  // Compute A*x using optimized BLAS GEMV
  arma::vec Ax = A * x;

  // Compute gradient: A^T * (A*x - b) + eta * x
  arma::vec half_grad = A.t() * (Ax - b) + eta * x;

  // Find simplex vertex minimizing <gradient, vertex>
  arma::uword i_min = half_grad.index_min();

  // Compute direction: dx = e_{i_min} - x
  arma::vec dx = -x;
  dx(i_min) = 1.0 - x(i_min);

  // Check for convergence: if already at vertex, return
  if (arma::all(dx == 0.0)) {
    return x;
  }

  // Compute change in prediction error along direction dx
  arma::vec d_err = A.col(i_min) - Ax;

  // Compute squared norms for line search
  double sum_err_sq = arma::dot(d_err, d_err);
  double sum_dx_sq = arma::dot(dx, dx);

  // Line search: optimal step t* = -<gradient, dx> / (<A*dx, A*dx> + eta*<dx, dx>)
  double num = arma::dot(half_grad, dx);
  double step = -num / (sum_err_sq + eta * sum_dx_sq);

  // Project step to [0, 1]
  double constrained_step = std::min(1.0, std::max(0.0, step));

  return x + constrained_step * dx;
}

/**
 * @brief R-callable Frank-Wolfe step (Rcpp export wrapper)
 *
 * Thin wrapper around fw_step_internal for R callers.
 * Supports an optional fixed step size alpha.
 *
 * @param A n_rows x n_cols matrix of control unit outcomes
 * @param x Current weight vector (simplex: x >= 0, sum(x) = 1)
 * @param b Target vector to approximate
 * @param eta Ridge regularization parameter
 * @param alpha Optional fixed step size (if provided, skip line search)
 * @return Updated weight vector after one Frank-Wolfe step
 */
// [[Rcpp::export(rng = false)]]
NumericVector fw_step_cpp(const arma::mat& A,
                          const arma::vec& x,
                          const arma::vec& b,
                          const double eta,
                          const Rcpp::Nullable<double> alpha = R_NilValue) {
  // If fixed step size provided, use it directly (no line search)
  if (alpha.isNotNull()) {
    const double a = Rcpp::as<double>(alpha);
    arma::vec half_grad = A.t() * (A * x - b) + eta * x;
    arma::uword i_min = half_grad.index_min();
    arma::vec out = (1.0 - a) * x;
    out(i_min) += a;
    return wrap(out);
  }

  return wrap(fw_step_internal(A, x, b, eta));
}

/**
 * @brief Compute synthetic control weights using Frank-Wolfe algorithm
 *
 * Solves the penalized least squares problem:
 *   min_(lambda >= 0, sum(lambda) = 1) ||Y * [lambda; -1]||^2 / N0 + zeta^2 * ||lambda||^2
 *
 * Uses fw_step_internal directly (no R round-trip per iteration).
 *
 * @param Y N0 x (T0+1) matrix of control outcomes (last column is target period)
 * @param zeta Ridge penalty parameter (larger = more regularization)
 * @param intercept If true, demean columns before fitting
 * @param lambda Initial weights (if empty, starts with uniform weights)
 * @param min_decrease Convergence threshold for objective function decrease
 * @param max_iter Maximum number of Frank-Wolfe iterations
 * @return List with components:
 *   - lambda: Optimal time weights
 *   - vals: Objective function values at each iteration
 *   - converged: Boolean indicating if optimization converged before max_iter
 *   - iterations: Number of iterations performed
 */
// [[Rcpp::export(rng = false)]]
List sc_weight_fw_cpp(arma::mat Y,
                      const double zeta,
                      const bool intercept,
                      arma::vec lambda,
                      const double min_decrease,
                      const int max_iter) {
  const int N0 = Y.n_rows;
  const int T0 = Y.n_cols - 1;
  const double min_dec_sq = min_decrease * min_decrease;

  // Initialize lambda with uniform weights if not provided
  if (lambda.n_elem == 0) {
    lambda = arma::vec(T0, arma::fill::ones) / static_cast<double>(T0);
  }

  // Remove column means if intercept requested
  if (intercept) {
    Y.each_row() -= arma::mean(Y, 0);
  }

  // Extract design matrix A (first T0 columns) and target vector b (last column)
  arma::mat A = Y.cols(0, T0 - 1);
  arma::vec b = Y.col(T0);

  // Scaled regularization parameter
  const double eta = N0 * (zeta * zeta);

  // Allocate storage for objective values
  arma::vec vals(max_iter);
  vals.fill(arma::datum::nan);

  arma::vec lambda_work = lambda;

  int final_iter = 0;
  bool converged = false;

  // Frank-Wolfe main loop — calls C++ directly, no R round-trip
  for (int t = 0; t < max_iter; ++t) {
    lambda_work = fw_step_internal(A, lambda_work, b, eta);

    // Compute objective function value
    arma::vec err = A * lambda_work - b;
    double sum_err_sq = arma::dot(err, err);
    double sum_lambda_sq = arma::dot(lambda_work, lambda_work);
    vals(t) = (zeta * zeta) * sum_lambda_sq + sum_err_sq / static_cast<double>(N0);

    final_iter = t + 1;

    // Early stopping for simplex corner convergence
    double max_weight = arma::max(lambda_work);
    if (max_weight > 0.99) {
      converged = true;
      break;
    }

    // Check convergence: stop if objective decrease is below threshold
    if (t >= 1) {
      double prev = vals(t - 1);
      double curr = vals(t);
      if (!std::isnan(prev) && !std::isnan(curr)) {
        if (!((prev - curr) > min_dec_sq)) {
          converged = true;
          break;
        }
      }
    }
  }

  return List::create(
    _["lambda"] = wrap(lambda_work),
    _["vals"] = wrap(vals),
    _["converged"] = converged,
    _["iterations"] = final_iter
  );
}
