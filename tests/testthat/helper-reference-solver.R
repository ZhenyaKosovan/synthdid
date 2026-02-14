#' Solve a ridge-regularized simplex-constrained least squares problem
#'
#' Minimizes \eqn{||Ax - b||^2 + zeta^2 n ||x||^2} subject to \eqn{x} lying in
#' the simplex. When \code{intercept = TRUE}, includes a free intercept term
#' \eqn{x0} in the objective.
#'
#' @param A Design matrix.
#' @param b Response vector.
#' @param zeta Non-negative ridge penalty multiplier.
#' @param intercept Logical flag indicating whether to include an intercept.
#'
#' @return Numeric vector of simplex-constrained coefficients.
#' @keywords internal
simplex.least.squares =  function(A, b, zeta = 0, intercept = FALSE) {
    x = CVXR::Variable(ncol(A))
    constraints = list(sum(x) == 1, x >= 0)
    if(intercept) {
	x0 = CVXR::Variable(1)
	objective = sum((A %*% x + x0 - b)^2) + zeta^2 * length(b) * sum(x^2)
    } else {
	objective = sum((A %*% x - b)^2) + zeta^2 * length(b) * sum(x^2)
    }
    cvx.problem = CVXR::Problem(CVXR::Minimize(objective), constraints)
    cvx.output = CVXR::solve(cvx.problem, solver = 'ECOS')
    as.numeric(cvx.output$getValue(x))
}


#' Default noise scale estimate
#'
#' Computes the standard deviation of first differences for pre-treatment
#' control units, used as a plug-in noise estimate.
#'
#' @inheritParams synthdid.reference
#'
#' @return Scalar noise level.
#' @keywords internal
sigma.default = function(Y, N0, T0) { sd(apply(Y[1:N0,1:T0], 1, diff)) }
epsilon = 1e-6

#' Reference implementation of the synthdid estimator
#'
#' Computes the synthetic difference-in-differences estimate using convex
#' optimization for the reference implementation described in the paper.
#'
#' @param Y Outcome matrix with controls first and pre-treatment periods first.
#' @param N0 Number of control units.
#' @param T0 Number of pre-treatment periods.
#' @param zeta.omega Ridge penalty for the unit weights \code{omega}.
#'
#' @return Scalar treatment effect estimate.
#' @keywords internal
synthdid.reference = function(Y, N0, T0, zeta.omega=((nrow(Y)-N0)*(ncol(Y)-T0))^(1/4) * sigma.default(Y, N0, T0)) {
    N = nrow(Y); T=ncol(Y); N1 = N-N0; T1=T-T0;
    lambda = simplex.least.squares(Y[1:N0, 1:T0],    rowMeans(Y[1:N0, (T0+1):T, drop=FALSE]), zeta=epsilon*sigma.default(Y,N0,T0), intercept=TRUE)
    omega =  simplex.least.squares(t(Y[1:N0, 1:T0]), colMeans(Y[(N0+1):N, 1:T0, drop=FALSE]), zeta=zeta.omega, intercept=TRUE)
    estimate = t(c(-omega, rep(1/N1, N1))) %*% Y  %*% c(-lambda, rep(1/T1, T1))
}

#' Reference synthetic control estimator
#'
#' Computes the SC analogue of \code{synthdid.reference} with an optional ridge
#' penalty on unit weights.
#'
#' @inheritParams synthdid.reference
#' @param zeta.omega Ridge penalty for the unit weights \code{omega}.
#'
#' @return Scalar synthetic control estimate.
#' @keywords internal
sc.reference = function(Y, N0, T0, zeta.omega=1e-6 * sigma.default(Y, N0, T0)) {
    N = nrow(Y); T=ncol(Y); N1 = N-N0; T1=T-T0;
    omega =  simplex.least.squares(t(Y[1:N0, 1:T0]), colMeans(Y[(N0+1):N, 1:T0, drop=FALSE]), zeta=zeta.omega, intercept=FALSE)
    estimate = t(c(-omega, rep(1 / N1, N1))) %*% Y  %*% c(-rep(0, T0), rep(1/T1, T1))
}

#' Reference difference-in-differences estimator
#'
#' Computes the canonical DiD estimate using uniform weights on controls and
#' pre-treatment periods.
#'
#' @inheritParams synthdid.reference
#'
#' @return Scalar DiD estimate.
#' @keywords internal
did.reference = function(Y, N0, T0) {
    N = nrow(Y); T=ncol(Y); N1 = N-N0; T1=T-T0;
    estimate = t(c(-rep(1/N0, N0), rep(1 / N1, N1))) %*% Y  %*% c(-rep(1/T0, T0), rep(1 / T1, T1))
}
