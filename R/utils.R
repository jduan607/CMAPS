#' Inverse link function
#'
#' Applies the inverse of a specified link function to a numeric input.
#'
#' This function is used internally to transform linear predictors back to
#' the response scale for different generalized linear model link functions.
#'
#' @param x A numeric vector or matrix of linear predictors.
#' @param link A character string specifying the link function. Supported
#'   options include \code{"logit"}, \code{"identity"}, \code{"log"},
#'   \code{"inverse"}, \code{"neginverse"}, \code{"sqroot"},
#'   \code{"invsq"}, \code{"loglog"}, and \code{"cloglog"}.
#'
#' @return A numeric vector or matrix on the response scale.
#'
#' @keywords internal
#' @noRd

inv.link = function(x, link) {
  switch(link, logit = 1/(1 + exp(-x)), identity = x, log = exp(x),
         inverse = x^(-1), neginverse = -(x^(-1)), sqroot = x^2,
         invsq = x^(-1/2), loglog = exp(-exp(-x)), cloglog = 1 - exp(-exp(x)))
}

#' Jacobian and Hessian matrix for link functions
#'
#' Computes the first- or second-order derivative of the mean function with
#' respect to the parameter vector under a specified link function.
#'
#' This function is used internally for calculating gradients and Hessians
#' in the semiparametric estimation procedure.
#'
#' @param x A numeric vector of length \eqn{p} containing model parameters.
#' @param constant An \eqn{n \times p} design matrix, where \eqn{n} is the
#'   number of observations.
#' @param link A character string specifying the link function. Supported
#'   options include \code{"logit"}, \code{"identity"}, and \code{"log"}.
#' @param w Optional numeric vector of observation weights. Defaults to 1.
#' @param order Integer indicating the derivative order. Use \code{1} for
#'   the Jacobian (first derivative) and \code{2} for the Hessian
#'   (second derivative).
#'
#' @return
#' If \code{order = 1}, returns an \eqn{n \times p} Jacobian matrix.
#' If \code{order = 2}, returns a \eqn{p \times p} Hessian matrix.
#'
#' @keywords internal
#' @noRd
jacobian <- function(x, constant, link, w = 1, order = 1) {
  if (order==1) {
    switch(link,
           logit = c(exp(constant %*% x)  / (exp(constant %*% x) + 1)^2) * constant,
           identity = constant,
           log = c(exp(constant %*% x)) * constant )
  } else if (order==2) {
    switch(link,
           identity = matrix(0, nrow=length(x), ncol=length(x)),
           logit = t(c(-(exp(constant %*% x)-1) * exp(constant %*% x) / (exp(constant %*% x)+1)^3) * constant * w) %*% constant,
           log = t(c(exp(constant %*% x)) * constant * w) %*% constant )
  }
}

#' Semiparametric estimating functions
#'
#' Computes the estimating functions used in the semiparametric mediation
#' model. The function evaluates moment conditions for treatment,
#' mediator, and outcome models, as well as their associated nuisance
#' parameters.
#'
#' The estimating equations are evaluated either at the individual level
#' or as sample averages, depending on the value of \code{ind}.
#'
#' @param par Numeric vector of model parameters ordered as
#'   \eqn{(\eta_1, \theta_A, \theta_M, \eta_2 / \theta_C,
#'   \beta_A, \eta_3 / \beta_C, \eta_4)}.
#' @param Y Numeric vector of outcomes.
#' @param M Numeric vector of mediators.
#' @param A Numeric vector of treatment indicators.
#' @param C Numeric matrix of confounders.
#' @param pi.link Character string specifying the link function for the
#'   treatment model.
#' @param rho.link Character string specifying the link function for the
#'   residual correlation model.
#' @param ind Logical indicator. If \code{TRUE}, returns individual-level
#'   estimating functions; if \code{FALSE}, returns their sample averages.
#'
#' @return
#' If \code{ind = TRUE}, returns an \eqn{n \times q} matrix of estimating
#' functions evaluated at the individual level, where \eqn{n} is the number
#' of observations and \eqn{q} is the total number of moment conditions.
#'
#' If \code{ind = FALSE}, returns a numeric vector of length \eqn{q}
#' corresponding to the averaged estimating equations.
#'
#' @details
#' The estimating functions consist of two components:
#' \itemize{
#'   \item \eqn{\phi}: moment conditions for the treatment, mediator,
#'   and outcome models.
#'   \item \eqn{\gamma}: score-type equations for nuisance parameters,
#'   including the treatment and residual correlation models.
#' }
#'
#' @keywords internal
#' @noRd
SPmoment <- function(par, Y, M, A, C, pi.link, rho.link, ind = FALSE) {
  n = length(A); nC = ncol(C)
  eta1 = par[1:nC]; theta = par[nC+1:(nC+2)]; beta = par[2*nC+2+1:(nC+1)]; eta4 = par[3*nC+3+1:nC]

  RA = c(A - inv.link(C %*% eta1, pi.link))
  RY = c(Y - cbind(A,M,C) %*% theta)
  RM = c(M - cbind(A,C) %*% beta)
  RYM = c(RY * RM - inv.link(C %*% eta4, rho.link))
  pi.grad = jacobian(eta1, C, pi.link)
  rho.grad = jacobian(eta4, C, rho.link)

  if (ind) {
    phi.val = RA * cbind(RY, RM, RYM)
    gamma.val = cbind(pi.grad * RA, C * RY, C * RM, rho.grad * RYM)
    return(cbind(phi.val, gamma.val))
  } else {
    phi.val = t(RA) %*% cbind(RY, RM, RYM)
    gamma.val = c(t(pi.grad) %*% RA, t(C) %*% RY, t(C) %*% RM, t(rho.grad) %*% RYM)
    return(c(phi.val, gamma.val)/n)
  }
}

#' Gradient matrix of semiparametric estimating equations
#'
#' Computes the Jacobian (gradient) matrix of the estimating equations
#' used in the semiparametric mediation model with respect to all model
#' parameters.
#'
#' The gradient matrix is required for variance estimation and for the
#' adaptive bootstrap procedure.
#'
#' @param par Numeric vector of model parameters ordered as
#'   \eqn{(\eta_1, \theta_A, \theta_M, \eta_2 / \theta_C,
#'   \beta_A, \eta_3 / \beta_C, \eta_4)}.
#' @param Y Numeric vector of outcomes.
#' @param M Numeric vector of mediators.
#' @param A Numeric vector of treatment indicators.
#' @param C Numeric matrix of confounders.
#' @param pi.link Character string specifying the link function for the
#'   treatment model.
#' @param rho.link Character string specifying the link function for the
#'   residual correlation model.
#'
#' @return
#' A square numeric matrix of dimension \eqn{p \times p}, where \eqn{p}
#' is the total number of parameters. Each entry corresponds to the
#' partial derivative of an estimating function with respect to a model
#' parameter.
#'
#' @details
#' The gradient matrix is constructed by analytically differentiating the
#' estimating equations with respect to all model parameters. The resulting
#' matrix is scaled by the sample size \eqn{n}.
#'
#' @keywords internal
#' @noRd
SPgradient <- function(par, Y, M, A, C, pi.link, rho.link) {
  n = length(A); nC = ncol(C)
  eta1 = par[1:nC]; theta = par[nC+1:(nC+2)]; beta = par[2*nC+2+1:(nC+1)]; eta4 = par[3*nC+3+1:nC]

  RA = c(A - inv.link(C %*% eta1, pi.link))
  RY = c(Y - cbind(A,M,C) %*% theta)
  RM = c(M - cbind(A,C) %*% beta)
  RYM = c(RY * RM - inv.link(C %*% eta4, rho.link))
  pi.grad = jacobian(eta1, C, pi.link)
  rho.grad = jacobian(eta4, C, rho.link)

  indices = list(eta1 = 1:nC, theta = nC+1:(nC+2), beta = 2*nC+2+1:(nC+1), eta4 = 3*nC+3+1:nC,
                 gamma1 = 3+1:nC, gamma2 = 3+nC+1:nC, gamma3 = 3+2*nC+1:nC, gamma4 = 3+3*nC+1:nC)
  J = matrix(0, nrow=length(par), ncol=length(par))

  J[1:3, indices$eta1] = -t(cbind(RY, RM, RYM)) %*% pi.grad
  J[c(1,3), indices$theta] = -t(cbind(RA, RA*RM)) %*% cbind(A,M,C)
  J[c(2,3), indices$beta] = -t(cbind(RA, RA*RY)) %*% cbind(A,C)
  J[3, indices$eta4] = -t(RA) %*% rho.grad

  J[indices$gamma1, indices$eta1] = jacobian(eta1, C, pi.link, w=RA, order=2) - t(pi.grad) %*% pi.grad
  J[c(indices$gamma2,indices$gamma4), indices$theta] = -t(cbind(C, rho.grad*RM)) %*% cbind(A,M,C)
  J[c(indices$gamma3,indices$gamma4), indices$beta] = -t(cbind(C, rho.grad*RY)) %*% cbind(A,C)
  J[indices$gamma4, indices$eta4] = jacobian(eta4, C, rho.link, w=RYM, order=2) - t(rho.grad) %*% rho.grad
  return(J / n)
}

