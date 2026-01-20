#' Semiparametric causal mediation analysis for a single mediator
#'
#' Fits the proposed semiparametric mediation model for a single mediator
#' using estimating equations and adaptive bootstrap inference.
#'
#' The method estimates the treatment effect on the outcome, the effect of
#' the mediator, and the natural indirect effect (NIE).
#'
#' @param dataset A matrix or data frame containing treatment, mediator,
#'   outcome, and covariates.
#' @param treatment Character string specifying the column name of the
#'   treatment variable (\eqn{A}).
#' @param mediator Character string specifying the column name of the
#'   mediator (\eqn{M}).
#' @param outcome Character string specifying the column name of the
#'   outcome (\eqn{Y}).
#' @param covariates Optional character vector of column names corresponding
#'   to baseline confounders (\eqn{C}). Default is \code{NULL}.
#' @param pi.link Character string specifying the link function for the
#'   treatment model. Default is \code{"logit"}.
#' @param rho.link Character string specifying the link function for the
#'   residual correlation model. Default is \code{"identity"}.
#' @param intercept Logical indicator for including an intercept term in
#'   the confounder matrix. Default is \code{TRUE}.
#' @param fast Logical indicator. If \code{TRUE}, a one-step estimator is
#'   used; otherwise, the estimating equations are solved directly using
#'   nonlinear equation solvers.
#' @param detail Logical indicator. If \code{TRUE}, returns additional
#'   variance components.
#'
#' @return
#' A \code{data.table} containing point estimates and standard errors for:
#' \itemize{
#'   \item \code{thetaA}: effect of treatment on outcome,
#'   \item \code{thetaM}: effect of mediator on outcome,
#'   \item \code{betaA}: effect of treatment on mediator,
#'   \item \code{NIE}: natural indirect effect.
#'   \item \code{SE}: standard error using delta method
#'   \item \code{t}: t-statistic using Sobel test
#'   \item \code{sobel.pval}: p-value using Sobel test
#' }
#'
#' @details
#' The estimator is based on semiparametric estimating equations that jointly
#' model the treatment, mediator, and outcome processes. Variance estimation
#' is obtained using the sandwich formula derived from the estimating
#' equations.
#'
#' The model is fitted one mediator at a time.
#'
#' @references
#' Sun, B. and Ye, T. (2023).
#' Semiparametric causal mediation analysis with unmeasured mediator--outcome
#' confounding.
#' \emph{Statistica Sinica}, 33, 2593--2612.
#'
#' @export

SPsingle = function(dataset, treatment, mediator, outcome, covariates = NULL,
                    pi.link='logit', rho.link='identity', intercept=TRUE, fast=FALSE, detail=FALSE) {
  #### Setup ####
  dataset = as.matrix(dataset)[,c(treatment,mediator,outcome,covariates)]
  dataset = dataset[!apply(dataset, 1, function(i) any(is.na(i)|i==Inf|i==-Inf)),] # omit rows with NA, Inf, or -Inf

  A = dataset[,treatment] # vector of A
  M = dataset[,mediator,drop=FALSE] # matrix of M
  Y = dataset[,outcome] # vector of Y
  C = cbind(intercept=1, dataset[,covariates,drop=FALSE]) # matrix of C

  #### Compute estimates ####
  fitAC = glm(A ~ C-1, family=ifelse(pi.link=='logit','binomial','gaussian'))
  fitYM = lm(Y ~ A + M + C - 1)
  fitMA = lm(M ~ A + C -1 )
  fitRho = lm(fitYM$residuals * fitMA$residuals ~ C - 1)

  init = c(fitAC$coef, fitYM$coef, fitMA$coef, fitRho$coef) # preliminary estimate

  if (fast) { # one-step estimator
    est = c(init - solve(SPgradient(init, Y, M, A, C, pi.link, rho.link)) %*% SPmoment(init, Y, M, A, C, pi.link, rho.link))
  } else { # solve directly
    est = nleqslv::nleqslv(init, SPmoment, Y=Y, M=M, A=A, C=C, pi.link=pi.link, rho.link=rho.link)$x
  }

  #### Compute asymptotic variances ####
  vM = SPmoment(est, Y, M, A, C, pi.link, rho.link, ind=TRUE)
  dM = SPgradient(est, Y, M, A, C, pi.link, rho.link)

  n = length(A); nC = ncol(C)
  thetaInd = c(nC+1:2, 2*nC+3) # indices for thetaA, thetaM, betaA

  mat = dM[1:3, -thetaInd] %*% solve(dM[-c(1:3),-thetaInd])
  matA = solve(dM[1:3,thetaInd] - mat %*% dM[-c(1:3),thetaInd])
  matB = t(vM[,1:3]) - mat %*% t(vM[,-c(1:3)])
  Sigma = matA %*% matB %*% t(matB) %*% t(matA) / n

  #dMinv = solve(dM)
  #Sigma = dMinv %*% (t(vM) %*% vM / n) %*% t(dMinv)

  theta = est[thetaInd] # esitmates for thetaA, thetaM, betaA
  se = sqrt(c(diag(Sigma), t(theta[3:2]) %*% Sigma[-1,-1] %*% theta[3:2]) / n) # SE for thetaA, thetaM, betaA

  output = data.table(variable = c('thetaA','thetaM','betaA','NIE'),
                      estimate = c(theta, theta[2]*theta[3]),
                      SE = se)
  #output[,'t'] = output[,estimate] / output[,SE]
  #output[,'sobel.pval'] = pnorm(-abs(output[,t])) * 2

  if (detail) {
    colnames(Sigma) = rownames(Sigma) = c('thetaA','thetaM','betaA')
    output = cbind(output, rbind(Sigma,NA))
  }
  return(output)
}
