#' Adaptive bootstrap inference for mediation effects
#'
#' Performs adaptive bootstrap inference for mediation effects based on
#' semiparametric estimates obtained from \code{SPsingle()}.
#'
#' The method implements an adaptive bootstrap procedure designed to improve
#' finite-sample inference when mediation effects are weak or partially
#' unidentified.
#'
#' @param sp A \code{data.table} returned by \code{SPsingle()}, containing
#'   parameter estimates, standard errors, and covariance components.
#' @param N.boot Number of bootstrap replicates. Default is \code{1e6}.
#' @param lambda Tuning parameter controlling the adaptive truncation region.
#'   Default is \code{1}.
#' @param lambda.n Optional numeric value specifying the truncation threshold.
#'   If \code{NULL}, it is computed as
#'   \eqn{\lambda \sqrt{n} / \log(n)}.
#' @param test Character string specifying the test statistic.
#'   Options include \code{"poc"} (product-of-coefficients)
#'   and \code{"maxp"} (minimum absolute t-statistic).
#'
#' @return
#' A named numeric vector containing:
#' \itemize{
#'   \item \code{AB.pval}: adaptive bootstrap p-value,
#'   \item \code{B.pval}: standard bootstrap p-value,
#'   \item \code{rej.prop}: proportion of truncated bootstrap samples,
#'   \item \code{lambda.n}: truncation threshold used.
#' }
#'
#' @details
#' The adaptive bootstrap modifies the standard parametric bootstrap by
#' truncating extreme realizations of weakly identified parameters. When
#' the parameter estimates fall within a neighborhood of zero, inference
#' is adjusted using an alternative resampling scheme to stabilize test
#' statistics.
#'
#' Two test statistics are supported:
#' \itemize{
#'   \item \code{"poc"}: product-of-coefficients statistic,
#'   \item \code{"maxp"}: minimum absolute t-statistic.
#' }
#'
#' @references
#' He, Y., Song, P. X.-K., and Xu, G. (2024).
#' Adaptive bootstrap tests for composite null hypotheses in the mediation
#' pathway analysis.
#' \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)},
#' 86(2), 411--434.
#'
#' @export

ABfast = function(sp, N.boot=1e6, lambda=1, lambda.n=NULL, test='poc') {
  index = match(c("thetaM", "betaA"), sp$variable)

  estimate = sp$estimate[index]
  se = sp$SE[index]
  covariance = as.matrix(sp[index, c("thetaM", "betaA")])
  covariance[2,1] = covariance[1,2]
  n = round(sp$thetaM[index[1]] / (sp$SE[index[1]])^2)

  t = estimate/se

  if (is.null(lambda.n)) {
    lambda.n = lambda * sqrt(n) / log(n)
  }

  ## Parametric bootstrap (sample from bivariate normal)
  Zstar = mvtnorm::rmvnorm(N.boot, mean=c(0,0), sigma=covariance)
  Estar = sweep(Zstar / sqrt(n), 2, estimate, '+') # sampled estimates
  Tstar = sweep(Estar, 2, se, '/') # corresponding t-statistics
  Istar = all(abs(t)<=lambda.n)  * (abs(Tstar[,1])<=lambda.n) * (abs(Tstar[,2])<=lambda.n) # corresponding indicator

  testFn = function(test) {
    if (tolower(test)=='poc') {
      obs = prod(estimate) # observed data
      boot = Estar[,1] * Estar[,2]
      Rstar = Zstar[,1] * Zstar[,2]
      Ustar = (boot - obs) * (1-Istar) + Rstar * Istar / n
    } else if (tolower(test)=='maxp') {
      obs = t[which.min(abs(t))] # observed data
      boot = apply(Tstar, 1, function(x) x[which.min(abs(x))])
      Kstar = sweep(Zstar / sqrt(n), 2, se, '/')
      Rstar = apply(Kstar, 1, function(x) x[which.min(abs(x))])
      Ustar = (boot - obs) * (1-Istar) + Rstar * Istar
    }

    AB.p = mean(Ustar < obs) # adaptive bootstrap
    B.p = mean((boot-obs) < obs) # nonparametric bootstrap
    c(AB.pval=min(AB.p,1-AB.p)*2, B.pval=min(B.p,1-B.p)*2, rej.prop=mean(Istar), lambda.n=lambda.n)
  }
  return(testFn(test))
}
