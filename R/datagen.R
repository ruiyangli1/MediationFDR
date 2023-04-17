#'
#' @title Data generation function
#' @description This function generates the simulation data (exposure X, outcome Y, and potential mediators M).
#' @description For each subject i = 1, ..., n:
#' - \eqn{X_i \stackrel{i.i.d}{\sim} Bernoulli(prob = 0.5)}
#' - \eqn{M_{i,k} = \alpha_k X_i + e_{1_{i,k}}}, where \eqn{e_{1_{i,k}} \stackrel{i.i.d}{\sim} N(0, \Sigma), k = 1, ..., p, \Sigma} is the variance-covariance matrix with diagonal entries 1 and off-diagonal entries rho
#' - \eqn{Y_i = X_i + \sum_k b_{k} M_{i,k} + e_{2_{i}}}, where \eqn{e_{2_{i}} \stackrel{i.i.d}{\sim} N(0, 1)}
#' @description The 6th-15th M variables (M6,M7,M8,M9,M10,M11,M12,M13,M14,M15) are set to be the true mediators (i.e., having non-zero alpha and beta coefficients).
#' @description Specifically, the first 15 coefficients of alpha are set to take some non-zero value a, the 6th to 20th coefficients of beta are set to take some non-zero value b, and all other alpha and beta coefficients are set to be 0.
#'
#' @param n sample size
#' @param p number of potential mediators
#' @param a effect size for X-M path (i.e., non-zero alpha coefficient)
#' @param b effect size for M-Y path (i.e., non-zero beta coefficient)
#' @param rho effect size for the correlation structure of M
#'
#' @return The resulting data has 3 lists: X, Y and M
#' @export
#'
#' @examples
#' set.seed(20205); datagen(n = 1000, p = 100, a = 0.5, b = 0.3, rho = 0.1)
datagen <- function(n, p, a, b, rho){
  alpha = rep(0,p); alpha[1:15] = a
  beta = rep(0,p); beta[6:20] = b
  gamma = 1
  covmat = matrix(rho, nrow = p, ncol = p); diag(covmat) = 1

  X = rbinom(n,1,0.5)
  M = X %o% alpha + mvtnorm::rmvnorm(n, mean = rep(0,p), sigma = covmat) ## sigma: covariance matrix
  Y = X*gamma + M %*% beta + rnorm(n, mean = 0, sd = 1)
  colnames(M) = paste("M", as.character(1:p), sep = "")

  list(X = X, Y = Y, M = M)
}
