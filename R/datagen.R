#'
#' @title Data generation function
#' @description This function generates the simulation data (exposure X, outcome Y, and potential mediators M).
#' @description For each subject i = 1, ..., n:
#' - \eqn{X_i ~ i.i.d Bernoulli(p = 0.5)}
#' - \eqn{M_{i,k} = \alpha_k X_i + e_{1_{i,k}}}, where \eqn{e_{1_{i,k}} ~ i.i.d N(0, 1), k = 1, ..., K}
#' - \eqn{Y_i = X_i + \sum_k b_{i,k} M_{i,k} + e_{2_{i,k}}}, where \eqn{e_{2_{i,k}} ~ i.i.d. N(0, 1)}
#' @description The 6th-10th M variables (M6,M7,M8,M9,M10) are set to be the true mediators (i.e., having non-zero alpha and beta coefficients).
#' @description Specifically, the first 10 coefficients of alpha are set to take some non-zero value a, the 6th to 15th coefficients of beta are set to take some non-zero value b, and all other alpha and beta coefficients are set to be 0.
#'
#' @param n sample size
#' @param K number of potential mediators
#' @param a effect size for X-M path (i.e., non-zero alpha coefficient)
#' @param b effect size for M-Y path (i.e., non-zero beta coefficient)
#'
#' @return The resulting dataset has 3 lists: X, Y and M
#' @export
#'
#' @examples
#' datagen(n = 100, K = 50, a = 0.5, b = 0.3)
datagen <- function(n, K, a, b){
  alpha = rep(0,K); alpha[1:10] = a
  beta = rep(0,K); beta[6:15] = b
  gamma = 1

  X = rbinom(n,1,0.5)
  M = X %o% alpha + mvrnorm(n, mu = rep(0,K), Sigma = diag(rep(1,K)))
  Y = X*gamma + M %*% beta + rnorm(n,1)
  colnames(M) = paste("M", as.character(1:K), sep = "")

  list(X = X, Y = Y, M = M)
}
