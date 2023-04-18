#'
#' @title Data generation function for counter example
#' @description This function generates the simulation data (exposure X, outcome Y, potential mediators M, and covariate V1) for counter example.
#' @description The number of potential mediators is set to be p = 19.
#' @description For each subject i = 1, ..., n:
#' - \eqn{X_i \stackrel{i.i.d}{\sim} N(0, 1)}
#' - \eqn{{V_1}_i \stackrel{i.i.d}{\sim} N(0, 1)}
#' - \eqn{cov(X, V_1) = 0.8}
#' - \eqn{M_{i,k} = \alpha_k X_i + {V_1}_i {\eta_1}_k + e_{1_{i,k}}}, where \eqn{e_{1_{i,k}} \stackrel{i.i.d}{\sim} N(0, 1), k = 1, ..., p}, and \eqn{cov(e_{1_{i,k}},e_{1_{i,j}} = 0)} for k not equal j
#' - \eqn{Y_i = X_i + \sum_k \beta_{k} M_{i,k} + e_{2_{i}}}, where \eqn{e_{2_{i}} \stackrel{i.i.d}{\sim} N(0, 1)}
#' @description We let \eqn{\alpha_1, \cdots, \alpha_{10}} be randomly from the uniform distribution from -1 to 1 and \eqn{\alpha_{11}, \cdots, \alpha_{19}} be 0.
#' @description We let \eqn{{\eta_1}_{11} = 1} and the rest of the \eqn{\eta_1} values be 0.
#' @description We let \eqn{\beta_1, \cdots, \beta{9} = 0} and \eqn{\beta_{10}, \cdots, \beta{19}} be either 1 or -1, randomly.
#' @description M10 is the true mediators (i.e., having non-zero alpha and beta coefficients).
#'
#' @param n sample size
#'
#' @return The resulting data has 4 lists: X, Y, M, V1
#' @export
#'
#' @examples
#' set.seed(23)
#' data = datagen.counter(n = 1000)
#' head(data$X)
#' head(data$Y)
#' head(data$M)
#' head(data$V1)
datagen.counter <- function(n){
  p = 19
  XV1 = mvtnorm::rmvnorm(n, mean = rep(0,2), sigma = matrix(c(1,.8,.8,1),2,2)) # sigma is covariance matrix
  X = as.matrix(XV1[,1])
  V1 = as.matrix(XV1[,2])
  alpha = rep(0,p); alpha[1:10] = runif(10, min = -1, max = 1)
  eta1 = rep(0,p); eta1[11] = 1
  e1 = matrix(rnorm(n*p, 0, 1), nrow = n, ncol = p)

  M = X %*% alpha + V1 %*% eta1 + e1
  colnames(M) = paste("M", as.character(1:p), sep = "")

  beta = rep(0,p); beta[10:19] = sample(c(-1,1),size=10,replace=TRUE)
  gamma = 1
  e2 = matrix(rnorm(n, 0, 1), nrow = n)

  Y = X %*% gamma + M %*% beta + e2

  list(X = X, Y = Y, M = M, V1 = V1, alpha = alpha, eta1 = eta1, beta = beta)
}
