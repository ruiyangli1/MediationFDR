#'
#' @title Data generation function
#' @description This function generates the simulation data (exposure X, outcome Y, and potential mediators M).
#'
#' @param n sample size
#' @param K number of potential mediators
#' @param a effect size for X-M path
#' @param b effect size for M-Y path
#'
#' @return The resulting dataset has 3 lists: X, Y and M
#' @export
#'
#' @examples
#' datagen(n = 100, K = 50, a = 0.5, b = 0.3)
datagen<-function(n,K,a,b){
  alpha = rep(0,K); alpha[1:10] = a
  beta = rep(0,K); beta[6:15] = b
  gamma = 1

  X=rbinom(n,1,0.5)
  M=X%o%alpha+mvrnorm(n,mu=rep(0,K),Sigma=diag(rep(1,K)))
  Y=X*gamma+M%*%beta+rnorm(n,1)
  colnames(M) = paste("M",as.character(1:K),sep="")

  list(X=X,Y=Y,M=M)
}
