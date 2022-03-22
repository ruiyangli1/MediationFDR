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
#' datagen(100, 50, 0.5, 0.3)
datagen<-function(n,K,a,b){
  alpha = rep(0,K); alpha[1:10] = a
  gamma = 1
  beta = rep(0,K); beta[6:15] = b

  X=rbinom(n,1,0.5)
  M=X%o%alpha+mvrnorm(n,mu=rep(0,K),Sigma=diag(rep(1,K)))
  Y=X*gamma+M%*%beta+rnorm(n,1)
  #Z=NULL
  data=data.frame(cbind(X,Y,M,Z))
  names(data)=c("X","Y",paste("M",as.character(1:K),sep=""))
  list(M=M,X=X,Y=Y#,Z=Z
       )
}
