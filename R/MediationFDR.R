#'
#' @title Mediator selection function
#'
#' @description This function selects the mediators while controlling for FDR.
#'
#' @param X exposure of interest
#' @param Y outcome of interest
#' @param M potential mediators
#' @param C1 confounding variables in X-M path
#' @param C2 confounding variables in M-Y path
#' @param q1 FWER level in step 1
#' @param q2 FDR controlled level in step 2
#' @param f_func function form of W in step 2
#' @param correction_method multiple comparison correction method in step 1
#' @param weighted whether we weight the Z statistics in step 2
#' @param binaryOutcome whether the outcome is binary
#'
#' @return selected mediators
#' @export
#'
#' @examples
#' data = datagen(n = 200, K = 50, a = 0.8, b = 0.5)
#' X = data$X; Y = data$Y; M = data$M
#' MediationFDR(X,Y,M)
#'
MediationFDR <- function(X, Y, M,
                         C1 = NULL, C2 = NULL,
                         q1 = 0.05, q2 = 0.05,
                         f_func = "Max", correction_method = "BH",
                         weighted = FALSE,
                         binaryOutcome = FALSE){
  M = scale(M)
  X = scale(X)
  p = ncol(M)

  ID_SIS <- 1:p
  d <- p

  a_z <- matrix(0,1,d)
  a_p <- matrix(0,1,d)
  XC1 <- as.matrix(cbind(X,C1))

  # Step 1
  for (i in 1:d) {
    fit_a  <- lm(M[,ID_SIS[i]]~XC1)
    a_z[i] <- summary(fit_a)$coef[2,1]
    a_p[i] <- anova(fit_a)[5][1,1]
  }
  a_p <- p.adjust(a_p, method = correction_method)
  a_z[which(a_p > q1)] <- 0


  # Step 2
  X1 = as.matrix(cbind(M[,ID_SIS], X, C2))
  X1s = X1

  if (nrow(X1s) >= (2*ncol(X1s)) & !binaryOutcome) {
    Scalemat = diag(1/sqrt(colSums(X1s^2)))
    X1c = X1s %*% Scalemat
    X1Knock = create.fixed(X1c)$Xk
    X1Knock = X1Knock %*% solve(Scalemat)
  }else{
    X1Knock = create.second_order(X1s)
  }
  MKnock = X1Knock[,1:d]
  XC2 = as.matrix(cbind(X,C2))

  MMat = as.matrix(cbind(M[,ID_SIS], MKnock, XC2))
  if (binaryOutcome) {
    fit1 = cv.glmnet(MMat, c(Y), family = "binomial", standardize = FALSE, intercept = TRUE)
  }else{
    fit1 = cv.glmnet(MMat, c(scale(Y)), family = "gaussian", standardize = FALSE, intercept = FALSE)
  }

  Z1 = abs(coef(fit1)[1 + (1:d)])
  Z1tilde = abs(coef(fit1)[1 + d + (1:d)])

  if (weighted) {
    Z1 = Z1*abs(a_z)
    Z1tilde = Z1tilde*abs(a_z)
  } else {
    Z1[which(a_z == 0)] <- 0
    Z1tilde[which(a_z == 0)] <- 0
  }


  if (f_func == "Max") {
    W1 = pmax(Z1, Z1tilde)*(-1)^(Z1 <= Z1tilde)
  }
  if (f_func == "Diff") {
    W1 = Z1 - Z1tilde
  }
  mythred0 = knockoff.threshold(W1, fdr = q2, offset = 0)
  myselect0 = ID_SIS[which(W1 >= mythred0)]

  list(med_select = myselect0)
}
