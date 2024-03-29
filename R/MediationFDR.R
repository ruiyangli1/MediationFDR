#'
#' @title Mediator selection function
#'
#' @description This function selects the mediators while controlling for FDR.
#'
#' @param X exposure of interest
#' @param Y outcome of interest
#' @param M potential mediators
#' @param V1 confounding variables in X-M path
#' @param V2 confounding variables in M-Y path
#' @param q1 FWER level in step 1
#' @param q2 FDR controlled level in step 4
#' @param f_func function form of W in step 4
#' @param correction_method multiple comparison correction method in step 1
#' @param weighted whether we weight the Z statistics in step 3
#' @param binaryOutcome whether the outcome is binary
#'
#' @return selected mediators
#' @export
#'
#' @examples
#' set.seed(20205)
#' data = datagen(n = 1000, p = 100, a = 0.5, b = 0.3, rho = 0.1)
#' MediationFDR(X = data$X, Y = data$Y, M = data$M, V1 = NULL, V2 = NULL, q1 = 0.025, q2 = 0.025, f_func = "Diff", correction_method = "holm", weighted = FALSE)
MediationFDR <- function(X, Y, M,
                         V1 = NULL, V2 = NULL,
                         q1 = 0.025, q2 = 0.025,
                         f_func = "Diff", correction_method = "holm",
                         weighted = FALSE,
                         binaryOutcome = FALSE){

  M = scale(M)
  X = scale(X)
  if (!is.null(V1)) {
    V1 = scale(V1)
  }
  if (!is.null(V2)) {
    V2 = scale(V2)
  }

  p = ncol(M)

  ID_SIS <- 1:p
  d <- p

  a_z <- matrix(0,1,d)
  a_p <- matrix(0,1,d)
  XV1 <- as.matrix(cbind(X,V1))

  # Step 1
  for (i in 1:d) {
    fit_a  <- lm(M[,ID_SIS[i]]~XV1)
    a_z[i] <- summary(fit_a)$coef[2,1]
    a_p[i] <- summary(fit_a)$coef[2,4]
  }
  a_p <- p.adjust(a_p, method = correction_method)
  a_z[which(a_p > q1)] <- 0

  patha = which(!a_z == 0)

  # Step 2
  X1 = as.matrix(cbind(M[,ID_SIS], X, V2))
  X1s = X1

  if (nrow(X1s) >= (2*ncol(X1s)) & !binaryOutcome) {
    Scalemat = diag(1/sqrt(colSums(X1s^2)))
    X1c = X1s %*% Scalemat
    X1Knock = create.fixed(X1c)$Xk # create.fixed() needs normalized input
    X1Knock = X1Knock %*% solve(Scalemat)
  }else{
    X1Knock = create.second_order(X1s) # create.second_order() needs original input
  } # X1Knock: knockoff copy for M,X,V2

  MKnock = X1Knock[,1:d] # knockoff copy for M
  XV2 = as.matrix(cbind(X,V2))

  MMat = as.matrix(cbind(M[,ID_SIS], MKnock, XV2))
  penalty = rep(0, ncol(MMat))
  penalty[1:(2*d)] = 1

  if (binaryOutcome) {
    fit1 = cv.glmnet(MMat, c(Y), family = "binomial", standardize = FALSE, intercept = TRUE, penalty.factor = penalty)
  }else{
    fit1 = cv.glmnet(MMat, c(scale(Y)), family = "gaussian", standardize = FALSE, intercept = FALSE, penalty.factor = penalty)
  } ## not standardize MMat variables (notice: they are on the same scale, they all have been SCALED)
  ## want intercept = 0

  Z1 = abs(coef(fit1)[1 + (1:d)])
  Z1tilde = abs(coef(fit1)[1 + d + (1:d)])

  pathb = which(!Z1 == 0)
  pathb.tilde = which(!Z1tilde == 0)

  # Step 3
  if (weighted) {
    Z1 = Z1*abs(a_z)
    Z1tilde = Z1tilde*abs(a_z)
  } else {
    Z1[which(a_z == 0)] <- 0
    Z1tilde[which(a_z == 0)] <- 0
  }

  pathb.a = which(!Z1 == 0)
  pathb.tilde.a = which(!Z1tilde == 0)

  # Step 4
  if (f_func == "Max") {
    W1 = pmax(Z1, Z1tilde)*(-1)^(Z1 <= Z1tilde)
  }
  if (f_func == "Diff") {
    W1 = Z1 - Z1tilde
  }
  mythred0 = knockoff.threshold(W1, fdr = q2, offset = 0)
  myselect0 = ID_SIS[which(W1 >= mythred0)]

  list(med_select = myselect0,
       path.a = patha,
       path.b = pathb, path.b.tilde = pathb.tilde,
       path.b.weighta = pathb.a, path.b.tilde.weighta = pathb.tilde.a
  )
}
