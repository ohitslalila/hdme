#' Generalized Matrix Uncertainty Selector for logistic regression
#'
#' @description Internal function.
#' @param W Design matrix, measured with error.
#' @param y Vector of the binomial response value.
#' @param lambda Regularization parameter due to model error.
#' @param delta Regularization parameter due to measurement error.
#' @param family "binomial" or "poisson"
#' @return Intercept and coefficients at the values of lambda and delta specified.
#'
mus_glm <- function(W, y, lambda, delta, family = c("binomial", "poisson"), alternative = F){

  family <- match.arg(family)

  if(family == "binomial") {
    mu <- logit
    dmu <- dlogit
  } else if(family == "poisson") {
    mu <- pois
    dmu <- dpois
  }
    

  n <- dim(W)[1]
  p <- dim(W)[2]
  
  if(T){
  #if(alternative == F){
      W <- scale(W[,2:p])
      scales <- attr(W, "scaled:scale")
      W <- cbind(rep(1,n), W) 
  }


  bOld <- stats::rnorm(p)/p
  bNew <- stats::rnorm(p)/p
  IRLSeps <- 1e-7
  maxit <- 100
  count <- 1
  Diff1 <- 1
  Diff2 <- 1

  while(Diff1 > IRLSeps & Diff2 > IRLSeps & count < maxit){
      bOlder <- bOld
      bOld <- bNew
      V <- dmu(W%*%bOld)
      z <- W%*%bOld + (y - mu(W%*%bOld))/dmu(W%*%bOld)
      Wtilde <- c(sqrt(V)) * W
      ztilde <- c(sqrt(V)) * c(z)
    
      if(alternative == F){
      #if(alternative == F | (Diff1 > 10 ^ (-2) & count < 10)){
      bNew <- musalgorithm(Wtilde, ztilde, lambda, delta * sqrt(sum((V)^2)) / sqrt(n))
        }
      if(alternative == T){ 
    #if(alternative == T & (Diff1 <= 10 ^ (-2) | count >= 10)){
        bNew <- musalgorithm_alt(Wtilde, ztilde, lambda, delta * abs(t(V)) %*% abs(W[,-1]) / n)
       #bNew <- musalgorithm_alt(Wtilde, ztilde, lambda, delta * sqrt(sum((V)^2)) / sqrt(n))
        }

      count <- count+1
      print(count)
      Diff1 <- sum(abs(bNew - bOld))
      Diff2 <- sum(abs(bNew - bOlder))
  }
  if(count >= maxit) print(paste("Did not converge"))
  
  if(alternative == F){
      return(c(bNew[1], bNew[2:p] / scales))
    }
  if(alternative == T){
      #return(bNew)
      return(c(bNew[1], bNew[2:p] / scales))
  }
}
