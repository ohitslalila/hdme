else
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
mus_glm_nb <- function(W, y, lambda, delta, theta, family = "nb", alternative = F){

  family <- match.arg(family)

  if(family == "nb") {
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
    Utilde <- theta / (theta + V)
    Q <- 1 / (theta + V)
    delta1 <- delta * sum(bOld^2) / n * abs(t(Q) %*% W[, -1])
    bNew <- musalgorithm_nb(Wtilde, Utilde, ztilde, lambda, delta, delta1)
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
