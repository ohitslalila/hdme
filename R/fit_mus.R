#' Matrix Uncertainty Selector
#' @description Matrix Uncertainty Selector
#' @param W Design matrix, measured with error. Must be a numeric matrix.
#' @param y Vector of responses.
#' @param lambda Regularization parameter.
#' @param delta Additional regularization parameter, bounding the measurement error.
#' @return Intercept and coefficients at the values of lambda and delta specified.
#' @references \insertRef{rosenbaum2010}{hdme}
#'
#' \insertRef{sorensen2018}{hdme}
#' @examples
#' # Example with logistic regression
#' set.seed(1)
#' # Number of samples
#' n <- 100
#' # Number of covariates
#' p <- 50
#' # True (latent) variables
#' X <- matrix(rnorm(n * p), nrow = n)
#' # Measurement matrix (this is the one we observe)
#' W <- X + matrix(rnorm(n*p, sd = 1), nrow = n, ncol = p)
#' # Coefficient vector
#' beta <- c(seq(from = 0.1, to = 1, length.out = 5), rep(0, p-5))
#' # Response
#' y <- X %*% beta + rnorm(n, sd = 1)
#' # Run the MU Selector
#' mus1 <- fit_mus(W, y)
#' # Draw an elbow plot to select delta
#' plot(mus1)
#'
#' # Now, according to the "elbow rule", choose the final delta where the curve has an "elbow".
#' # In this case, the elbow is at about delta = 0.08, so we use this to compute the final estimate:
#' mus2 <- fit_mus(W, y, delta = 0.08)
#' plot(mus2) # Plot the coefficients
#'
#' @export
fit_mus <- function(W, y, lambda = NULL, delta = NULL) {
  fit <- fit_gmus(W, y, lambda = lambda, delta = delta, family = "gaussian")
  return(fit)
}



