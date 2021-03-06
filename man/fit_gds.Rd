% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_gds.R
\name{fit_gds}
\alias{fit_gds}
\title{Generalized Dantzig Selector}
\usage{
fit_gds(X, y, lambda = NULL, family = c("gaussian", "binomial"))
}
\arguments{
\item{X}{Design matrix.}

\item{y}{Vector of the continuous response value.}

\item{lambda}{Regularization parameter. Only a single value is supported.}

\item{family}{Use "gaussian" for linear regression and "binomial" for logistic regression.}
}
\value{
Intercept and coefficients at the values of lambda specified.
}
\description{
Generalized Dantzig Selector
}
\examples{
# Example with logistic regression
n <- 1000  # Number of samples
p <- 10 # Number of covariates
X <- matrix(rnorm(n * p), nrow = n) # True (latent) variables # Design matrix
beta <- c(seq(from = 0.1, to = 1, length.out = 5), rep(0, p-5)) # True regression coefficients
y <- rbinom(n, 1, (1 + exp(-X \%*\% beta))^(-1)) # Binomially distributed response
gds <- fit_gds(X, y, family = "binomial")

}
\references{
\insertRef{candes2007}{hdme}

\insertRef{james2009}{hdme}
}
