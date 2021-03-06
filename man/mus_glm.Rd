% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mus_glm.R
\name{mus_glm}
\alias{mus_glm}
\title{Generalized Matrix Uncertainty Selector for logistic regression}
\usage{
mus_glm(W, y, lambda, delta, family = c("binomial", "poisson"))
}
\arguments{
\item{W}{Design matrix, measured with error.}

\item{y}{Vector of the binomial response value.}

\item{lambda}{Regularization parameter due to model error.}

\item{delta}{Regularization parameter due to measurement error.}

\item{family}{"binomial" or "poisson"}
}
\value{
Intercept and coefficients at the values of lambda and delta specified.
}
\description{
Internal function.
}
