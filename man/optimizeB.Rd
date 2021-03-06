% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{optimizeB}
\alias{optimizeB}
\title{Maximum likelihood optimization}
\usage{
optimizeB(B, S, C = diag(nrow(B)), ...)
}
\arguments{
\item{B}{initial coefficient matrix}

\item{S}{empirical covariance matrix}

\item{C}{positive definite term of the Lyapunov equation}

\item{...}{additional parameters passed to \code{optim}}
}
\value{
The MLE estimation of the coefficient matrix
}
\description{
Maximum likelihood optimization
}
\details{
MLE for coefficient matrix with the same 0s pattern of \code{B}
}
