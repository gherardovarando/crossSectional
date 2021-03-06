% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/triangular.R
\name{l1Bfast}
\alias{l1Bfast}
\title{Estimate triangular OU from cross sectional data}
\usage{
l1Bfast(S, D = diag(ncol(S)), alpha = 0.1, B = NULL, ...)
}
\arguments{
\item{S}{the sample covariance matrix}

\item{D}{noise matrix}

\item{alpha}{penalization coefficient}

\item{B}{initial estimate for the coefficient matrix}

\item{...}{additional parameters passed to \code{optim}}
}
\description{
Use l1
}
