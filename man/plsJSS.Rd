\name{plsJSS}
\alias{plsJSS}
\title{Linear mixed model deviance function as it
appears in the pseudocode of the JSS article}
\usage{
plsJSS(X, y, Zt, Lambdat, mapping, weights, offset = numeric(n),
  REML = TRUE, ...)
}
\arguments{
  \item{X}{fixed effects model matrix}

  \item{y}{response}

  \item{Zt}{transpose of the sparse model matrix for the
  random effects}

  \item{Lambdat}{upper triangular sparse Cholesky factor of
  the relative covariance matrix of the random effects}

  \item{mapping}{a function that takes a value of
  \code{theta} and produces the non-zero elements of
  \code{Lambdat}.  The structure of \code{Lambdat} cannot
  change, only the numerical values}

  \item{weights}{prior weights}

  \item{offset}{offset}

  \item{REML}{calculate REML deviance?}

  \item{...}{additional arguments}
}
\value{
a function that evaluates the deviance or REML criterion
}
\description{
A pure \code{R} implementation of the penalized least
squares (PLS) approach for computing linear mixed model
deviances. The purpose is to clarify how PLS works without
having to read through C++ code, and as a sandbox for
trying out modifications to PLS.
}
\keyword{models}

