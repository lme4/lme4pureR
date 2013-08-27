\name{plsform}
\alias{plsform}
\title{Create linear mixed model deviance function from formula/data specification}
\usage{
  plsform(formula, data, REML = TRUE, weights, offset, ...)
}
\arguments{
  \item{formula}{a two-sided model formula with
  random-effects terms and, optionally, fixed-effects
  terms.}

  \item{data}{a data frame in which to evaluate the
  variables from \code{form}}

  \item{weights}{prior weights}

  \item{offset}{offset}

  \item{REML}{calculate REML deviance?}
}
\value{
  a function that evaluates the deviance or REML criterion
}
\description{
  A pure \code{R} implementation of the penalized least
  squares (PLS) approach to evaluation of the deviance or
  the REML criterion for linear mixed-effects models.
}
\keyword{models}

