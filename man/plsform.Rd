\name{plsform}
\alias{plsform}
\title{Create the structure of a linear mixed model from formula/data specification}
\usage{
  plsform(formula, data, REML = TRUE, weights, offset,
    sparseX = FALSE, ...)
}
\arguments{
  \item{formula}{a two-sided model formula with
  random-effects terms and, optionally, fixed-effects
  terms.}

  \item{data}{a data frame in which to evaluate the
  variables from \code{form}}

  \item{REML}{calculate REML deviance?}

  \item{weights}{prior weights}

  \item{offset}{offset}

  \item{sparseX}{should X, the model matrix for the
  fixed-effects coefficients be sparse?}

  \item{...}{additional arguments}
}
\value{
  a \code{list} with: \itemize{ \item \code{X} Fixed
  effects design (model) matrix \item \code{y} Observed
  response vector \item \code{fr} Model frame \item
  \code{call} Matched call \item \code{REML} Logical
  indicating REML or not \item \code{Zt} Transposed random
  effects design matrix \item \code{Lambdat} Transposed
  relative covariance factor \item \code{theta} Vector of
  covariance parameters \item \code{lower} Vector of lower
  bounds for \code{theta} \item \code{upper} Vector of
  upper bounds for \code{theta} \item \code{thfun} A
  function that maps \code{theta} into the structural
  non-zero elements of \code{Lambdat}, which are stored in
  \code{slot(Lambdat, 'x')} }
}
\description{
  A pure \code{R} implementation of the penalized least
  squares (PLS) approach to evaluation of the deviance or
  the REML criterion for linear mixed-effects models.
}
\examples{
form <- Reaction ~ Days + (Days|Subject)
data(sleepstudy, package="lme4")
ll <- plsform(form, sleepstudy, REML=FALSE)
names(ll)
}
\keyword{models}

