\name{mkRanefRepresentation}
\alias{mkRanefRepresentation}
\title{Make random effects representation}
\usage{
  mkRanefRepresentation(grps, mms)
}
\arguments{
  \item{grps}{List of factor vectors of length n indicating
  groups.  Each element corresponds to a random effects
  term.}

  \item{mms}{List of model matrices.  Each element
  corresponds to a random effects term.}
}
\value{
  A \code{list} with: \itemize{ \item \code{Lambdat}
  Transformed relative covariance factor \item \code{Zt}
  Transformed random effects design matrix \item
  \code{theta} Vector of covariance parameters \item
  \code{lower} Vector of lower bounds on the covariance
  parameters \item \code{upper} Vector of upper bounds on
  the covariance parameters \item \code{thfun} A function
  that maps \code{theta} onto the non-zero elements of
  \code{Lambdat} }
}
\description{
  Create all of the elements required to specify the
  random-effects structure of a mixed effects model.
}
\details{
  The basic idea of this function is to call
  \code{\link{Zsection}} and \code{\link{LambdatBlock}}
  once for each random effects term (ie. each list element
  in \code{grps} and \code{mms}). The results of
  \code{\link{Zsection}} for each term are \code{rBind}ed
  together. The results of \code{\link{LambdatBlock}} are
  \code{bdiag}ed together, unless all terms have only a
  single column ('predictor') in which case a diagonal
  matrix is created directly.
}

