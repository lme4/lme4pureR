\name{Zsection}
\alias{blockLambdat}
\alias{mkRanefRepresentation}
\alias{Zsection}
\title{Create a section of a transposed random effects design matrix}
\usage{
  Zsection(grp, mm)

  blockLambdat(nl, nc)

  mkRanefRepresentation(grps, mms)
}
\arguments{
  \item{grp}{Grouping factor for a particular random
  effects term.}

  \item{mm}{Dense model matrix for a particular random
  effects term.}

  \item{nl}{Number of levels in a grouping factor for a
  particular random effects term (the number of levels in
  the \code{grp} argument in \code{\link{Zsection}}).}

  \item{nc}{Number of columns in a dense model matrix for a
  particular random effects term (the number of columns in
  the \code{mm} argument in \code{\link{Zsection}}).}

  \item{grps}{List of factor vectors of length n indicating
  groups.  Each element corresponds to a random effects
  term.}

  \item{mms}{List of model matrices.  Each element
  corresponds to a random effects term.}
}
\value{
  Dense section of a random effects design matrix

  A \code{list} with: \itemize{ \item the block \item
  ititial values of theta for the block \item lower bounds
  on these initial theta values \item a function that
  updates the block given the section of theta for this
  block }

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
  Create a section of a transposed random effects design
  matrix

  Each random-effects term is represented by diagonal block
  on the transposed relative covariance factor.
  \code{blockLambdat} creates such a block, and returns
  related information along with it.

  Create all of the elements required to specify the
  random-effects structure of a mixed effects model.
}
\details{
  The basic idea of this function is to call
  \code{\link{Zsection}} and \code{\link{blockLambdat}}
  once for each random effects term (ie. each list element
  in \code{grps} and \code{mms}). The results of
  \code{\link{Zsection}} for each term are \code{rBind}ed
  together. The results of \code{\link{blockLambdat}} are
  \code{bdiag}ed together, unless all terms have only a
  single column ('predictor') in which case a diagonal
  matrix is created directly.
}
\examples{
## consider a term (x | g) with:
## number of observations, n = 6
## number of levels, nl = 3
## number of columns ('predictors'), nc = 2
(X <- cbind("(Intercept)"=1,x=1:6)) # an intercept in the first column and 1:6 predictor in the other
(g <- gl(3,2,labels=letters[1:3])) # grouping factor
stopifnot(nrow(X) == length(g))
ncol(X) # nc = 2
nlevels(g) # nl = 3
Zsection(g, X)
(l <- blockLambdat(2, 3))
within(l, slot(Lambdat, 'x') <- updateLambdatx(as.numeric(10:12)))
}

