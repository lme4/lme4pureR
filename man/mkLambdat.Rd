\name{mkLambdat}
\alias{mkLambdat}
\title{Make transposed random effects design matrix and relative covariance factor}
\usage{
  mkLambdat(grps, mms)
}
\arguments{
  \item{grps}{List of factor vectors of length n indicating
  groups.}

  \item{mms}{List of model matrices.}
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
  Create Lambdat as a ddiMatrix or a dgCMatrix and Zt as a
  dgCMatrix.
}
\details{
  Each list element in \code{grps} and \code{mms}
  corresponds to a random effects term.
}

