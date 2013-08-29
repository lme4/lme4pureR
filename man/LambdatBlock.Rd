\name{LambdatBlock}
\alias{LambdatBlock}
\title{Create digonal block on transposed relative covariance factor}
\usage{
  LambdatBlock(nl, nc)
}
\arguments{
  \item{nl}{Number of levels in a grouping factor for a
  particular random effects term (the number of levels in
  the \code{grp} argument in \code{\link{Zsection}}).}

  \item{nc}{Number of columns in a dense model matrix for a
  particular random effects term (the number of columns in
  the \code{mm} argument in \code{\link{Zsection}}).}
}
\value{
  A \code{list} with: \itemize{ \item the block \item
  ititial values of theta for the block \item lower bounds
  on these initial theta values \item a function that
  updates the block given the section of theta for this
  block }
}
\description{
  Each random-effects term is represented by diagonal block
  on the transposed relative covariance factor.
  \code{LambdatBlock} creates such a block, and returns
  related information along with it.
}
\details{
  FIXME: change the name of this function to have proper
  camelCase
}
\examples{
(l <- LambdatBlock(2, 3))
within(l, slot(Lambdat, 'x') <- updateLambdatx(as.numeric(10:12)))
}

