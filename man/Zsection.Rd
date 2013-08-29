\name{Zsection}
\alias{Zsection}
\title{Create a section of a transposed random effects design matrix}
\usage{
  Zsection(grp, mm)
}
\arguments{
  \item{grp}{Grouping factor for a particular random
  effects term.}

  \item{mm}{Dense model matrix for a particular random
  effects term.}
}
\value{
  Dense section of a random effects design matrix
}
\description{
  Create a section of a transposed random effects design
  matrix
}
\examples{
## consider a term (x | g) with:
## number of observations, n = 6
## number of levels, nl = 3
## number of columns ('predictors'), nc = 2
(X <- cbind("(Intercept)"=1,x=1:6)) # an intercept in the first column and 1:6 predictor in the other
(g <- as.factor(letters[rep(1:3,2)])) # grouping factor
nrow(X) # n = 6
nrow(X) == length(g) # and consistent n between X and g
ncol(X) # nc = 2
nlevels(g) # nl = 3
Zsection(g, X)
}

