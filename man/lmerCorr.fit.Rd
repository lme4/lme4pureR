\name{lmerCorr.fit}
\alias{lmerCorr.fit}
\title{lmer.fit function for the single correlation template model}
\usage{
lmerCorr.fit(y, mmFE, corr, grp, weights, offset = numeric(n), REML = TRUE)
}
\arguments{
  \item{y}{response vector}

  \item{mmFE}{model matrix for the fixed effects}

  \item{corr}{template correlation matrix for the single
  scalar random effect}

  \item{grp}{grouping factor for the random effect (levels
  correspond to \code{dimnames} of \code{corr}}

  \item{weights}{weights}

  \item{offset}{offset}

  \item{REML}{should restricted maximum likelihood be
  used?}
}
\description{
lmer.fit function for the single correlation template model
}
\examples{
library(lme4pureR)
library(lme4)
library(minqa)
library(subscript)
library(rmv)
set.seed(1)
n <- 100
x <- rnorm(n)
X <- cbind(1, x)
grp <- gl(n/5,5)
ugrps <- unique(grp)
q <- length(ugrps)
corr <- rcov(q+1,q)
dimnames(corr) <- rep(list(as.character(ugrps)), 2)
b <- as.numeric(rmv(1,corr) \%*\% as(grp, "sparseMatrix"))
y <- as.numeric(X\%*\%rnorm(ncol(X)) + b + rnorm(n))

m <- lmerCorr.fit(y, X, corr, grp)
m$par
}

