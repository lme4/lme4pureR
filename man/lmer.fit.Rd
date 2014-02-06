\name{lmer.fit}
\alias{lmer.fit}
\title{Fit a mixed effects from raw matrices, vectors and grouping factors}
\usage{
lmer.fit(y, mmFE, mmRE, grp, weights, offset = numeric(n), REML = TRUE)
}
\arguments{
  \item{y}{response vector}

  \item{mmFE}{model matrix for the fixed effects}

  \item{mmRE}{template model matrix for the random effects
  (or optionally a list of such matrices)}

  \item{grp}{grouping factor for the random effects (or
  optionally a list of such factors)}

  \item{weights}{weights}

  \item{offset}{offset}

  \item{REML}{should restricted maximum likelihood be
  used?}
}
\description{
The only output from this function is the result of an
optimization over the covariance parameters.
}
\examples{
library(lme4pureR)
library(lme4)
library(minqa)
set.seed(1)
n <- 1000
x <- rnorm(n)
z <- rnorm(n)
X <- cbind(1, x)
ZZ <- cbind(1, z)
grp <- gl(n/5,5)
RE <- mkRanefStructures(list(grp), list(ZZ))
Z <- t(RE$Zt)
y <- as.numeric(X\%*\%rnorm(ncol(X)) + Z\%*\%rnorm(ncol(Z)) + rnorm(n))
m <- lmer.fit(y,X,ZZ,grp)
m$par
Lambdat <- RE$Lambdat
Lambdat
}

