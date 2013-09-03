##' Create an approximate deviance evaluation function for GLMMs using Laplace
##'
##' A pure \code{R} implementation of the penalized iteratively
##' reweighted least squares (PIRLS) algorithm for computing
##' generalized linear mixed model deviances. The purpose is to
##' clarify how PIRLS works without having to read through C++ code,
##' and as a sandbox for trying out modified versions of PIRLS.
##'
##' @param glmod output of \code{glFormula}
##' @param y response
##' @param eta linear predictor
##' @param family a \code{glm} family object
##' @param weights prior weights
##' @param offset offset
##' @param tol convergence tolerance
##' @param npirls maximum number of iterations
##' @param nAGQ either 0 (PIRLS for \code{u} and \code{beta}) or 1 (\code{u} only).
##'     currently no quadature is available
##' @param verbose verbose
##'
##' @details \code{pirls1} is a convenience function for optimizing
##' \code{pirls} under \code{nAGQ = 1}. In particular, it wraps
##' \code{theta} and \code{beta} into a single argument
##' \code{thetabeta}.
##'
##' @return A function for evaluating the GLMM Laplace approximated deviance
##' @export
pirls <- function(X, y, Zt, Lambdat, thfun, theta, nobs, n, fr, family,
                  call, mustart, weights, etastart, offset=numeric(nobs),
                  tol = 10^-6, npirls = 30,nstephalf = 10,nAGQ = 1,verbose=0L,
                  ...){
    n <- length(y); p <- ncol(X); q <- nrow(Zt)
    stopifnot(nrow(X) == n, ncol(Zt) == n, inherits(family, "family"),
              nrow(Lambdat) == q, ncol(Lambdat) == q, is.function(thfun))
    local({    
        nth <- max(thfun(theta))
        betaind <- -seq_len(nth) # indices to drop 1:nth
        linkinv <- family$linkinv
        variance <- family$variance
        muEta <- family$mu.eta
        aic <- family$aic
        sqDevResid <- family$dev.resid
        eta <- etastart
        mu <- linkinv(eta)
        beta <- numeric(p)
        u <- numeric(q)
        L <- Cholesky(tcrossprod(Lambdat %*% Zt), perm=FALSE, LDL=FALSE, Imult=1)
        if (nAGQ > 0L) {
                                        # create function for conducting PIRLS
            function(thetabeta) {
                                        # initialize
                Lambdat@x[] <<- thfun(thetabeta[-betaind])
                LtZt <- Lambdat %*% Zt
                beta[] <<- thetabeta[betaind]
                offb <- offset + X %*% beta
                updatemu <- function(uu) {
                    eta[] <<- offb + as.vector(crossprod(LtZt, uu))
                    mu[] <<- linkinv(eta)
                    sum(sqDevResid(y, mu, weights)) + sum(uu^2)
                }
                u[] <<- numeric(q)
                olducden <- updatemu(u)
                cvgd <- FALSE
                for(i in 1:npirls){
                                        # update w and muEta
                    Whalf <- Diagonal(x=sqrt(weights/variance(mu)))
                                        # update weighted model matrix
                    LtZtMWhalf <- LtZt %*% (Diagonal(x=muEta(eta)) %*% Whalf)
                                        # update Cholesky decomposition
                    L <<- update(L, LtZtMWhalf, 1)
                                        # alternative (more explicit but slower)
                                        # Cholesky update
                    # L <- Cholesky(tcrossprod(LtZtMWhalf), perm=FALSE, LDL=FALSE, Imult=1)
                                        # update weighted residuals
                    wtres <- Whalf %*% (y - mu)
                                        # solve for the increment
                    delu <- as.vector(solve(L, LtZtMWhalf %*% wtres - u))
                    if (verbose > 0L) {
                        cat(sprintf("inc: %12.4g", delu[1]))
                        nprint <- min(5,length(delu))
                        for (j in 2:nprint) cat(sprintf(" %12.4g", delu[j]))
                        cat("\n")
                    }
                                        # update mu and eta and calculate
                                        # new unscaled conditional log density
                    ucden <- updatemu(u + delu)
                    if (verbose > 1L) {
                        cat(sprintf("%6.4f: %10.3f\n", 1, ucden))
                    }

                    if(abs((olducden - ucden) / ucden) < tol){
                        cvgd <- TRUE
                        break
                    }
                                        # step-halving
                    if(ucden > olducden){
                        for(j in 1:nstephalf){
                            ucden <- updatemu(u + (delu <- delu/2))
                            if (verbose > 1L) {
                                cat(sprintf("%6.4f: %10.3f\n", 1/2^j, ucden))
                            }
                            if(ucden < olducden) break
                        }
                        if(ucden > olducden) stop("Step-halving failed")
                    }
                                        # set old unscaled conditional log density
                                        # to the new value
                    olducden <- ucden
                                        # update the conditional modes (take a step)
                    u[] <<- u + delu
                }
                if(!cvgd) stop("PIRLS failed to converge")
                                        # create Laplace approx to -2log(L)
                ldL2 <- 2*determinant(L, logarithm = TRUE)$modulus
                attributes(ldL2) <- NULL
                # FIXME: allow for quadrature approximations too
                Lm2ll <- aic(y,n,mu,weights,NULL) + sum(u^2) + ldL2

                if (verbose > 0L) {
                    cat(sprintf("%10.3f: %12.4g", Lm2ll, thetabeta[1]))
                    for (j in 2:length(thetabeta)) cat(sprintf(" %12.4g", thetabeta[j]))
                    cat("\n")
                }

                Lm2ll
            }
        } else stop("code for nAGQ == 0 needs to be added")
    })
}

## All the silly dancing around that must be done because of the brain-dead
## design of R's family class, especially the initialize component.
initializeResp <- function(fr, family){
    y <- model.response(fr)
    if(length(dim(y)) == 1) {           # handle this obscure case
        nm <- rownames(y)
        dim(y) <- NULL
        if(!is.null(nm)) names(y) <- nm
    }
    nobs <- nrow(fr)
    etastart <- model.extract(fr, "etastart")
    mustart <- model.extract(fr, "mustart")
    weights <- if (is.null(wts <- model.weights(fr))) rep.int(1,nobs) else wts
    eval(family$initialize)
    rm(wts); family$initialize <- NULL # remove clutter from str output
    if (is.null(etastart)) etastart <- family$linkfun(mustart)
    as.list(environment())
}

pirlsform <- function(formula, data, family=binomial, weights, offset, sparseX = FALSE, ...)  {
    mc <- match.call()
    if (!inherits(family, "family")) family <- as.function(family)()
    stopifnot(inherits(family, "family"),
              inherits(formula, "formula"),
              length(formula) == 3L,
              length(rr <- findbars(formula[[3]])) > 0L)
    fr <- eval(mkMFCall(mc, formula), parent.frame())         # evaluate the model frame
    fr1 <- eval(mkMFCall(mc, formula, TRUE), parent.frame())  # fixed-effects model frame
    trms <- attr(fr, "terms") <- attr(fr1, "terms")
    c(initializeResp(fr, family),
      list(X = if (sparseX) sparse.model.matrix(trms,fr) else model.matrix(trms,fr),
           call = mc),
      if (is.null(off <- model.offset(fr))) off else list(offset=off),
      mkRanefRepresentation(lapply(rr, function(t) as.factor(eval(t[[3]], fr))),
                            lapply(rr, function(t)
                                        model.matrix(eval(substitute( ~ foo,
                                                                     list(foo = t[[2]]))), fr))))
}
