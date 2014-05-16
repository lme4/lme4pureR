##' Linear mixed model deviance function as it
##' appears in the pseudocode of the JSS article
##' 
##' A pure \code{R} implementation of the
##' penalized least squares (PLS) approach for computing
##' linear mixed model deviances. The purpose
##' is to clarify how PLS works without having
##' to read through C++ code, and as a sandbox for
##' trying out modifications to PLS.
##'
##' @param X fixed effects model matrix
##' @param y response
##' @param Zt transpose of the sparse model matrix for the random effects
##' @param Lambdat upper triangular sparse Cholesky factor of the
##'    relative covariance matrix of the random effects
##' @param mapping a function that takes a value of \code{theta} and produces
##'    the non-zero elements of \code{Lambdat}.  The structure of \code{Lambdat}
##'    cannot change, only the numerical values
##' @param weights prior weights
##' @param offset offset
##' @param REML calculate REML deviance?
##' @param ... additional arguments
##' @keywords models
##'
##' @return a function that evaluates the deviance or REML criterion
##' @export
plsJSS <- function(X, y, Zt, Lambdat, mapping, weights,
                   offset = numeric(n), REML = TRUE, ...)
{
    # SW: how to test for sparse matrices, without specifying the specific class?
    stopifnot(is.matrix(X)) #  is.matrix(Zt), is.matrix(Lambdat))
    n <- length(y); p <- ncol(X); q <- nrow(Zt)
    stopifnot(nrow(X) == n, ncol(Zt) == n,
              nrow(Lambdat) == q, ncol(Lambdat) == q)
                                        # calculate weighted products
    sqrtW <- if (missing(weights)) Diagonal(n=n) else Diagonal(x=sqrt(as.numeric(weights)))
    WX <- sqrtW %*% X
    Wy <- sqrtW %*% y
    ZtW <- Zt %*% sqrtW
    XtWX <- crossprod(WX)
    XtWy <- crossprod(WX, Wy)
    ZtWX <- ZtW %*% WX
    ZtWy <- ZtW %*% Wy
    rm(WX,Wy)
    local({                             # mutable values stored in local environment
        b <- numeric(q)                 # conditional mode of random effects
        beta <- numeric(p)              # conditional estimate of fixed-effects
        cu <- numeric(q)                # intermediate solution
        RXtRX <- XtWX                   # down-dated XtWX
        L <- Cholesky(tcrossprod(Lambdat %*% ZtW), LDL = FALSE, Imult=1)
        Lambdat <- Lambdat              # stored here b/c x slot will be updated
        mu <- numeric(n)                # conditional mean of response
        RZX <- matrix(0,nrow=q,ncol=p)  # intermediate matrix in solution
        u <- numeric(q)                 # conditional mode of spherical random effects
        degFree <- as.numeric(n)        # degrees of freedom (depends on REML)
        if(REML) degFree <- degFree - as.numeric(p)
        function(theta) {

            ##################################################
            # Step I: update covariance parameters
            ##################################################
                                        # update relative covariance factor
                                        # by placing the new values of theta
                                        # in the appropriate positions
            Lambdat@x[] <<- mapping(theta)
                                        # update random-effects Cholesky factor
            L <<- update(L, Lambdat %*% ZtW, mult = 1)
            
            ##################################################
            # Step II: solve normal equations
            ##################################################
                                        # solve eqn. ??
            cu[] <<- as.vector(solve(L, solve(L, Lambdat %*% ZtWy,
                                              system="P"), system="L"))
                                        # solve eqn. ??
            RZX[] <<- as.vector(solve(L, solve(L, Lambdat %*% ZtWX,
                                               system="P"), system="L"))
                                        # downdate XtWX and form Cholesky
                                        # factor (eqn. ??)
            RXtRX <<- as(XtWX - crossprod(RZX), "dpoMatrix")
                                        # conditional estimate of fixed-effects
                                        # coefficients (solve eqn. ??)
            beta[] <<- as.vector(solve(RXtRX, XtWy - crossprod(RZX, cu)))
                                        # conditional mode of the spherical
                                        # random-effects coefficients (eqn. ??)
            u[] <<- as.vector(solve(L, solve(L, cu - RZX %*% beta,
                                             system = "Lt"), system="Pt"))
                                        # update conditional model of the
                                        # non-spherical random-effects
                                        # coefficients
            b[] <<- as.vector(crossprod(Lambdat,u))

            
            ##################################################
            # Step III: update linear predictor and residuals
            ##################################################
                                        # update linear predictor
            mu[] <<- as.vector(crossprod(Zt,b) + X %*% beta + offset)
                                        # weighted residuals
            wtres <- sqrtW*(y-mu)
                                        # penalized, weighted residual
                                        # sum-of-squares


            ##################################################
            # Step IV: compute profiled deviance
            ##################################################
            pwrss <- sum(wtres^2) + sum(u^2)
                                        # log determinant (depends on
                                        # whether REML or ML is used)
            logDet <- 2*determinant(L, logarithm = TRUE)$modulus 
            if (REML) logDet <- logDet + determinant(RXtRX,
                                                     logarithm = TRUE)$modulus
            attributes(logDet) <- NULL
                                        # profiled deviance or REML criterion
            profDev <- logDet + degFree*(1 + log(2*pi*pwrss) - log(degFree))
            return(profDev)
        }
    })
}
