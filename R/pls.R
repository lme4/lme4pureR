##' @importMethodsFrom Matrix t %*% crossprod diag tcrossprod solve determinant update
##' @importFrom Matrix bdiag rBind Diagonal Cholesky sparse.model.matrix
##' @importFrom lme4 findbars nobars subbars
NULL

##' Create linear mixed model deviance function
##'
##' A pure \code{R} implementation of the penalized least squares
##' (PLS) approach for computing linear mixed model deviances. The
##' purpose is to clarify how PLS works without having to read through
##' C++ code, and as a sandbox for trying out modifications to PLS.
##'
##' @param X fixed effects design (model) matrix
##' @param y response
##' @param Zt transpose of the sparse model matrix for the random effects
##' @param Lambdat upper triangular sparse Cholesky factor of the
##'    relative covariance matrix of the random effects
##' @param thfun a function that takes a value of \code{theta} and produces
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
pls <- function(X,y,Zt,Lambdat,thfun,weights,
                offset = numeric(n),REML = TRUE,...)
{
    stopifnot(is.matrix(X), is.matrix(Zt), is.matrix(Lambdat))
    n <- length(y); p <- ncol(X); q <- nrow(Zt)
    stopifnot(nrow(X) == n, ncol(Zt) == n,
              nrow(Lambdat) == q, ncol(Lambdat) == q, is.function(thfun))
                                        # calculate weighted products
    Whalf <- if (missing(weights)) Diagonal(n=n) else Diagonal(x=sqrt(as.numeric(weights)))
    WX <- Whalf %*% X
    Wy <- Whalf %*% y
    ZtW <- Zt %*% Whalf
    XtWX <- crossprod(WX)
    XtWy <- crossprod(WX, Wy)
    ZtWX <- ZtW %*% WX
    ZtWy <- ZtW %*% Wy
    rm(WX,Wy,ZtW)
    local({                             # mutable values stored in local environment
        b <- numeric(q)                 # conditional mode of random effects
        beta <- numeric(p)              # conditional estimate of fixed-effects
        cu <- numeric(q)                # intermediate solution
        DD <- XtWX                      # down-dated XtWX
        L <- Cholesky(tcrossprod(Lambdat %*% Zt), LDL = FALSE, Imult=1)
        Lambdat <- Lambdat              # stored here b/c x slot will be updated
        mu <- numeric(n)                # conditional mean of response
        RZX <- matrix(0,nrow=q,ncol=p)  # intermediate matrix in solution
        u <- numeric(q)                 # conditional mode of spherical random effects
        function(theta) {
            Lambdat@x[] <<- thfun(theta)
            L <<- update(L, Lambdat %*% Zt, mult = 1)
                                        # solve eqn. 30
            cu[] <<- as.vector(solve(L, solve(L, Lambdat %*% ZtWy, system="P"),
                                     system="L"))
                                        # solve eqn. 31
            RZX[] <<- as.vector(solve(L, solve(L, Lambdat %*% ZtWX, system="P"),
                                      system="L"))
            ## downdate XtWX and form Cholesky factor (eqn. 32)
            DD <<- as(XtWX - crossprod(RZX), "dpoMatrix")
            ## conditional estimate of fixed-effects coefficients (solve eqn. 33)
            beta[] <<- as.vector(solve(DD, XtWy - crossprod(RZX, cu)))
            ## conditional mode of the spherical random-effects coefficients (eqn. 34)
            u[] <<- as.vector(solve(L, solve(L, cu - RZX %*% beta, system = "Lt"),
                                    system="Pt"))
            b[] <<- as.vector(crossprod(Lambdat,u))
                                        # conditional mean of the response
            mu[] <<- as.vector(crossprod(Zt,b) + X %*% beta + offset)
            wtres <- Whalf*(y-mu)       # weighted residuals
            pwrss <- sum(wtres^2) + sum(u^2) # penalized, weighted residual sum-of-squares
            fn <- as.numeric(length(mu))
            ld <- 2*determinant(L,logarithm=TRUE)$modulus # log determinant
            if (REML) {
                ld <- ld + determinant(DD,logarithm=TRUE)$modulus
                fn <- fn - length(beta)
            }
            attributes(ld) <- NULL
                                        # profiled deviance or REML criterion
            ld + fn*(1 + log(2*pi*pwrss) - log(fn))
        }
    })
}

##' Create the structure of a linear mixed model from formula/data specification
##'
##' A pure \code{R} implementation of the penalized least squares
##' (PLS) approach to evaluation of the deviance or the REML criterion
##' for linear mixed-effects models.
##'
##' @param formula a two-sided model formula with random-effects terms
##'   and, optionally, fixed-effects terms.
##' @param data a data frame in which to evaluate the variables from \code{form}
##' @param REML calculate REML deviance?
##' @param weights prior weights
##' @param offset offset
##' @param sparseX should X, the model matrix for the fixed-effects coefficients be sparse?
##' @param ... additional arguments
##' @keywords models
##'
##' @return a \code{list} with:
##' \itemize{
##' \item \code{X} Fixed effects design (model) matrix
##' \item \code{y} Observed response vector
##' \item \code{fr} Model frame
##' \item \code{call} Matched call
##' \item \code{REML} Logical indicating REML or not
##' \item \code{weights} Prior weights or \code{NULL}
##' \item \code{offset} Prior offset term or \code{NULL}
##' \item \code{Zt} Transposed random effects design matrix
##' \item \code{Lambdat} Transposed relative covariance factor
##' \item \code{theta} Vector of covariance parameters
##' \item \code{lower} Vector of lower bounds for \code{theta}
##' \item \code{upper} Vector of upper bounds for \code{theta}
##' \item \code{thfun} A function that maps \code{theta} into the structural non-zero
##' elements of \code{Lambdat}, which are stored in \code{slot(Lambdat, 'x')}
##' }
##' @export
##' @examples
##' form <- Reaction ~ Days + (Days|Subject)
##' data(sleepstudy, package="lme4")
##' ll <- plsform(form, sleepstudy, REML=FALSE)
##' names(ll)
plsform <- function(formula, data, REML=TRUE, weights, offset, sparseX = FALSE, ...)  {
    stopifnot(inherits(formula, "formula"), length(formula) == 3L,
              length(rr <- findbars(formula[[3]])) > 0L)
    mc <- match.call()
    fr <- eval(mkMFCall(mc, formula), parent.frame())         # evaluate the model frame
    fr1 <- eval(mkMFCall(mc, formula, TRUE), parent.frame())  # fixed-effects model.frame
    trms <- attr(fr, "terms") <- attr(fr1, "terms")
    c(list(X = if (sparseX) sparse.model.matrix(trms,fr) else model.matrix(trms,fr),
           y = model.response(fr,type="numeric"),
           fr = fr, call = mc,
           REML = as.logical(REML)[1]),
      if (is.null(wts <- model.weights(fr))) wts else list(weights=wts),
      if (is.null(off <- model.offset(fr))) off else list(offset=off),
      mkRanefRepresentation(lapply(rr, function(t) as.factor(eval(t[[3]], fr))),
                            lapply(rr, function(t)
                                        model.matrix(eval(substitute( ~ foo,
                                                                     list(foo = t[[2]]))), fr))))
}
 
