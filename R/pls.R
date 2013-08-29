##' @importMethodsFrom Matrix t %*% crossprod diag tcrossprod solve determinant update
##' @importFrom Matrix bdiag rBind Diagonal Cholesky sparse.model.matrix
##' @importFrom lme4 findbars nobars subbars
NULL

##' Create linear mixed model deviance function
##'
##' A pure \code{R} implementation of the
##' penalized least squares (PLS) approach for computing
##' linear mixed model deviances. The purpose
##' is to clarify how PLS works without having
##' to read through C++ code, and as a sandbox for
##' trying out modifications to PLS.
##'
##' @param X output of \code{lFormula} or a model matrix
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
                offset = numeric(n),REML = TRUE,...) {
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

##' Create linear mixed model deviance function from formula/data specification
##'
##' A pure \code{R} implementation of the
##' penalized least squares (PLS) approach to evaluation of the deviance or the
##' REML criterion for linear mixed-effects models.
##'
##' @param formula a two-sided model formula with random-effects terms
##'   and, optionally, fixed-effects terms.
##' @param data a data frame in which to evaluate the variables from \code{form}
##' @param weights prior weights
##' @param offset offset
##' @param REML calculate REML deviance?
##' @param sparseX should X, the model matrix for the fixed-effects coefficients be sparse?
##' @param ... additional arguments
##' @keywords models
##'
##' @return a function that evaluates the deviance or REML criterion
##' @export
plsform <- function(formula, data, REML=TRUE, weights, offset, sparseX = FALSE, ...)  {
    stopifnot(inherits(formula, "formula"), length(formula) == 3L,
              length(rr <- findbars(formula[[3]])) > 0L)
    mc <- match.call()
    fr <- eval(mkMFCall(mc, formula), parent.frame())      # evaluate the model frame
    fr1 <- eval(mkMFCall(mc, formula, TRUE), parent.frame())
    trms <- attr(fr, "terms") <- attr(fr1, "terms")
    c(list(X = if (sparseX) sparse.model.matrix(trms,fr) else model.matrix(trms,fr),
           y = model.response(fr, type="numeric"),
           fr = fr, call = mc,
           REML = as.logical(REML)[1]),
      if (is.null(wts <- model.weights(fr))) wts else list(weights=wts),
      if (is.null(off <- model.offset(fr))) off else list(offset=off),
      mkLambdat(lapply(rr, function(t) as.factor(eval(t[[3]], fr))),
                lapply(rr, function(t)
                       model.matrix(eval(substitute( ~ foo, list(foo = t[[2]]))), fr))))
}

## Create a section of Zt from a sparse matrix of indicators, zt, and
## a dense model matrix, mm. When mm has a single column, just
## multiply zt by a diagonal matrix. For multiple columns, rBind the
## products of zt and a diagonal matrix for each column, then
## rearrange the order of the rows.
Zsection <- function(zt,mm) {
    if ((m <- ncol(mm)) == 1L) return(zt %*% Diagonal(x=mm))
    ## Row indices.  If m = 2, nrow(zt) = 10 we want the order 1,11,2,12,3,13,...,20
    rinds <- as.vector(matrix(seq_len(m*nrow(zt)), nrow=m, byrow=TRUE))
    do.call(rBind,lapply(seq_len(m), function(j) zt %*% Diagonal(x=mm[,j])))[rinds,]
}

## Create the diagonal block on Lambdat for a random-effects term with
## nc columns and nl levels.  The value is a list with the starting
## value of theta for the block, the lower bounds, the block of
## Lambdat and the function that updates the block given the section
## of theta for this block.
Lambdatblock <- function(nc, nl) {
    if (nc == 1L)
        return(list(theta = 1,
                    lower = 0,
                    Lambdat = Diagonal(x = rep(1, nl)),
                    updateLambdatx=local({nl <- nl;function(theta) rep.int(theta[1],nl)})))
    m <- diag(nrow=nc, ncol=nc); ut <- upper.tri(m, diag=TRUE)
    theta <- m[ut]; m[ut] <- seq_along(theta)
    Lambdat <- do.call(bdiag, lapply(seq_len(nl), function(i) m))
    list(theta=theta,
         lower=ifelse(theta,0,-Inf),
         Lambdat=Lambdat,
         updateLambdatx = local({
             Lind <- Lambdat@x
             function(theta) theta[Lind]
         })
         )
}

## Create the call to model.frame from the matched call with the
## appropriate substitutions and eliminations
mkMFCall <- function(mc, form, nobars=FALSE) {
    m <- match(c("data", "subset", "weights", "na.action", "offset"), names(mc), 0)
    mc <- mc[c(1, m)]                   # retain only a subset of the arguments
    mc$drop.unused.levels <- TRUE       # ensure unused levels of factors are dropped
    mc[[1]] <- quote(model.frame)       # change the call to model.frame
    form[[3]] <- if (nobars) {
        if(is.null(nb <- nobars(form[[3]]))) 1 else nb
    } else subbars(form[[3]])
    mc$formula <- form
    mc
}

##' Make transposed random effects design matrix and relative covariance factor
##'
##' Create Lambdat as a ddiMatrix or a dgCMatrix and Zt as a dgCMatrix.
##' 
##' @param grps List of factor vectors of length n indicating groups.
##' @param mms List of model matrices.
##' @details Each list element in \code{grps} and \code{mms} corresponds to
##'   a random effects term.
##' @return A \code{list} with:
##' \itemize{
##' \item \code{Lambdat} Transformed relative covariance factor
##' \item \code{Zt} Transformed random effects design matrix
##' \item \code{theta} Vector of covariance parameters
##' \item \code{lower} Vector of lower bounds on the covariance parameters
##' \item \code{upper} Vector of upper bounds on the covariance parameters
##' \item \code{thfun} A function that maps \code{theta} onto the non-zero
##'   elements of \code{Lambdat}
##' }
##' @export
mkLambdat <- function(grps, mms) {
    ll <- list(Zt = do.call(rBind, mapply(Zsection, lapply(grps, as, Class="sparseMatrix"), mms)))
    nl <- sapply(grps, function(g) length(levels(g)))
    nc <- sapply(mms, ncol)
    if (all(nc == 1L)) {  # for scalar models use a diagonal Lambdat and a simpler thfun
        nth <- length(nc)
        ll$lower <- numeric(nth)
        ll$theta <- rep.int(1, nth)
        ll$upper <- rep.int(Inf, nth)
        ll$Lambdat <- Diagonal(x=rep(1,sum(nl)))
        ll$thfun <- local({
            nlvs <- nl
            function(theta) rep.int(theta,nlvs)})
    } else {
        zz <- mapply(Lambdatblock, nc, nl, SIMPLIFY=FALSE)
        ll$Lambdat <- do.call(bdiag, lapply(zz, "[[", "Lambdat"))
        th <- lapply(zz, "[[", "theta")
        ll$theta <- unlist(th)
        ll$lower <- sapply(zz, "[[", "lower")
        ll$upper <- rep.int(Inf, length(ll$theta))
        ll$thfun <- local({
            splits <- rep.int(seq_along(th), sapply(th, length))
            thfunlist <- lapply(zz,"[[","updateLambdatx")
            function (theta)
                unlist(mapply(do.call,thfunlist,lapply(split(theta,splits),list)))
        })
    }
    ll
}
