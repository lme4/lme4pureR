##' @importMethodsFrom Matrix t %*% crossprod diag tcrossprod solve determinant update
##' @importFrom Matrix bdiag rBind Diagonal Cholesky
##' @importFrom minqa bobyqa
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
pls <- function(X,y,Zt,Lambdat,thfun,
                weights,
                offset = numeric(n),
                REML = TRUE,...){
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
                                        # other values in the environment
    L <- Cholesky(tcrossprod(Lambdat %*% Zt), LDL = FALSE, Imult=1)
    beta <- numeric(p)
    u <- numeric(q)
    mu <- numeric(n)
    DD <- NULL
    RZX <- matrix(0,nrow=q,ncol=p)
    cu <- numeric(q)
    rm(WX, Wy, ZtW)
    function(theta) {
        Lambdat@x[] <<- thfun(theta)
        L <<- update(L, Lambdat %*% Zt, mult = 1)
        ## solve system from equation 30
        cu[] <<- as.vector(solve(L, solve(L, Lambdat %*% ZtWy, system="P"),
                                system="L"))
        ## solve system from eqn. 31
        RZX[] <<- as.vector(solve(L, solve(L, Lambdat %*% ZtWX, system="P"),
                                  system="L"))
        ## downdate XtWX and form Cholesky factor (eqn. 32)
        DD <<- as(XtWX - crossprod(RZX), "dpoMatrix")
        ## conditional estimate of fixed-effects coefficients (solve eqn. 33)
        beta[] <<- as.vector(solve(DD, XtWy - crossprod(RZX, cu)))
        ## conditional mode of the spherical random-effects coefficients (eqn. 34)
        u[] <<- as.vector(solve(L, solve(L, cu - RZX %*% beta, system = "Lt"), system="Pt"))
                                        # conditional mean of the response
        ## crossprod(Zt,crossprod(Lambdat,u))  == Z Lambda u == Z b
        mu[] <<- as.vector(crossprod(Zt,crossprod(Lambdat,u)) + X %*% beta + offset)

        wtres <- Whalf*(y-mu)      # weighted residuals
        wrss <- sum(wtres^2)       # weighted residual sums of squares
        pwrss <- wrss + sum(u^2)   # plus penalty
        ldL2 <- 2*determinant(L,logarithm=TRUE)$modulus # log determinants
        ldRX2 <- determinant(DD,logarithm=TRUE)$modulus
        attributes(ldL2) <- attributes(ldRX2) <- NULL
                                        # profiled deviance or REML criterion
        if (REML)
            ldL2 + ldRX2 + (n-p)*(1 + log(2*pi*pwrss) - log(n-p))
        else
            ldL2 + n*(1 + log(2*pi*pwrss) - log(n))
    }
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
##' @param ... additional arguments
##' @keywords models
##'
##' @return a function that evaluates the deviance or REML criterion
##' @export
plsform <- function(formula, data, REML=TRUE, weights, offset, ...)  {
    stopifnot(inherits(formula, "formula"), length(formula) == 3L,
              length(rr <- findbars(formula[[3]])) > 0L)
    mf <- mc <- mcout <- match.call()
    m <- match(c("data", "subset", "weights", "na.action", "offset"), names(mf), 0)
    mf <- mf[c(1, m)]; mf$drop.unused.levels <- TRUE; mf[[1]] <- quote(model.frame)
    fr.form <- subbars(formula)         # substitute "+" for "|"
    environment(fr.form) <- environment(formula)
    mf$formula <- fr.form
    fr <- eval(mf, parent.frame())      # evaluate the model frame
    ff <- formula
    ff[[3]] <- if(is.null(nb <- nobars(ff[[3]]))) 1 else nb
    mf$formula <- ff
    fr1 <- eval(mf,parent.frame())
    trms <- attr(mf, "terms") <- attr(fr1, "terms")
    ll <- list(X = model.matrix(trms,fr1),
               y = model.response(fr, type="numeric"),
               fr = fr, REML = as.logical(REML)[1])
    n <- length(ll$y)
    if (!is.null(wts <- model.weights(fr1))) {
        stopifnot((lwts <- length(wts)) %in% c(1L,n))
        ll$weights <- rep.int(wts, n %/% lwts)
    }
    if (!is.null(off <- model.offset(fr1))) {
        stopifnot((loff <- length(off)) %in% c(1L,n))
        ll$weights <- rep.int(off, n %/% loff)
    }
    grps <- lapply(rr, function(t) as.factor(eval(t[[3]], fr)))
    nl <- sapply(grps, function(g) length(levels(g)))
    zsl <- lapply(grps, as, Class="sparseMatrix")
    mms <- lapply(rr, function(t) model.matrix(eval(substitute( ~ foo, list(foo = t[[2]]))), fr))
    ll$Zt <- do.call(rBind, mapply(Zsection, zsl, mms))
    nc <- sapply(mms, ncol)
    if (all(nc == 1L)) {  # for scalar models use a diagonal Lambdat and a simpler thfun
        nth <- length(zsl)
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

## Create a section of Zt from a sparse matrix of indicators, zt, and
## a dense model matrix, mm When mm has a single column it is trivial
## - multiply zt by a diagonal matrix For multiple columns, rBind the
## products of zt and a diagonal matrix for each column, then
## rearrange the order of the rows.
Zsection <- function(zt,mm) {
    if ((m <- ncol(mm)) == 1L) return(zt %*% Diagonal(x=mm))
    ## Row indices.  If m = 2 and nrow(zt) = 10 we want the order 1,11,2,12,3,13,...,20
    rinds <- as.vector(matrix(seq_len(m*nrow(zt)), nrow=m, byrow=TRUE))
    do.call(rBind,lapply(seq_len(m), function(j) zt %*% Diagonal(x=mm[,j])))[rinds,]
}

## Create the diagonal block on Lambdat for a random-effects term with nc columns and nl
## levels.  The value is a list with the starting value of theta for the block, the lower
## bounds, the block of Lambdat and the function that updates the block given the section
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
