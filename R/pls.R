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
      mkRanefRepresentation(lapply(rr, function(t) as.factor(eval(t[[3]], fr))),
                lapply(rr, function(t)
                       model.matrix(eval(substitute( ~ foo, list(foo = t[[2]]))), fr))))
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

##' Create a section of a transposed random effects design matrix
##'
##' @param grp Grouping factor for a particular random effects term.
##' @param mm Dense model matrix for a particular random effects term.
##' @return Dense section of a random effects design matrix
##' @export
##' @examples
##' ## number of observations, n = 6
##' ## number of levels, nl = 3
##' ## number of columns ('predictors'), nc = 2
##' (X <- cbind("(Intercept)"=1,x=1:6)) # an intercept in the first column and
##' 1:6 predictor in the other
##' (g <- as.factor(letters[rep(1:3,2)])) # grouping factor
##' nrow(X) # n = 6
##' nrow(X) == length(g) # and consistent n between X and g
##' ncol(X) # nc = 2
##' nlevels(g) # nl = 3
##' Zsection(g, X)
Zsection <- function(grp,mm) {
                                        # zt is a sparse matrix of indicators to groups.
                                        # this is an interesting feature of coercing a
                                        # factor to a sparseMatrix.
    zt <- as(as.factor(grp), Class="sparseMatrix")
                                        # if mm has one column, multiply zt by a diagonal
                                        # matrix.
    if ((m <- ncol(mm)) == 1L) return(zt %*% Diagonal(x=mm))
                                        # if mm has more than one column, carry on.
                                        # figure out how to rearrange the order of the
                                        # rows by calculating row indices (rinds)
                                        # eg: if m = 2, nrow(zt) = 10, we want the order:
                                        #     1,11,2,12,3,13,...,20
    rinds <- as.vector(matrix(seq_len(m*nrow(zt)), nrow=m, byrow=TRUE))
                                        # rBind products of zt and a diagonal matrix
                                        # for each column, then rearrange rows.
    do.call(rBind,lapply(seq_len(m), function(j) zt %*% Diagonal(x=mm[,j])))[rinds,]
}

## Create the diagonal block on Lambdat for a random-effects term with
## nc columns and nl levels.  The value is a list with the starting
## value of theta for the block, the lower bounds, the block of
## Lambdat and the function that updates the block given the section
## of theta for this block.

##' Create digonal block on transposed relative covariance factor
##'
##' Each random-effects term is represented by diagonal block on
##' the transposed relative covariance factor. \code{Lambdatblock}
##' creates such a block, and returns related information along
##' with it.
##'
##' @param nc Number of columns in a dense model matrix for a particular
##' random effects term.
##' @param nl Number of levels in a grouping factor for a particular
##' random effects term.
##' @return A \code{list} with:
##' \itemize{
##' \item the block
##' \item ititial values of theta for the block
##' \item lower bounds on these initial theta values
##' \item a function that updates the block given the section
##'   of theta for this block
##' }
##' @details
##' FIXME:  change the name of this function to have proper camelCase
##' FIXME:  change the order of the arguments to better match \code{\link{Zsection}}
##' @export
##' @examples
##' (l <- Lambdatblock(2, 3))
##' within(l, slot(Lambdat, 'x') <- updateLambdatx(as.numeric(10:12)))
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

##' Make random effects representation
##'
##' Create all of the elements required to specify the random-effects
##' structure of a mixed effects model.
##' 
##' @param grps List of factor vectors of length n indicating groups.  Each
##' element corresponds to a random effects term.
##' @param mms List of model matrices.  Each
##' element corresponds to a random effects term.
##' @details
##' The basic idea of this function is to call \code{\link{Zsection}} and
##' \code{\link{Lambdatblock}} once for each random effects term (ie.
##' each list element in \code{grps} and \code{mms}). The results of
##' \code{\link{Zsection}} for each term are \code{rBind}ed together.
##' The results of \code{\link{Lambdatblock}} are \code{bdiag}ed
##' together, unless all terms have only a single column ('predictor')
##' in which case a diagonal matrix is created directly.
##' 
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
mkRanefRepresentation <- function(grps, mms) {
                                        # compute transposed random effects design
                                        # matrix, Zt (Class="dgCMatrix"), by
                                        # rBinding the sections for each term.
    ll <- list(Zt = do.call(rBind, mapply(Zsection, grps, mms)))
                                        # number of levels in each grouping factor
    nl <- sapply(grps, function(g) length(levels(g)))
                                        # number of columns in each model matrix
    nc <- sapply(mms, ncol)
                                        # for scalar models use a diagonal Lambdat
                                        # (Class="ddiMatrix") and a simpler thfun. 
    if (all(nc == 1L)) {  
        nth <- length(nc)
        ll$lower <- numeric(nth)
        ll$theta <- rep.int(1, nth)
        ll$upper <- rep.int(Inf, nth)
        ll$Lambdat <- Diagonal(x=rep(1,sum(nl)))
        ll$thfun <- local({
            nlvs <- nl
            function(theta) rep.int(theta,nlvs)})
                                        # for vector models bdiag the lambdat
                                        # blocks together (Class="dgCMatrix")
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
