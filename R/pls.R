##' @importMethodsFrom Matrix t %*% crossprod diag tcrossprod solve determinant bdiag
##' @importFrom minqa bobyqa
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
##' @param obj output of \code{lFormula} or a model matrix
##' @param y response
##' @param ... Arguments to pass to other functions
##'
##' @return a function that evaluates the deviance or REML criterion
##' @export
pls <- function(obj,y,...) UseMethod("pls")

##' @rdname pls
##' @param Zt transpose of the sparse model matrix for the random effects
##' @param Lambdat upper triangular sparse Cholesky factor of the
##'    relative covariance matrix of the random effects
##' @param thfun a function that takes a value of \code{theta} and produces
##'    the non-zero elements of \code{Lambdat}.  The structure of \code{Lambdat}
##'    cannot change, only the numerical values
##' @param weights prior weights
##' @param offset offset
##' @param REML calculate REML deviance?
##' @keywords models
##' @examples
##' library(lme4)
##' library(nloptwrap)
##' lmod <- lFormula(Reaction ~ Days + (Days|Subject), sleepstudy)
##' devf <- pls(lmod,sleepstudy$Reaction)
##' devf(c(1,0,1))             # starting value
##' bobyqa(c(1, 0, 1), devf, lower=c(0,-Inf,0))[c("par","value")]
##' @method pls matrix
##' @S3method pls matrix
pls.matrix <- function(obj,y,Zt,Lambdat,thfun,
                       weights = rep(1, n),
                       offset = numeric(n),
                       REML = TRUE,...){
    n <- length(y)
    p <- ncol(obj)
    q <- nrow(Zt)
    stopifnot(nrow(obj) == n,
              ncol(Zt) == n,
              nrow(Lambdat) == q,
              ncol(Lambdat) == q,
              is.function(thfun))
                                        # calculate weighted products
    W <- Diagonal(x = weights)
    L <- Cholesky(tcrossprod(Lambdat %*% Zt), LDL = FALSE, Imult=1)
    XtWX <- crossprod(obj, W %*% obj)
    XtWy <- crossprod(obj, W %*% y)
    ZtWX <- Zt %*% (W %*% obj)
    ZtWy <- Zt %*% (W %*% y)
    beta <- numeric(p)
    u <- numeric(q)
    mu <- numeric(n)
    DD <- NULL
    RZX <- matrix(0,nrow=q,ncol=p)
    cu <- numeric(q)
    function(theta) {
        Lambdat@x[] <<- thfun(theta) 
        L <- update(L, Lambdat %*% Zt, mult = 1)
        ## solve system from equation 30
        cu[] <- as.vector(solve(L, solve(L, Lambdat %*% ZtWy, system="P"),
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
        mu[] <<- as.vector(crossprod(Zt,crossprod(Lambdat,u)) + obj %*% beta + offset)
                                        # weighted residuals
        wtres <- sqrt(weights)*(y-mu)
                                        # weighted residual sums of squares
        wrss <- sum(wtres^2)
        pwrss <- wrss + sum(u^2)        # penalize
                                        # log determinants
        ldL2 <- 2*determinant(L,logarithm=TRUE)$modulus
#        ldRX2 <- determinant(DD,logarithm=TRUE)$modulus # need method in Matrix
        attributes(ldL2) <- NULL #attributes(ldRX2) <- NULL
                                        # profiled deviance or REML criterion
        ## if (REML) 
        ##     ldL2 + ldRX2 + (n-p)*(1 + log(2*pi*pwrss) - log(n-p))
        ## else
        ldL2 + n*(1 + log(2*pi*pwrss) - log(n))
    }
}
##' @rdname pls

##' @method pls list
##' @S3method pls list
pls.list <- function(obj,y,weights=rep(1,length(y)),
                     offset=rep(0,length(y)),REML=TRUE,...) {
    retrm <- obj$reTrms
    thfun <- retrm$thfun
    if (is.null(thfun)) thfun <- local({
        Lind <- retrm$Lind
        function(theta) theta[Lind]
    })
    pls.matrix(obj$X,y,retrm$Zt,retrm$Lambdat,thfun,weights,offset,REML)
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
##' @param ... Arguments to pass to other functions
##' @param weights prior weights
##' @param offset offset
##' @param REML calculate REML deviance?
##' @keywords models
##'
##' @return a function that evaluates the deviance or REML criterion
##' @export
plsform <- function(formula, data, REML=TRUE, weights = rep(1, n), offset = numeric(n), ...)
    UseMethod("plsform")

##' @rdname plsform
##' @method plsform formula
##' @S3method plsform formula
plsform.formula <- function(formula, data, REML=TRUE, weights, offset, ...)  {
    stopifnot(length(formula) == 3L,
              length(rr <- lme4::findbars(formula[[3]])) > 0L)
    mf <- mc <- mcout <- match.call()
    m <- match(c("data", "subset", "weights", "na.action", "offset"),
               names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- quote(model.frame)
    fr.form <- lme4::subbars(formula)   # substitute "+" for "|"
    environment(fr.form) <- environment(formula)
    mf$formula <- fr.form
    fr <- eval(mf, parent.frame())
    ff <- formula
    ff[[3]] <- if(is.null(nb <- lme4::nobars(ff[[3]]))) 1 else nb
    mf$formula <- ff
    fr1 <- eval(mf,parent.frame())
    trms <- attr(mf, "terms") <- attr(fr1, "terms")
    ll <- list(obj = model.matrix(trms,fr1),
               y = model.response(fr, type="numeric"),
               REML=as.logical(REML)[1])
    n <- length(ll$y)
    if (!is.null(wts <- model.weights(fr1))) { # FIXME: create a utility function for this
        stopifnot((lwts <- length(wts)) %in% c(1L,n))
        ll$weights <- rep.int(wts, n %/% lwts)
    }
    if (!is.null(off <- model.offset(fr1))) {
        stopifnot((loff <- length(off)) %in% c(1L,n))
        ll$weights <- rep.int(off, n %/% loff)
    }
    grps <- lapply(rr, function(t) as.factor(eval(t[[3]], fr)))
    nlvs <- sapply(grps, function(g) length(levels(g)))
    zsl <- lapply(grps, as, Class="sparseMatrix")
    mms <- lapply(rr, function(t) model.matrix(eval(substitute( ~ foo, list(foo = t[[2]]))), fr))
    nth <- choose(sapply(mms, ncol) + 1L, 2L) # length of theta for each r.e. term
    scalar <- all(nth == 1L)
    if (scalar) {
        ll$Zt <- do.call(Matrix::rBind, mapply(function(zt,mm) zt %*% Matrix::Diagonal(x=mm[,1]),
                                               zsl, mms))
        nthtot <- sum(nth)
        lower <- numeric(nthtot)
        theta <- rep(1, nthtot)
        upper <- rep(Inf, nthtot)
        ll$Lambdat <- do.call(Matrix::bdiag, lapply(nlvs, function(n) Diagonal(x=rep.int(1,n))))
        thfun <- function(theta) as.numeric(unlist(mapply(do.call,thfunlist,
                                                          lapply(split(theta,splits),
                                                                 function(x)list(theta=x)))))
        rho <- new.env()
        rho$thfunlist <- lapply(nlvs, function(x) function(theta) rep.int(theta[1], x))
        rho$splits <- seq_along(nth)
        environment(thfun) <- rho
        ll$thfun <- thfun
    } else stop("code for vector-valued random effects terms not yet written")
    do.call(pls, ll)
}
