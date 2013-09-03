## Create the call to model.frame from the matched call with the
## appropriate substitutions and eliminations
mkMFCall <- function(mc, form, nobars=FALSE) {
    m <- match(c("data", "subset", "weights", "na.action", "offset",
                 "etastart", "mustart"), names(mc), 0)
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
##' ## consider a term (x | g) with:
##' ## number of observations, n = 6
##' ## number of levels, nl = 3
##' ## number of columns ('predictors'), nc = 2
##' (X <- cbind("(Intercept)"=1,x=1:6)) # an intercept in the first column and 1:6 predictor in the other
##' (g <- gl(3,2,labels=letters[1:3])) # grouping factor
##' stopifnot(nrow(X) == length(g))
##' ncol(X) # nc = 2
##' nlevels(g) # nl = 3
##' Zsection(g, X)
##' @rdname utilities
Zsection <- function(grp,mm) {
    zt <- as(as.factor(grp), Class="sparseMatrix") # convert factor to sparse indicators
    unname(do.call(rBind, lapply(seq_len(nrow(zt)), function(i) t(mm) %*% Diagonal(x=zt[i,]))))
}

## Create the diagonal block on Lambdat for a random-effects term with
## nc columns and nl levels.  The value is a list with the starting
## value of theta for the block, the lower bounds, the block of
## Lambdat and the function that updates the block given the section
## of theta for this block.

##' Create digonal block on transposed relative covariance factor
##'
##' Each random-effects term is represented by diagonal block on
##' the transposed relative covariance factor. \code{blockLambdat}
##' creates such a block, and returns related information along
##' with it.
##'
##' @param nl Number of levels in a grouping factor for a particular
##' random effects term (the number of levels in the \code{grp} argument
##' in \code{\link{Zsection}}).
##' @param nc Number of columns in a dense model matrix for a particular
##' random effects term (the number of columns in the \code{mm} argument
##' in \code{\link{Zsection}}).
##' @return A \code{list} with:
##' \itemize{
##' \item the block
##' \item ititial values of theta for the block
##' \item lower bounds on these initial theta values
##' \item a function that updates the block given the section
##'   of theta for this block
##' }
##' @export
##' @examples
##' (l <- blockLambdat(2, 3))
##' within(l, slot(Lambdat, 'x') <- updateLambdatx(as.numeric(10:12)))
##' @rdname utilities
blockLambdat <- function(nl, nc) {
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
##' \code{\link{blockLambdat}} once for each random effects term (ie.
##' each list element in \code{grps} and \code{mms}). The results of
##' \code{\link{Zsection}} for each term are \code{rBind}ed together.
##' The results of \code{\link{blockLambdat}} are \code{bdiag}ed
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
##' @rdname utilities
mkRanefRepresentation <- function(grps, mms) {
    ## List to be returned starts with Zt, the transposed model matrix.
    ll <- list(Zt = do.call(rBind, mapply(Zsection, grps, mms)))

    nl <- sapply(grps, nlevels) # number of levels in each grouping factor
    nc <- sapply(mms, ncol) # number of columns in each model matrix

    if (all(nc == 1L)) {    # for scalar models Lambdat is diagonal ("ddiMatrix")
        nth <- length(nc)
        ll$lower <- numeric(nth)
        ll$theta <- rep.int(1, nth)
        ll$upper <- rep.int(Inf, nth)
        ll$Lambdat <- Diagonal(x=rep.int(1, sum(nl)))
        ll$thfun <- local({
            nlvs <- nl
            function(theta) rep.int(theta,nlvs)
        })
    } else { # otherwise Lambdat is block-diagonal in "dgCMatrix" format
        zz <- mapply(blockLambdat, nl, nc, SIMPLIFY=FALSE)
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
