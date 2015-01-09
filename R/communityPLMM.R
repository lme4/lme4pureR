##' Phylogenetic Generalized Linear Mixed Model for Community Data
##'
##' This function performs Generalized Linear Mixed Models for binary
##' and continuous phylogenetic data, estimating regression
##' coefficients with approximate standard errors. It is a modeled
##' after \code{lmer} but is more general by allowing correlation
##' structure within random effects; these correlations can be
##' phylogenetic among species, or any other correlation structure,
##' such as geographical correlations among sites. It is, however,
##' much more specific than \code{lmer} in that it can only analyze a
##' subset of the types of model designed handled by \code{lmer}. It
##' is also much slower than \code{lmer} and requires users to specify
##' correlation structures as covariance matrices.
##' \code{communityPGLMM} can analyze models in Ives and Helmus
##' (2011). It can also analyze bipartite phylogenetic data, such as
##' that analyzed in Rafferty and Ives (2011), by giving sites
##' phylogenetic correlations.
##'
##' @param formula a two-sided linear formula object describing the
##' fixed-effects of the model; for example, \code{Y ~ X}.
##' @param data data frame containing the variables named in formula. The
##' data frame should have long format with factors specifying species
##' and sites.  \code{communityPGLMM} will reorder rows of the data
##' frame so that species are nested within sites.
##' @param family either \code{gaussian} for a Linear Mixed Model, or
##' \code{binomial} for binary dependent data.
##' @param sp a factor variable that identifies species
##' @param site a factor variable that identifies sites
##' @param random.effects a list that contains, for non-nested random
##' effects, lists of triplets of the form \code{list(X, group =
##' group, covar = V)}. This is modeled after the \code{lmer} formula
##' syntax \code{(X | group)} where \code{X} is a variable and group
##' is a grouping factor. Note that group should be either your sp or
##' site variable specified in sp and site. The additional term
##' \code{V} is a covariance matrix of rank equal to the number of
##' levels of group that specifies the covariances among groups in the
##' random effect \code{X}. For nested variable random effects,
##' random.effects contains lists of quadruplets of the form
##' \code{list(X, group1 = group1, covar = V, group2 = group2)} where
##' \code{group1} is nested within \code{group2}.
##' @param REML whether REML or ML is used for model fitting. For the
##' generalized linear mixed model for binary data, these don't have
##' standard interpretations, and there is no log likelihood function
##' that can be used in likelihood ratio tests.
##' @param s2.init an array of initial estimates of \code{s2} for each
##' random effect that scales the variance. If s2.init is not provided
##' for \code{family="gaussian"}, these are estimated using in a
##' clunky way using \code{lm} assuming no phylogenetic signal.  A
##' better approach is to run \code{lmer} and use the output random
##' effects for \code{s2.init}. If \code{s2.init} is not provided for
##' \code{family="binomial"}, these are set to 0.25.
##' @param B.init initial estimates of B, a matrix containing
##' regression coefficients in the model for the fixed effects. This
##' matrix must have dim(B.init)=c(p+1,1), where p is the number of
##' predictor (independent) variables; the first element of B
##' corresponds to the intercept, and the remaining elements
##' correspond in order to the predictor (independent) variables in
##' the formula.  If B.init is not provided, these are estimated using
##' in a clunky way using lm() or glm() assuming no phylogenetic
##' signal.  A better approach is to run lmer() and use the output
##' fixed effects for B.init.
##' @param reltol a control parameter dictating the relative tolerance
##' for convergence in the optimization; see optim().
##' @param maxit a control parameter dictating the maximum number of
##' iterations in the optimization; see optim().
##' @param reltol.pgl a control parameter dictating the tolerance for
##' convergence in the PQL estimates of the mean components of the
##' binomial GLMM.
##' @param maxit.pgl a control parameter dictating the maximum number
##' of iterations in the PQL estimates of the mean components of the
##' binomial GLMM.
##' @param verbose if TRUE, the model deviance and running estimates
##' of s2 and B are plotted each iteration during optimization.
##'
##' @details For linear mixed models (family = "gaussian"), the
##' function estimates parameters for the model of the form FIXME
communityPGLMM <- function(formula, data = list(), family = "gaussian",
                           sp = NULL, site = NULL, random.effects = list(), REML =
                           TRUE, s2.init = NULL, B.init = NULL, reltol = 10^-6, maxit = 500,
                           tol.pql = 10^-6, maxit.pql = 200, verbose = FALSE) {

    if (family == "gaussian") 
        z <- communityPGLMM.gaussian(formula = formula, data = data, sp = sp,
                                     site = site, random.effects = random.effects, 
                                     REML = REML, s2.init = s2.init, B.init = B.init,
                                     reltol = reltol, maxit = maxit, verbose = verbose)
    if (family == "binomial") {
        s2.init <- 0.25
        z <- communityPGLMM.binary(formula = formula, data = data, sp = sp, site = site,
                                   random.effects = random.effects, REML = REML,
                                   s2.init = s2.init, B.init = B.init, reltol = reltol,
                                   maxit = maxit, tol.pql = tol.pql, 
                                   maxit.pql = maxit.pql, verbose = verbose)
    }
    if (!is.element(family, c("gaussian", "binomial"))) 
        cat("\nSorry, but only binomial (binary) and gaussian options exist at this time")
    return(z)
}

######################################################
######################################################
# communityPLMM.gaussian
######################################################
######################################################
communityPGLMM.gaussian <- function(formula, data = list(), family = "gaussian",
                                    sp = NULL, site = NULL, random.effects = list(),
                                    REML = TRUE, s2.init = NULL, B.init = NULL,
                                    reltol = 10^-8, maxit = 500, verbose = FALSE) {

                                        # Begin pglmm.LL
    plmm.LL <- function(par, X, Y, Zt, St, nestedsp = NULL, nestedsite = NULL, REML, verbose) {
        n <- dim(X)[1]
        p <- dim(X)[2]
        
        if (!is.null(St)) {
            q.nonNested <- dim(St)[1]
            sr <- Re(par[1:q.nonNested])
            iC <- sr[1] * St[1, ]
            if (length(sr) > 1) 
                for (i in 2:q.nonNested) {
                    iC <- iC + sr[i] * St[i, ]
                }
            iC <- as(diag(iC), "dsCMatrix")
            Ut <- iC %*% Zt
            U <- t(Ut)
        } else {
            q.nonNested <- 0
            sr <- NULL
        }
        if (is.null(nestedsp[[1]])) {
            q.Nested <- 0
        } else {
            q.Nested <- length(nestedsp)
        }
        
        if (q.Nested == 0) {
            sn <- NULL
        } else {
            sn <- Re(par[(q.nonNested + 1):(q.nonNested + q.Nested)])
        }
        
        if (q.Nested == 0) {
            iA <- as(diag(n), "dsCMatrix")
            Ishort <- as(diag(nrow(Ut)), "dsCMatrix")
            Ut.iA.U <- Ut %*% U
                                        # Woodbury identity
            iV <- iA - U %*% solve(Ishort + Ut.iA.U) %*% Ut
        } else {
            A <- as(diag(n), "dsCMatrix")
            for (j in 1:q.Nested) {
                A <- A + sn[j]^2 * kronecker(nestedsite[[j]], nestedsp[[j]])
            }
            iA <- solve(A)
            if (q.nonNested > 0) {
                Ishort <- as(diag(nrow(Ut)), "dsCMatrix")
                Ut.iA.U <- Ut %*% iA %*% U
                iV <- iA - iA %*% U %*% solve(Ishort + Ut.iA.U) %*% Ut %*% iA
            } else {
                iV <- iA
            }
        }
        
        denom <- t(X) %*% iV %*% X
        num <- t(X) %*% iV %*% Y
        B <- solve(denom, num)
        B <- as.matrix(B)
        H <- Y - X %*% B
        
        if (q.Nested == 0) {
                                        # Sylvester identity
            logdetV <- determinant(Ishort + Ut.iA.U)$modulus[1]
            if (is.infinite(logdetV)) 
                logdetV <- 2 * sum(log(diag(chol(Ishort + Ut.iA.U))))
        } else {
            logdetV <- -determinant(iV)$modulus[1]
            if (is.infinite(logdetV)) 
                logdetV <- -2 * sum(log(diag(chol(iV, pivot = T))))
            if (is.infinite(logdetV)) 
                return(10^10)
        }
        
        
        if (REML == TRUE) {
                                        # concentrated REML likelihood function
            s2.conc <- t(H) %*% iV %*% H/(n - p)
            LL <- 0.5 * ((n - p) * log(s2.conc) + logdetV + (n - p) + log(det(t(X) %*% iV %*% X)))
        } else {
                                        # concentrated ML likelihood function
            s2.conc <- t(H) %*% iV %*% H/n
            LL <- 0.5 * (n * log(s2.conc) + logdetV + n)
        }
        
        if (verbose == T) 
            show(c(as.numeric(LL), par))
        return(as.numeric(LL))
    }
                                        # End plmm.LL
    
                                        # Main program
    if (is.null(sp) | is.null(site)) 
        stop("Categorical variables for 'sp' and 'site' must be specified")
    nspp <- nlevels(sp)
    nsite <- nlevels(site)
    
                                        # order data first by site, second by species
    sp.order <- order(sp)
    data <- data[sp.order, ]
    sp <- sp[sp.order]
    site <- site[sp.order]

    site.order <- order(site)
    data <- data[site.order, ]
    sp <- sp[site.order]
    site <- site[site.order]
    
    mf <- model.frame(formula = formula, data = data)
    X <- model.matrix(attr(mf, "terms"), data = mf)
    Y <- model.response(mf)
    
    re <- random.effects
    q <- length(re)
    
    Ztt <- list(NULL)
    St.lengths <- array(0, q)
    nestedsp <- list(NULL)
    nestedsite <- list(NULL)
    ii <- 0
    jj <- 0
    
    for (i in 1:q) {
        re.i <- re[[i]]
                                        # non-nested terms
        if (length(re.i) == 3) {
            if (setequal(levels(re.i[[2]]), levels(sp)) && all(re.i[[2]] == sp)) {
                Zt.i <- kronecker(matrix(1, nrow = 1, ncol = nsite), chol(re.i[[3]]))
                if (length(re.i[[1]]) > 1) {
                    Zt.i <- Zt.i * kronecker(t(re.i[[1]]), matrix(1, nrow = nspp, ncol = 1))
                }
                ii <- ii + 1
                Ztt[[ii]] <- Zt.i
                St.lengths[ii] <- nspp
            }
            if (setequal(levels(re.i[[2]]), levels(site)) && all(re.i[[2]] == site)) {
                Zt.i <- kronecker(chol(re.i[[3]]), matrix(re.i[[1]], nrow = 1, ncol = nspp))
                if (length(re.i[[1]]) > 1) {
                    Zt.i <- Zt.i * kronecker(t(re.i[[1]]), matrix(1, nrow = nspp, ncol = 1))
                }
                ii <- ii + 1
                Ztt[[ii]] <- Zt.i
                St.lengths[ii] <- nsite
            }
        }
        
                                        # nested terms
        if (length(re.i) == 4) {
            if (setequal(levels(re.i[[2]]), levels(sp)) && all(re.i[[2]] == sp)) {
                if (length(re.i[[1]]) > 1) 
                    stop("Nested terms can only be for intercepts")
                nestedsp.j <- re.i[[3]]
                nestedsite.j <- diag(nsite)
            }
            if (setequal(levels(re.i[[2]]), levels(site)) && all(re.i[[2]] == site)) {
                if (length(re.i[[1]]) > 1) 
                    stop("Nested terms can only be for intercepts")
                nestedsp.j <- diag(nspp)
                nestedsite.j <- re.i[[3]]
            }
            jj <- jj + 1
            nestedsp[[jj]] <- nestedsp.j
            nestedsite[[jj]] <- nestedsite.j
        }
    }
    q.nonNested <- ii
    q.Nested <- jj
    
    if (q.nonNested > 0) {
        St <- matrix(0, nrow = q.nonNested, ncol = sum(St.lengths))
        Zt <- matrix(0, nrow = sum(St.lengths), ncol = nspp * nsite)
        count <- 1
        for (i in 1:q.nonNested) {
            St[i, count:(count + St.lengths[i] - 1)] <- matrix(1, nrow = 1, ncol = St.lengths[i])
            Zt[count:(count + St.lengths[i] - 1), ] <- Ztt[[i]]
            count <- count + St.lengths[i]
        }
        Zt <- as(Zt, "dgTMatrix")
        St <- as(St, "dgTMatrix")
    } else {
        Zt <- NULL
        St <- NULL
    }
    
    p <- ncol(X)
    n <- nrow(X)
    
                                        # Compute initial estimates
                                        # assuming no phylogeny if not
                                        # provided
    if (!is.null(B.init) & length(B.init) != p) {
        warning("B.init not correct length, so computed B.init using glm()")
    }
    if ((is.null(B.init) | (!is.null(B.init) & length(B.init) != p)) & !is.null(s2.init)) {
        B.init <- t(matrix(lm(formula = formula, data = data)$coefficients, ncol = p))
    }
    if (!is.null(B.init) & is.null(s2.init)) {
        s2.init <- var(lm(formula = formula, data = data)$residuals)/q
    }
    if ((is.null(B.init) | (!is.null(B.init) & length(B.init) != p)) & is.null(s2.init)) {
        B.init <- t(matrix(lm(formula = formula, data = data)$coefficients, ncol = p))
        s2.init <- var(lm(formula = formula, data = data)$residuals)/q
    }
    B <- B.init
    s <- as.vector(array(s2.init^0.5, dim = c(1, q)))
    
    if (q > 1) {
        opt <- optim(fn = plmm.LL, par = s, X = X, Y = Y, Zt = Zt, St = St,
                     nestedsp = nestedsp, nestedsite = nestedsite, 
                     REML = REML, verbose = verbose, method = "Nelder-Mead",
                     control = list(maxit = maxit, reltol = reltol))
    } else {
        opt <- optim(fn = plmm.LL, par = s, X = X, Y = Y, Zt = Zt, St = St,
                     nestedsp = nestedsp, nestedsite = nestedsite, 
                     REML = REML, verbose = verbose, method = "L-BFGS-B",
                     control = list(maxit = maxit))

    }
                                        # Extract parameters
    par <- abs(Re(opt$par))
    LL <- opt$value
    if (!is.null(St)) {
        q.nonNested <- dim(St)[1]
        sr <- Re(par[1:q.nonNested])
        iC <- sr[1] * St[1, ]
        if (length(sr) > 1) 
            for (i in 2:q.nonNested) {
                iC <- iC + sr[i] * St[i, ]
            }
        iC <- as(diag(iC), "dsCMatrix")
        Ut <- iC %*% Zt
        U <- t(Ut)
    } else {
        q.nonNested <- 0
        sr <- NULL
    }
    if (is.null(nestedsp[[1]])) {
        q.Nested <- 0
    } else {
        q.Nested <- length(nestedsp)
    }
    if (q.Nested == 0) {
        sn <- NULL
    } else {
        sn <- Re(par[(q.nonNested + 1):(q.nonNested + q.Nested)])
    }
    if (q.Nested == 0) {
        iA <- as(diag(n), "dsCMatrix")
        Ishort <- as(diag(nrow(Ut)), "dsCMatrix")
        Ut.iA.U <- Ut %*% U
                                        # Woodbury identity
        iV <- iA - U %*% solve(Ishort + Ut.iA.U) %*% Ut
    } else {
        A <- as(diag(n), "dsCMatrix")
        for (j in 1:q.Nested) {
            A <- A + sn[j]^2 * kronecker(nestedsite[[j]], nestedsp[[j]])
        }
        iA <- solve(A)
        if (q.nonNested > 0) {
            Ishort <- as(diag(nrow(Ut)), "dsCMatrix")
            Ut.iA.U <- Ut %*% iA %*% U
            iV <- iA - iA %*% U %*% solve(Ishort + Ut.iA.U) %*% Ut %*% iA
        } else {
            iV <- iA
        }
    }
    
    denom <- t(X) %*% iV %*% X
    num <- t(X) %*% iV %*% Y
    B <- solve(denom, num)
    B <- as.matrix(B)
    H <- Y - X %*% B
    
    if (q.Nested == 0) {
                                        # Sylvester identity
        logdetV <- determinant(Ishort + Ut.iA.U)$modulus[1]
        if (is.infinite(logdetV)) {
            logdetV <- 2 * sum(log(diag(chol(Ishort + Ut.iA.U))))
        }
    } else {
        logdetV <- -determinant(iV)$modulus[1]
        if (is.infinite(logdetV)) {
            logdetV <- -2 * sum(log(diag(chol(iV, pivot = T))))
        }
        if (is.infinite(logdetV)) {
            return(10^10)
	}
    }

    if (REML == TRUE) {
        s2resid <- as.numeric(t(H) %*% iV %*% H/(n - p))
    } else {
        s2resid <- as.numeric(t(H) %*% iV %*% H/n)
    }
        
    s2r <- s2resid * sr^2
    s2n <- s2resid * sn^2
    ss <- c(sr, sn, s2resid^0.5)
    
    iV <- iV/s2resid
    
    B.cov <- solve(t(X) %*% iV %*% X)
    B.se <- as.matrix(diag(B.cov))^0.5
    B.zscore <- B/B.se
    B.pvalue <- 2 * pnorm(abs(B/B.se), lower.tail = FALSE)

    if (REML == TRUE) {
        logLik <- -0.5 * (n - p) * log(2 * pi) + 0.5 * log(det(t(X) %*% X)) - LL
    } else {
        logLik <- -0.5 * n * log(2 * pi) - LL
    }
    k <- p + q + 1
    AIC <- -2 * logLik + 2 * k
    BIC <- -2 * logLik + k * (log(n) - log(pi))
        
    results <- list(formula = formula, data = data, family = family,
                    random.effects = random.effects, B = B, 
                    B.se = B.se, B.cov = B.cov, B.zscore = B.zscore,
                    B.pvalue = B.pvalue, ss = ss, s2r = s2r, s2n = s2n, 
                    s2resid = s2resid, logLik = logLik, AIC = AIC, BIC = BIC,
                    REML = REML, s2.init = s2.init, B.init = B.init, 
                    Y = Y, X = X, H = H, iV = iV, mu = NULL, nestedsp = nestedsp,
                    nestedsite = nestedsite, sp = sp, site = site, 
                    Zt = Zt, St = St, convcode = opt$convergence,
                    niter = opt$counts)

    class(results) <- "communityPGLMM"
    results
}

######################################################
######################################################
# communityPGLMM.binary
######################################################
######################################################
communityPGLMM.binary <- function(formula, data = list(), family = "binomial",
                                  sp = NULL, site = NULL, random.effects = list(), 
                                  REML = TRUE, s2.init = 0.25, B.init = NULL,
                                  reltol = 10^-5, maxit = 40, tol.pql = 10^-6,
                                  maxit.pql = 200, verbose = FALSE) {

    plmm.binary.V <- function(par, Zt, St, mu, nestedsp = NULL, nestedsite = NULL) {
            
        if (!is.null(St)) {
            q.nonNested <- dim(St)[1]
            sr <- Re(par[1:q.nonNested])
            iC <- sr[1] * St[1, ]
            if (length(sr) > 1) 
                for (i in 2:q.nonNested) {
                    iC <- iC + sr[i] * St[i, ]
                }
            iC <- as(diag(iC), "dsCMatrix")
            Ut <- iC %*% Zt
            U <- t(Ut)
        } else {
            q.nonNested <- 0
            sr <- NULL
        }
        if (is.null(nestedsp[[1]])) {
            q.Nested <- 0
        } else {
            q.Nested <- length(nestedsp)
        }
        
        if (q.Nested == 0) {
            sn <- NULL
        } else {
            sn <- Re(par[(q.nonNested + 1):(q.nonNested + q.Nested)])
        }
        
        iW <- diag(as.vector((mu * (1 - mu))^-1))
        if (q.Nested == 0) {
            A <- iW
        } else {
            A <- iW
            for (j in 1:q.Nested) {
                A <- A + sn[j]^2 * kronecker(nestedsite[[j]], nestedsp[[j]])
            }
        }
        if (q.nonNested > 0) {
            V <- A + U %*% Ut
        } else {
            V <- A
        }
        return(V)
    }
                                        # End plmm.binary.V
    
    plmm.binary.iV <- function(par, Zt, St, mu, nestedsp = NULL, nestedsite = NULL) {
        
        if (!is.null(St)) {
            q.nonNested <- dim(St)[1]
            sr <- Re(par[1:q.nonNested])
            iC <- sr[1] * St[1, ]
            if (length(sr) > 1) { 
                for (i in 2:q.nonNested) {
                    iC <- iC + sr[i] * St[i, ]
                }
            }
            iC <- as(diag(iC), "dsCMatrix")
            Ut <- iC %*% Zt
            U <- t(Ut)
        } else {
            q.nonNested <- 0
            sr <- NULL
        }
        if (is.null(nestedsp[[1]])) {
            q.Nested <- 0
        } else {
            q.Nested <- length(nestedsp)
        }
        
        if (q.Nested == 0) {
            sn <- NULL
        } else {
            sn <- Re(par[(q.nonNested + 1):(q.nonNested + q.Nested)])
        }
        
        if (q.Nested == 0) {
            iA <- diag(as.vector((mu * (1 - mu))))
            Ishort <- as(diag(nrow(Ut)), "dsCMatrix")
            Ut.iA.U <- Ut %*% iA %*% U
                                        # Woodbury identity
            iV <- iA - iA %*% U %*% solve(Ishort + Ut.iA.U) %*% Ut %*% iA
        } else {
            A <- diag(as.vector((mu * (1 - mu))^-1))
            for (j in 1:q.Nested) {
                A <- A + sn[j]^2 * kronecker(nestedsite[[j]], nestedsp[[j]])
            }
            iA <- solve(A)
            if (q.nonNested > 0) {
                Ishort <- as(diag(nrow(Ut)), "dsCMatrix")
                Ut.iA.U <- Ut %*% iA %*% U
                iV <- iA - iA %*% U %*% solve(Ishort + Ut.iA.U) %*% Ut %*% iA
            } else {
                iV <- iA
            }
        }
        return(iV)
    }
                                        # End plmm.binary.iV
    
    plmm.binary.logdetV <- function(par, Zt, St, mu, nestedsp = NULL, nestedsite = NULL) {
        n <- dim(X)[1]
        p <- dim(X)[2]
        
        if (!is.null(St)) {
            q.nonNested <- dim(St)[1]
            sr <- Re(par[1:q.nonNested])
            iC <- sr[1] * St[1, ]
            if (length(sr) > 1) 
                for (i in 2:q.nonNested) {
                    iC <- iC + sr[i] * St[i, ]
                }
            iC <- as(diag(iC), "dsCMatrix")
            Ut <- iC %*% Zt
            U <- t(Ut)
        } else {
            q.nonNested <- 0
            sr <- NULL
        }
        if (is.null(nestedsp[[1]])) {
            q.Nested <- 0
        } else {
            q.Nested <- length(nestedsp)
        }
        
        if (q.Nested == 0) {
            sn <- NULL
        } else {
            sn <- Re(par[(q.nonNested + 1):(q.nonNested + q.Nested)])
        }
        
        if (q.Nested == 0) {
            iA <- diag(as.vector((mu * (1 - mu))))
            Ishort <- as(diag(nrow(Ut)), "dsCMatrix")
            Ut.iA.U <- Ut %*% iA %*% U
            logdetV <- determinant(Ishort + Ut.iA.U)$modulus[1] - determinant(iA)$modulus[1]
            if (is.infinite(logdetV)) 
                logdetV <- 2 * sum(log(diag(chol(Ishort + Ut.iA.U)))) - determinant(iA)$modulus[1]
        } else {
            A <- diag(as.vector((mu * (1 - mu))^-1))
            for (j in 1:q.Nested) {
                A <- A + sn[j]^2 * kronecker(nestedsite[[j]], nestedsp[[j]])
            }
            iA <- solve(A)
            if (q.nonNested > 0) {
                Ishort <- as(diag(nrow(Ut)), "dsCMatrix")
                Ut.iA.U <- Ut %*% iA %*% U
                iV <- iA - iA %*% U %*% solve(Ishort + Ut.iA.U) %*% Ut %*% iA
            } else {
                iV <- iA
            }
            logdetV <- -determinant(iV)$modulus[1]
            if (is.infinite(logdetV)) 
                logdetV <- -2 * sum(log(diag(chol(iV, pivot = T))))
            if (is.infinite(logdetV)) 
                return(list(iV = NULL, logdetV = NULL))
        }
        
        return(logdetV)
    }
                                        # End plmm.binary.logdetV
	
                                        # Begin pglmm.binary.LL
    plmm.binary.LL <- function(par, H, X, Zt, St, mu, nestedsp = NULL, nestedsite = NULL,
                               REML = TRUE, verbose = FALSE) {
        par <- abs(par)
        n <- dim(H)[1]
        p <- dim(H)[2]
        
        iV <- plmm.binary.iV(par = par, Zt = Zt, St = St, mu = mu,
                             nestedsp = nestedsp, nestedsite = nestedsite)
        logdetV <- plmm.binary.logdetV(par = par, Zt = Zt, St = St, mu = mu,
                                       nestedsp = nestedsp, nestedsite = nestedsite)
        if (REML == TRUE) {
                                        # REML likelihood function
            LL <- 0.5 * (logdetV + t(H) %*% iV %*% H + log(det(t(X) %*% iV %*% X)))
        } else {
                                        # ML likelihood function
            LL <- 0.5 * (logdetV + t(H) %*% iV %*% H)
        }
        if (verbose == T) 
            show(c(as.numeric(LL), par))
        
        return(as.numeric(LL))
    }
                                        # End plmm.binary.LL
    

                                        # Begin main program
    if (is.null(sp) | is.null(site)) 
        stop("Categorical variables for 'sp' and 'site' must be specified")
    nspp <- nlevels(sp)
    nsite <- nlevels(site)
    
                                        # order data first by site, second by species
    sp.order <- order(sp)
    data <- data[sp.order, ]
    sp <- sp[sp.order]
    site <- site[sp.order]
    
    site.order <- order(site)
    data <- data[site.order, ]
    sp <- sp[site.order]
    site <- site[site.order]
    
    mf <- model.frame(formula = formula, data = data)
    X <- model.matrix(attr(mf, "terms"), data = mf)
    Y <- model.response(mf)
    
    re <- random.effects
    q <- length(re)
    
    Ztt <- list(NULL)
    St.lengths <- array(0, q)
    nestedsp <- list(NULL)
    nestedsite <- list(NULL)
    ii <- 0
    jj <- 0
    
    for (i in 1:q) {
        re.i <- re[[i]]
                                        # non-nested terms
        if (length(re.i) == 3) {
            if (setequal(levels(re.i[[2]]), levels(sp)) && all(re.i[[2]] == sp)) {
                Zt.i <- kronecker(matrix(1, nrow = 1, ncol = nsite), chol(re.i[[3]]))
                if (length(re.i[[1]]) > 1) {
                    Zt.i <- Zt.i * kronecker(t(re.i[[1]]), matrix(1, nrow = nspp, ncol = 1))
                }
                ii <- ii + 1
                Ztt[[ii]] <- Zt.i
                St.lengths[ii] <- nspp
            }
            if (setequal(levels(re.i[[2]]), levels(site)) && all(re.i[[2]] == site)) {
                Zt.i <- kronecker(chol(re.i[[3]]), matrix(re.i[[1]], nrow = 1, ncol = nspp))
                if (length(re.i[[1]]) > 1) {
                    Zt.i <- Zt.i * kronecker(t(re.i[[1]]), matrix(1, nrow = nspp, ncol = 1))
                }
                ii <- ii + 1
                Ztt[[ii]] <- Zt.i
                St.lengths[ii] <- nsite
            }
        }
        
                                        # nested terms
        if (length(re.i) == 4) {
            if (setequal(levels(re.i[[2]]), levels(sp)) && all(re.i[[2]] == sp)) {
                if (length(re.i[[1]]) > 1) 
                    stop("Nested terms can only be for intercepts")
                nestedsp.j <- re.i[[3]]
                nestedsite.j <- diag(nsite)
            }
            if (setequal(levels(re.i[[2]]), levels(site)) && all(re.i[[2]] == site)) {
                if (length(re.i[[1]]) > 1) 
                    stop("Nested terms can only be for intercepts")
                nestedsp.j <- diag(nspp)
                nestedsite.j <- re.i[[3]]
            }
            jj <- jj + 1
            nestedsp[[jj]] <- nestedsp.j
            nestedsite[[jj]] <- nestedsite.j
        }
    }
    q.nonNested <- ii
    q.Nested <- jj
    
    if (q.nonNested > 0) {
        St <- matrix(0, nrow = q.nonNested, ncol = sum(St.lengths))
        Zt <- matrix(0, nrow = sum(St.lengths), ncol = nspp * nsite)
        count <- 1
        for (i in 1:q.nonNested) {
            St[i, count:(count + St.lengths[i] - 1)] <- matrix(1, nrow = 1, ncol = St.lengths[i])
            Zt[count:(count + St.lengths[i] - 1), ] <- Ztt[[i]]
            count <- count + St.lengths[i]
        }
        Zt <- as(Zt, "dgTMatrix")
        St <- as(St, "dgTMatrix")
    } else {
        Zt <- NULL
        St <- NULL
    }
    
    p <- ncol(X)
    n <- nrow(X)
    
                                        # Compute initial estimates
                                        # assuming no phylogeny if not
                                        # provided
    if (!is.null(B.init) & length(B.init) != p) {
        warning("B.init not correct length, so computed B.init using glm()")
    }
    if ((is.null(B.init) | (!is.null(B.init) & length(B.init) != p))) {
        B.init <- t(matrix(glm(formula = formula, data = data,
                               family = binomial)$coefficients, ncol = p))
    } else {
        B.init <- matrix(B.init, ncol = 1)
    }
    B <- B.init
    ss <- as.vector(array(s2.init^0.5, dim = c(1, q)))
    
    b <- matrix(0, nrow = n)
    beta <- rbind(B, b)
    mu <- exp(X %*% B)/(1 + exp(X %*% B))
    XX <- cbind(X, diag(1, nrow = n, ncol = n))
    
    est.ss <- ss
    est.B <- B
    oldest.ss <- 10^6
    oldest.B <- matrix(10^6, nrow = length(est.B))
    
    iteration <- 0
    exitflag <- 0
    rcondflag <- 0
    while (((t(est.ss - oldest.ss) %*% (est.ss - oldest.ss) > tol.pql^2) |
            (t(est.B - oldest.B) %*% (est.B - oldest.B) > tol.pql^2)) &
           (iteration <= maxit.pql)) {

        iteration <- iteration + 1
        oldest.ss <- est.ss
        oldest.B <- est.B
        
        est.B.m <- B
        oldest.B.m <- matrix(10^6, nrow = length(est.B))
        
                                        # mean component
        while ((t(est.B.m - oldest.B.m) %*% (est.B.m - oldest.B.m) > tol.pql^2) &
               (iteration <= maxit.pql)) {
            
            oldest.B.m <- est.B.m
            
            iV <- plmm.binary.iV(par = ss, Zt = Zt, St = St, mu = mu,
                                 nestedsp = nestedsp, nestedsite = nestedsite)
            
            Z <- X %*% B + b + (Y - mu)/(mu * (1 - mu))
            denom <- t(X) %*% iV %*% X
            num <- t(X) %*% iV %*% Z
            B <- solve(denom, num)
            B <- as.matrix(B)
            
            V <- plmm.binary.V(par = ss, Zt = Zt, St = St, mu = mu,
                               nestedsp = nestedsp, nestedsite = nestedsite)
            iW <- diag(as.vector((mu * (1 - mu))^-1))
            C <- V - iW
            b <- C %*% iV %*% (Z - X %*% B)
            beta <- rbind(B, matrix(b))
            mu <- exp(XX %*% beta)/(1 + exp(XX %*% beta))
            
            est.B.m <- B
            if (verbose == TRUE) 
                show(c(iteration, B))
        }
        
                                        # variance component
        Z <- X %*% B + b + (Y - mu)/(mu * (1 - mu))
        H <- Z - X %*% B
        if (q > 1) {
            opt <- optim(fn = plmm.binary.LL, par = ss, H = H, X = X, Zt = Zt,
                         St = St, mu = mu, nestedsp = nestedsp, 
                         nestedsite = nestedsite, REML = REML, verbose = verbose,
                         method = "Nelder-Mead",
                         control = list(maxit = maxit, reltol = reltol))
        } else {
            opt <- optim(fn = plmm.binary.LL, par = ss, H = H, X = X, Zt = Zt,
                         St = St, mu = mu, nestedsp = nestedsp, nestedsite = nestedsite,
                         REML = REML, verbose = verbose, method = "L-BFGS-B",
                         control = list(maxit = maxit))
        }
        ss <- abs(opt$par)
        LL <- opt$value
        
        est.ss <- ss
        est.B <- B
    }
    
                                        # Extract parameters
    if (q.nonNested > 0) {
		sr <- ss[1:q.nonNested]
            } else {
		sr <- NULL
            }
    if (q.Nested > 0) {
        sn <- ss[(q.nonNested + 1):(q.nonNested + q.Nested)]
    } else {
        sn <- NULL
    }
    
    s2r <- sr^2
    s2n <- sn^2
    
    B.cov <- solve(t(X) %*% iV %*% X)
    B.se <- as.matrix(diag(B.cov))^0.5
    B.zscore <- B/B.se
    B.pvalue <- 2 * pnorm(abs(B/B.se), lower.tail = FALSE)
    
    results <- list(formula = formula, data = data, family = family,
                    random.effects = random.effects, B = B, 
                    B.se = B.se, B.cov = B.cov, B.zscore = B.zscore,
                    B.pvalue = B.pvalue, ss = ss, s2r = s2r, s2n = s2n, 
                    s2resid = NULL, logLik = NULL, AIC = NULL, BIC = NULL,
                    REML = REML, s2.init = s2.init, B.init = B.init, 
                    Y = Y, X = X, H = H, iV = iV, mu = mu, nestedsp = nestedsp,
                    nestedsite = nestedsite, sp = sp, site = site, 
                    Zt = Zt, St = St, convcode = opt$convergence, niter = opt$counts)
    class(results) <- "communityPGLMM"
    results
    return(results)
}

######################################################
######################################################
# communityPGLMM.binary.LRT
######################################################
######################################################
communityPGLMM.binary.LRT <- function(x, re.number = 0, ...) {

    plmm.binary.iV <- function(par, Zt, St, mu, nestedsp = NULL, nestedsite = NULL) {
        
        if (!is.null(St)) {
            q.nonNested <- dim(St)[1]
            sr <- Re(par[1:q.nonNested])
            iC <- sr[1] * St[1, ]
            if (length(sr) > 1) 
                for (i in 2:q.nonNested) {
                    iC <- iC + sr[i] * St[i, ]
                }
            iC <- as(diag(iC), "dsCMatrix")
			Ut <- iC %*% Zt
			U <- t(Ut)
        } else {
            q.nonNested <- 0
            sr <- NULL
        }
		if (is.null(nestedsp[[1]])) {
			q.Nested <- 0
                    } else {
			q.Nested <- length(nestedsp)
                    }
        
        if (q.Nested == 0) {
            sn <- NULL
        } else {
            sn <- Re(par[(q.nonNested + 1):(q.nonNested + q.Nested)])
        }
        
        if (q.Nested == 0) {
            iA <- diag(as.vector((mu * (1 - mu))))
            Ishort <- as(diag(nrow(Ut)), "dsCMatrix")
            Ut.iA.U <- Ut %*% iA %*% U
                                        # Woodbury identity
            iV <- iA - iA %*% U %*% solve(Ishort + Ut.iA.U) %*% Ut %*% iA
        } else {
            A <- diag(as.vector((mu * (1 - mu))^-1))
            for (j in 1:q.Nested) {
                A <- A + sn[j]^2 * kronecker(nestedsite[[j]], nestedsp[[j]])
            }
            iA <- solve(A)
            if (q.nonNested > 0) {
                Ishort <- as(diag(nrow(Ut)), "dsCMatrix")
                Ut.iA.U <- Ut %*% iA %*% U
                iV <- iA - iA %*% U %*% solve(Ishort + Ut.iA.U) %*% Ut %*% iA
            } else {
                iV <- iA
            }
        }
        return(iV)
    }
                                        # End plmm.binary.iV
    
    plmm.binary.logdetV <- function(par, Zt, St, mu, nestedsp = NULL, nestedsite = NULL) {
        n <- dim(X)[1]
        p <- dim(X)[2]
        
        if (!is.null(St)) {
            q.nonNested <- dim(St)[1]
            sr <- Re(par[1:q.nonNested])
            iC <- sr[1] * St[1, ]
            if (length(sr) > 1) 
                for (i in 2:q.nonNested) {
                    iC <- iC + sr[i] * St[i, ]
                }
            iC <- as(diag(iC), "dsCMatrix")
            Ut <- iC %*% Zt
            U <- t(Ut)
        } else {
            q.nonNested <- 0
            sr <- NULL
        }
        if (is.null(nestedsp[[1]])) {
            q.Nested <- 0
        } else {
            q.Nested <- length(nestedsp)
        }
        
        if (q.Nested == 0) {
            sn <- NULL
        } else {
            sn <- Re(par[(q.nonNested + 1):(q.nonNested + q.Nested)])
        }
        
        if (q.Nested == 0) {
            iA <- diag(as.vector((mu * (1 - mu))))
            Ishort <- as(diag(nrow(Ut)), "dsCMatrix")
            Ut.iA.U <- Ut %*% iA %*% U
            logdetV <- determinant(Ishort + Ut.iA.U)$modulus[1] - determinant(iA)$modulus[1]
            if (is.infinite(logdetV)) 
                logdetV <- 2 * sum(log(diag(chol(Ishort + Ut.iA.U)))) - determinant(iA)$modulus[1]
        } else {
            A <- diag(as.vector((mu * (1 - mu))^-1))
            for (j in 1:q.Nested) {
                A <- A + sn[j]^2 * kronecker(nestedsite[[j]], nestedsp[[j]])
            }
            iA <- solve(A)
            if (q.nonNested > 0) {
                Ishort <- as(diag(nrow(Ut)), "dsCMatrix")
                Ut.iA.U <- Ut %*% iA %*% U
                iV <- iA - iA %*% U %*% solve(Ishort + Ut.iA.U) %*% Ut %*% iA
            } else {
                iV <- iA
            }
            logdetV <- -determinant(iV)$modulus[1]
            if (is.infinite(logdetV)) 
                logdetV <- -2 * sum(log(diag(chol(iV, pivot = T))))
            if (is.infinite(logdetV)) 
                return(list(iV = NULL, logdetV = NULL))
        }
        
        return(logdetV)
    }
                                        # End plmm.binary.logdetV
    
                                        # Begin pglmm.binary.LL
    plmm.binary.LL <- function(par, H, X, Zt, St, mu, nestedsp = NULL,
                               nestedsite = NULL, REML = TRUE, verbose = FALSE) {
        par <- abs(par)
        n <- dim(H)[1]
        p <- dim(H)[2]
        
        iV <- plmm.binary.iV(par = par, Zt = Zt, St = St, mu = mu,
                             nestedsp = nestedsp, nestedsite = nestedsite)
        logdetV <- plmm.binary.logdetV(par = par, Zt = Zt, St = St, mu = mu,
                                       nestedsp = nestedsp, nestedsite = nestedsite)
        if (REML == TRUE) {
                                        # REML likelihood function
            LL <- 0.5 * (logdetV + t(H) %*% iV %*% H + log(det(t(X) %*% iV %*% X)))
        } else {
                                        # ML likelihood function
            s2.conc <- t(H) %*% iV %*% H/(n - p)
            LL <- 0.5 * (logdetV + t(H) %*% iV %*% H)
        }
        if (verbose == T) 
            show(c(as.numeric(LL), par))
        
        return(as.numeric(LL))
    }
                                        # End plmm.binary.LL
    

                                        # Begin main program
    
    n <- dim(x$X)[1]
    p <- dim(x$X)[2]
    par <- x$ss
    par[re.number] <- 0
    df <- length(re.number)
    
    LL <- plmm.binary.LL(par = x$ss, H = x$H, X = x$X, Zt = x$Zt, St = x$St,
                         mu = x$mu, nestedsp = x$nestedsp, 
                         nestedsite = x$nestedsite, REML = x$REML)
    if (x$REML == TRUE) {
        logLik <- -0.5 * (n - p - 1) * log(2 * pi) + 0.5 * log(det(t(x$X) %*% x$X)) - LL
    } else {
        logLik <- -0.5 * n * log(2 * pi) - LL
    }
    
    LL0 <- plmm.binary.LL(par = par, H = x$H, X = x$X, Zt = x$Zt, St = x$St,
                          mu = x$mu, nestedsp = x$nestedsp, 
                          nestedsite = x$nestedsite, REML = x$REML)
    if (x$REML == TRUE) {
        logLik0 <- -0.5 * (n - p - 1) * log(2 * pi) + 0.5 * log(det(t(x$X) %*% x$X)) - LL0
    } else {
        logLik0 <- -0.5 * n * log(2 * pi) - LL0
    }
    
    P.H0.s2 <- pchisq(2 * (logLik - logLik0), df = df, lower.tail = F)/2
    return(list(LR = logLik - logLik0, df = df, Pr = P.H0.s2))
}

######################################################
######################################################
# communityPGLMM.matrix.structure
######################################################
######################################################
communityPGLMM.matrix.structure <- function(formula, data = list(), family = "binomial",
                                            sp = NULL, site = NULL, random.effects = list(),
                                            ss = 1) {

    plmm.binary.V.test <- function(par, Zt, St, X, nestedsp = NULL, nestedsite = NULL) {
        n <- nrow(X)
        
        if (!is.null(St)) {
            q.nonNested <- dim(St)[1]
            sr <- Re(par[1:q.nonNested])
            iC <- sr[1] * St[1, ]
            if (length(sr) > 1) 
                for (i in 2:q.nonNested) {
                    iC <- iC + sr[i] * St[i, ]
                }
            iC <- as(diag(iC), "dsCMatrix")
            Ut <- iC %*% Zt
            U <- t(Ut)
        } else {
            q.nonNested <- 0
            sr <- NULL
        }
        if (is.null(nestedsp[[1]])) {
            q.Nested <- 0
        } else {
            q.Nested <- length(nestedsp)
        }
        
        if (q.Nested == 0) {
            sn <- NULL
        } else {
            sn <- Re(par[(q.nonNested + 1):(q.nonNested + q.Nested)])
        }
        
        iW <- 0 * diag(n)
        if (q.Nested == 0) {
            A <- iW
        } else {
            A <- iW
            for (j in 1:q.Nested) {
                A <- A + sn[j]^2 * kronecker(nestedsite[[j]], nestedsp[[j]])
            }
        }
        if (q.nonNested > 0) {
            V <- A + U %*% Ut
        } else {
            V <- A
        }
        return(V)
    }
                                        # End plmm.binary.V
    

                                        # Begin main program
    if (is.null(sp) | is.null(site)) 
        stop("Categorical variables for 'sp' and 'site' must be specified")
    nspp <- nlevels(sp)
    nsite <- nlevels(site)
    
                                        # order data first by site,
                                        # second by species
    sp.order <- order(sp)
    data <- data[sp.order, ]
    sp <- sp[sp.order]
    site <- site[sp.order]
    
    site.order <- order(site)
    data <- data[site.order, ]
    sp <- sp[site.order]
    site <- site[site.order]
    
    mf <- model.frame(formula = formula, data = data)
    X <- model.matrix(attr(mf, "terms"), data = mf)
    Y <- model.response(mf)
    
    re <- random.effects
    q <- length(re)
    
    Ztt <- list(NULL)
    St.lengths <- array(0, q)
    nestedsp <- list(NULL)
    nestedsite <- list(NULL)
    ii <- 0
    jj <- 0
    
    for (i in 1:q) {
        re.i <- re[[i]]
                                        # non-nested terms
        if (length(re.i) == 3) {
            if (setequal(levels(re.i[[2]]), levels(sp)) && all(re.i[[2]] == sp)) {
                Zt.i <- kronecker(matrix(1, nrow = 1, ncol = nsite), chol(re.i[[3]]))
                if (length(re.i[[1]]) > 1) {
                    Zt.i <- Zt.i * kronecker(t(re.i[[1]]), matrix(1, nrow = nspp, ncol = 1))
                }
                ii <- ii + 1
                Ztt[[ii]] <- Zt.i
                St.lengths[ii] <- nspp
            }
            if (setequal(levels(re.i[[2]]), levels(site)) && all(re.i[[2]] == site)) {
                Zt.i <- kronecker(chol(re.i[[3]]), matrix(re.i[[1]], nrow = 1, ncol = nspp))
                if (length(re.i[[1]]) > 1) {
                    Zt.i <- Zt.i * kronecker(t(re.i[[1]]), matrix(1, nrow = nspp, ncol = 1))
                }
                ii <- ii + 1
                Ztt[[ii]] <- Zt.i
                St.lengths[ii] <- nsite
            }
        }
        
                                        # nested terms
        if (length(re.i) == 4) {
            if (setequal(levels(re.i[[2]]), levels(sp)) && all(re.i[[2]] == sp)) {
                if (length(re.i[[1]]) > 1) 
                    stop("Nested terms can only be for intercepts")
                nestedsp.j <- re.i[[3]]
                nestedsite.j <- diag(nsite)
            }
            if (setequal(levels(re.i[[2]]), levels(site)) && all(re.i[[2]] == site)) {
                if (length(re.i[[1]]) > 1) 
                    stop("Nested terms can only be for intercepts")
                nestedsp.j <- diag(nspp)
                nestedsite.j <- re.i[[3]]
            }
            jj <- jj + 1
            nestedsp[[jj]] <- nestedsp.j
            nestedsite[[jj]] <- nestedsite.j
        }
    }
    q.nonNested <- ii
    q.Nested <- jj
    
    if (q.nonNested > 0) {
        St <- matrix(0, nrow = q.nonNested, ncol = sum(St.lengths))
        Zt <- matrix(0, nrow = sum(St.lengths), ncol = nspp * nsite)
        count <- 1
        for (i in 1:q.nonNested) {
            St[i, count:(count + St.lengths[i] - 1)] <- matrix(1, nrow = 1, ncol = St.lengths[i])
            Zt[count:(count + St.lengths[i] - 1), ] <- Ztt[[i]]
            count <- count + St.lengths[i]
        }
        Zt <- as(Zt, "dgTMatrix")
        St <- as(St, "dgTMatrix")
    } else {
        Zt <- NULL
        St <- NULL
    }
    
    V <- plmm.binary.V.test(par = array(ss, c(1, q)), Zt = Zt, St = St, X = X,
                            nestedsp = nestedsp, nestedsite = nestedsite)
    return(V)
}


######################################################
######################################################
# summary.communityPGLMM
######################################################
######################################################
summary.communityPGLMM <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    if (x$family == "gaussian") {
        if (x$REML == TRUE) 
            cat("Linear mixed model fit by restricted maximum likelihood")
        else cat("Linear mixed model fit by maximum likelihood")
    }
    if (x$family == "binomial") {
        if (x$REML == TRUE) 
            cat("Generalized linear mixed model for binary data fit by restricted maximum likelihood")
        else cat("Generalized linear mixed model for binary data fit by maximum likelihood")
    }
    
    cat("\n\nCall:")
    print(x$formula)
    cat("\n")
    
    if (x$family == "gaussian") {
        
        logLik = x$logLik
        AIC = x$AIC
        BIC = x$BIC
        
        names(logLik) = "logLik"
        names(AIC) = "AIC"
        names(BIC) = "BIC"
        print(c(logLik, AIC, BIC), digits = digits)
    }
    cat("\nRandom effects:\n")
    w <- data.frame(Variance = matrix(c(x$s2r, x$s2n, x$s2resid), ncol = 1),
                    Std.Dev = matrix(c(x$s2r^0.5, x$s2n^0.5, x$s2resid^0.5), ncol = 1))
    
    re.names <- NULL
    if (length(x$s2r) > 0) 
        for (i in 1:length(x$s2r)) re.names <- c(re.names, paste("non-nested ", i, sep = ""))
    
    if (length(x$s2n) > 0) 
        for (i in 1:length(x$s2n)) re.names <- c(re.names, paste("nested ", i, sep = ""))
    
    if (x$family == "gaussian") 
        re.names <- c(re.names, "residual")
    
    row.names(w) <- re.names
    print(w, digits = digits)
    
    cat("\nFixed effects:\n")
    coef <- data.frame(Value = x$B, Std.Error = x$B.se, Zscore = x$B.zscore, Pvalue = x$B.pvalue)
    printCoefmat(coef, P.values = TRUE, has.Pvalue = TRUE)
    cat("\n")
}

print.communityPGLMM <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    summary.communityPGLMM(x, digits = digits)
}

######################################################
######################################################
# plot.communityPGLMM
######################################################
######################################################
plot.communityPGLMM <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    if (!require(plotrix)) {
        stop("The 'plotrix' package is required to plot images from this function")
    }
    
    W <- data.frame(Y = x$Y, sp = x$sp, site = x$site)
    Y <- reshape(W, v.names = "Y", idvar = "sp", timevar = "site", direction = "wide")
    Y <- Y[, 2:dim(Y)[2]]
    
    par(mfrow = c(1, 1), las = 1, mar = c(4, 4, 2, 2) - 0.1)
    
    color2D.matplot(Y, ylab = "species", xlab = "sites", main = "Observed values")
}


######################################################
######################################################
# communityPGLMM.predicted.values
######################################################
######################################################
communityPGLMM.predicted.values <- function(x, show.plot = TRUE, ...) {
    
    if (x$family == "gaussian") {
        V <- solve(x$iV)
        h <- matrix(0, nrow = length(x$Y), ncol = 1)
        for (i in 1:length(x$Y)) {
            h[i] <- as.numeric(V[i, -i] %*% solve(V[-i, -i]) %*% matrix(x$H[-i]))
        }
        predicted.values <- h
    }
    
    if (x$family == "binomial") {
        h <- x$H + x$X %*% x$B
        predicted.values <- as.numeric(h)
    }
    
    if (show.plot == TRUE) {
        if (!require(plotrix)) {
            stop("The 'plotrix' package is required to plot images from this function")
        }
        
        W <- data.frame(Y = predicted.values, sp = x$sp, site = x$site)
        Y <- reshape(W, v.names = "Y", idvar = "sp", timevar = "site", direction = "wide")
        Y <- Y[, 2:dim(Y)[2]]
        par(mfrow = c(1, 1), las = 1, mar = c(4, 4, 2, 2) - 0.1)
        
        color2D.matplot(Y, ylab = "species", xlab = "sites", main = "Predicted values")
    }
    return(predicted.values)
}
