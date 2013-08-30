library(lme4pureR)
glmerReproduce <- FALSE
tol <- 1e-3

form <- cbind(incidence, size - incidence) ~ period + (1 | herd)
data(cbpp, package = 'lme4')
ll <- plsform(form, data = cbpp, family = binomial)
devf <- do.call(pirls, c(ll, list(family=binomial)))
opt <- minqa:::bobyqa(c(ll$theta, 0,0,0,0), devf)
if(glmerReproduce){
    mML <- lme4::glmer(form, data = cbpp, family = binomial)
    par <- unname(c(lme4::getME(mML, "theta"),getME(mML, "beta")))
    fval <- deviance(mML)
} else{
    par <- c(0.6420622, -1.3983429, -0.9919250, -1.1282162, -1.5797454)
    fval <- 184.0531
}
all.equal(par, opt$par, tolerance = tol)
all.equal(fval, opt$fval, tolerance = tol)
