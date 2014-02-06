library(lme4pureR)
library(minqa)
lmerReproduce <- FALSE
tol <- 1e-3

# sleepstudy
form <- Reaction ~ Days + (Days|Subject)
data(sleepstudy, package="lme4")
ll <- plsform(form, data = sleepstudy, REML=FALSE)
devf <- do.call(pls, ll)
opt <- minqa::bobyqa(ll$theta, devf, ll$lower, ll$upper)[c("par","fval")]
if(lmerReproduce){
  mML <- lme4::lmer(form, sleepstudy, REML=FALSE)
  par <- unname(lme4::getME(mML, "theta"))
  fval <- deviance(mML)
} else{
  par <- c(0.92922239, 0.01816504, 0.22264528)
  fval <- 1751.939
}
all.equal(par, opt$par, tolerance = tol)
all.equal(fval, opt$fval, tolerance = tol)


# Pastes
form <- strength ~ (1|sample) + (1|batch)
data(Pastes, package="lme4")
ll <- plsform(form, data = Pastes, REML=FALSE)
devf <- do.call(pls, ll)
opt <- minqa::bobyqa(ll$theta, devf, ll$lower, ll$upper)[c("par","fval")]
if(lmerReproduce){
  mML <- lme4::lmer(form, Pastes, REML=FALSE)
  par <- unname(lme4::getME(mML, "theta"))
  fval <- deviance(mML)
} else{
  par <- c(3.526904, 1.329914)
  fval <- 247.9945
}
all.equal(par, opt$par, tolerance = tol)
all.equal(fval, opt$fval, tolerance = tol)

# Dyestuff
form <- Yield ~ 1|Batch
data(Dyestuff,package="lme4")
ll <- plsform(form, data = Dyestuff, REML=FALSE)
devf <- do.call(pls, ll)
opt <- minqa::bobyqa(ll$theta, devf, ll$lower, ll$upper)[c("par","fval")]
if(lmerReproduce){
  mML <- lme4::lmer(form, Dyestuff, REML=FALSE)
  par <- unname(lme4::getME(mML, "theta"))
  fval <- deviance(mML)
} else{
  par <- 0.7525807
  fval <- 327.3271
}
all.equal(par, opt$par, tolerance = tol)
all.equal(fval, opt$fval, tolerance = tol)
