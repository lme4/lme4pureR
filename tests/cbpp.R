library(lme4)
library(minqa)
library(lme4pureR)

# model structure
form <- cbind(incidence, size - incidence) ~ period + (1 | herd)
glmod <- glFormula(form, cbpp, binomial)
do.call(mkGlmerDevfun, glmod)
data(cbpp, package = 'lme4')
ll <- plsform(form, data = cbpp, family = binomial)
devf <- do.call(pirls, c(ll, list(family=binomial)))
devf(c(1,0,0,0,0))
bobyqa(c(1,0,0,0,0), devf)

# get initial values
gm1 <- glm(nobars(form), binomial, cbpp)
weights <- weights(gm1)
beta0 <- coef(gm1)
ratios <- gm1$y
eta <- gm1$linear.predictors

# create deviance function with `lme4pureR`
devf <- pirls(glmod, ratios, binomial(), weights=weights, eta=eta, nAGQ=1)

# run `bobyqa`
bobyqa(c(1,beta0), devf, lower=c(0,rep.int(-Inf,length(beta0))))$par
