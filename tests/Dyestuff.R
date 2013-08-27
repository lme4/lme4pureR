library(lme4)   # avoid name conflict with lmer
library(lme4pureR)
library(minqa)

form <- Yield ~ 1|Batch
data(Dyestuff,package="lme4")
devf <- plsform(form, Dyestuff, REML=FALSE)
bobyqa(1, devf, lower=0, upper=Inf)[c("par","fval")]
mML <- lme4::lmer(form, Dyestuff, REML=FALSE)
getME(mML, "theta")
deviance(mML)
