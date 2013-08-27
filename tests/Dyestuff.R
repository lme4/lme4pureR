library(lme4pureR)
form <- Yield ~ 1|Batch
data(Dyestuff,package="lme4")
ll <- plsform(form, Dyestuff, REML=FALSE)
devf <- do.call(pls, ll)
dput(minqa::bobyqa(1, devf, lower=0, upper=Inf)[c("par","fval")])
mML <- lme4::lmer(form, Dyestuff, REML=FALSE)
dput(lme4::getME(mML, "theta"))
dput(deviance(mML))

