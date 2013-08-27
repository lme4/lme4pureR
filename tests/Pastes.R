library(lme4)         # don't attach because of name conflicts
library(lme4pureR)
library(minqa)

form <- strength ~ (1|sample) + (1|batch)
devf <- plsform(form, Pastes, REML=FALSE)
bobyqa(c(1,1), devf, lower=c(0,0))[c("par","fval")]
mML <- lme4::lmer(form, Pastes,REML=FALSE)
getME(mML, "theta")
deviance(mML)
