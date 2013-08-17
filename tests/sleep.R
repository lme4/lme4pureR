library(lme4)
library(lme4pureR)
library(minqa)
lmod <- lFormula(Reaction ~ Days + (Days|Subject), sleepstudy)
devfun <- function(theta) pls(lmod,theta,REML=FALSE)
bobyqa(c(1, 0, 1),
       function(theta) pls(lmod,theta,REML=FALSE),
       lower=c(0,-Inf,0))$par
mML <- lmer(Reaction ~ Days + (Days|Subject),
            sleepstudy, REML = FALSE)
getME(mML, "theta")


lme4pureR::lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
