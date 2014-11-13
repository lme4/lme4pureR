library(lme4pureR)

## only *one* of these two ?
library(nloptwrap)
library(minqa)

## library(lme4)
data(sleepstudy, package="lme4")
lmod <- lme4::lFormula(Reaction ~ Days + (Days|Subject), sleepstudy)

devf <- pls(lmod,sleepstudy$Reaction)
bobyqa(c(1, 0, 1), devf, lower=c(0,-Inf,0))[c("par","value")]
mML <- lme4::lmer(Reaction ~ Days + (Days|Subject),
                  sleepstudy, REML = FALSE)
lme4::getME(mML, "theta")
lme4::deviance(mML)
