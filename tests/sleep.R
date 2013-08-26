library(lme4)
#library(nloptwrap)
library(lme4pureR)
library(minqa)
form <- Reaction ~ Days + (Days|Subject)
if (FALSE) {            # vector-valued random effects not yet working
    devf <- plsform(form, sleepstudy, REML=FALSE)
    bobyqa(1, devf, lower=0, upper=Inf)[c("par","value")]
}
mML <- lmer(form, sleepstudy, REML=FALSE)
getME(mML, "theta")
deviance(mML)
