library(lme4pureR)
form <- Reaction ~ Days + (Days|Subject)
data(sleepstudy, package="lme4")
ll <- plsform(form, sleepstudy, REML=FALSE)
devf <- do.call(pls, ll)
dput(minqa::bobyqa(ll$theta, devf, ll$lower, ll$upper)[c("par","fval")])

runLmer <- FALSE
if(runLmer){
  mML <- lme4::lmer(form, sleepstudy, REML=FALSE)
  dput(lme4::getME(mML, "theta"))
  dput(deviance(mML))
}
