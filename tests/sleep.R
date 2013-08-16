library(lme4)
library(minqa)
lmod <- lFormula(Reaction ~ Days + (Days|Subject), sleepstudy)
bobyqa(c(1, 0, 1),
       function(theta) pls(lmod,theta,sleepstudy$Reaction,REML=FALSE),
       lower=c(0,-Inf,0))$par
mML <- lmer(Reaction ~ Days + (Days|Subject),
            sleepstudy, REML = FALSE)
getME(mML, "theta")
