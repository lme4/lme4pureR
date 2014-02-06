library(lme4pureR)
library(minqa)
glmerReproduce <- FALSE
tol <- 1e-3

form <- cbind(incidence, size - incidence) ~ period + (1 | herd)
data(cbpp, package = 'lme4')
ll <- plsform(form, data = cbpp, family = binomial)
devf <- do.call(pirls, c(ll, list(family=binomial)))
rho <- environment(devf)
opt <- minqa:::bobyqa(c(ll$theta, rho$beta), devf)
if(glmerReproduce){
    mML <- lme4::glmer(form, data = cbpp, family = binomial)
    par <- unname(c(lme4::getME(mML, "theta"), getME(mML, "beta")))
    fval <- deviance(mML)
} else{
    par <- c(0.6420622, -1.3983429, -0.9919250, -1.1282162, -1.5797454)
    fval <- 184.0531
}
all.equal(par, opt$par, tolerance = tol)
all.equal(fval, opt$fval, tolerance = tol)

options(show.signif.stars = FALSE)
form <- use ~ age + I(age^2) + ch + urban + (1|district)
data(Contraception, package = 'mlmRev')
Contraception <- within(Contraception, ch <- factor(as.numeric(as.integer(livch)>1L)))
glm0 <- glm(lme4:::nobars(form),binomial,Contraception)
ll <- plsform(form, data = Contraception, family = binomial)
devf <- do.call(pirls, c(ll, list(family=binomial,eta=glm0$linear.predictor, tol=1e-6)))
rho <- environment(devf)

if(FALSE){ # FIXME: why step-halving problem?
opt <- minqa:::bobyqa(c(ll$theta, rho$beta), devf)
opt <- minqa:::bobyqa(c(ll$theta, coef(glm0)), devf)
opt <- minqa:::bobyqa(opt$par, devf)
}

form <- use ~ age + I(age^2) + ch + urban + (1|district)
data(Contraception, package = 'mlmRev')
Contraception <- within(Contraception, ch <- factor(as.numeric(as.integer(livch)>1L)))
glm0 <- glm(lme4:::nobars(form),binomial,Contraception)
ll <- plsform(form, data = Contraception, family = binomial)
devf <- do.call(pirls, c(ll, list(family=binomial,eta=glm0$linear.predictor, tol=1e-6)))

#body(devf)[8] <- parse("olducden <- updatemu(u); print(ucden); print(olducden)")
#body(devf)[7:9] <- parse(text = "olducden <- updatemu(u); print(ucden); print(olducden)")
paropt <- c(0.474010082, -1.006445615,  0.006255540, -0.004635385,  0.860439478, 0.692959336)
devf(paropt)
if(FALSE) opt <- minqa:::bobyqa(paropt, devf) # FIXME: step-halving again
devf(paropt)

library(lme4)
glmer0 <- glmer(form, data = Contraception, family = binomial, nAGQ = 0)
ll <- plsform(form, data = Contraception, family = binomial)
devf <- do.call(pirls, c(ll, list(family=binomial,eta=qlogis(getME(glmer0, 'mu')), verbose=2L)))
if(FALSE) opt <- minqa:::bobyqa(c(glmer0@theta,glmer0@beta), devf) #FIXME



glmer0 <- glmer(form, data = Contraception, family = binomial, nAGQ = 0)
ll <- plsform(form, data = Contraception, family = binomial)
devf <- do.call(pirls, c(ll, list(family=binomial,eta=qlogis(getME(glmer0, 'mu')), tol=1e-6)))
if(FALSE) opt <- minqa:::bobyqa(c(glmer0@theta,glmer0@beta), devf) # FIXME

gmod <- glFormula(form, data = Contraception, family = binomial, nAGQ = 1)
devf <- do.call(mkGlmerDevfun, gmod)
devf <- updateGlmerDevfun(devf, gmod$reTrms)
optimizeGlmer(devf, stage=2)$par

form1 <- use ~ age + I(age^2) + ch + urban + (1|district) + (1|district:urban)

## > print(summary(gm2 <- glmer(form1,Contraception,binomial)), corr=FALSE)
## Generalized linear mixed model fit by maximum likelihood ['summary.merMod']
##  Family: binomial ( logit )
## Formula: use ~ age + I(age^2) + ch + urban + (1 | district) + (1 | district:urban) 
##    Data: Contraception 

##      AIC      BIC   logLik deviance 
##  2375.88  2414.85 -1180.94  2361.88 

## Random effects:
##  Groups         Name        Variance Std.Dev.
##  district:urban (Intercept) 0.31863  0.5645  
##  district       (Intercept) 0.00715  0.0845  
## Number of obs: 1934, groups: district:urban, 102; district, 60

## Fixed effects:
##              Estimate Std. Error z value Pr(>|z|)
## (Intercept) -1.032745   0.175593   -5.88  4.1e-09
## age          0.005876   0.007935    0.74     0.46
## I(age^2)    -0.004537   0.000724   -6.26  3.8e-10
## ch1          0.872691   0.149028    5.86  4.7e-09
## urbanY       0.764670   0.169725    4.51  6.6e-06

beta2 <- c(-1.03274535458238, 0.00587602227525672, -0.0045372320199417,
           0.872690858261632, 0.764670407439094)
theta2 <- c(0.564472021994911, 0.0845334785825269)

## dput(head(unname(getME(gm2,"u"))))
## c(-1.60622219325897, -0.981855964992823, -0.0339024704692014, 
##    0.505250674118111, -0.428344594466731, 1.10122436874433)

## dput(head(unname(getME(gm2,"mu"))))
## c(0.195157958203462, 0.26347290296241, 0.504168959587931, 0.436350391068322, 
## 0.145679693412778, 0.178092903689705)

form2 <- use ~ age + I(age^2) + ch + urban + (1 | district:urban)

## print(summary(gm3 <- glmer(form2,Contraception,binomial)), corr=FALSE)
## Generalized linear mixed model fit by maximum likelihood ['summary.merMod']
##  Family: binomial ( logit )
## Formula: use ~ age + I(age^2) + ch + urban + (1 | district:urban) 
##    Data: Contraception 

##      AIC      BIC   logLik deviance 
##  2373.88  2407.29 -1180.94  2361.88 

## Random effects:
##  Groups         Name        Variance Std.Dev.
##  district:urban (Intercept) 0.327    0.572   
## Number of obs: 1934, groups: district:urban, 102

## Fixed effects:
##              Estimate Std. Error z value Pr(>|z|)
## (Intercept) -1.032829   0.175651   -5.88  4.1e-09
## age          0.005862   0.007935    0.74     0.46
## I(age^2)    -0.004535   0.000724   -6.26  3.9e-10
## ch1          0.872716   0.149034    5.86  4.7e-09
## urbanY       0.766667   0.170601    4.49  7.0e-06

beta3 <- c(-1.03282920682185, 0.00586163322902896, -0.00453480196155543, 
           0.872715839072072, 0.766667032721774)
theta3 <- 0.571568547571509

## dput(head(unname(getME(gm3,"mu"))))
## c(0.195804532044805, 0.264187372529822, 0.505051933870551, 0.43723605274129, 
## 0.146198921135808, 0.178681337173397)

## dput(head(unname(getME(gm3,"u"))))
## c(-1.63882047908625, -1.02416872093547, -0.0343530888451908, 
## 0.510937117406129, -0.419566092717808, 1.10846888333543)
