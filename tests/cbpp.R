library(lme4)
library(minqa)
glmod <- glFormula(incidence / size ~ period + (1 | herd), weights = size,
             family = binomial, data = cbpp)
# todo: add actual pirls tests here
