# Fitting a muti-state model with unkown transition times
library("ggplot2")
library("msm")
library("cea")
theme_set(theme_bw())

## ---- FIT MULTI-STATE MODEL WITH PANEL DATA ----------------------------------
## @knitr observed_transitions
print(statetable.msm(state, PTNUM, data  = cav))
print(statetable.msm(state, ptnum, data  = aneur))

## @knitr fit model
cav.Q <- rbind(c(0, 1, 0, 1),
             c(1, 0, 1, 1),
             c(0, 1, 0, 1),
             c(0, 0, 0, 0))
cav.msm <- msm(state ~ years, subject = PTNUM, data = cav,
                qmatrix = cav.Q, deathexact = 4, gen.inits = TRUE,
               covariates = ~ sex)

## @knitr transition_probability
qmatrix.msm(cav.msm, covariates=list(sex=1))
phat <- pmatrix.msm(cav.msm, t=1)
x <- c(1, 0, 0, 0)
for (t in 1:10){
  x <- x %*% phat
}
print(x)
pmatrix.msm(cav.msm, t=10)[1, ]

## ---- SAMPLE FROM CUMULATIVE HAZARD ------------------------------------------