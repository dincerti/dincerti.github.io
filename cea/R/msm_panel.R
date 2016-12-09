# Fitting a muti-state model with unkown transition times
library("ggplot2")
library("msm")
library("cea")
library("MASS")
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
               covariates = list("1-2" = ~ sex), center = T)

## @knitr transition_intensity_matrix
# estimates
qmat <- qmatrix.msm(cav.msm, covariates=list(sex=1))
cav.msm$Qmatrices$baseline[2, 1] * exp(cav.msm$Qmatrices$sex[2, 1] * 
                                         (1 - cav.msm$qcmodel$covmeans))
exp(cav.msm$Qmatrices$logbaseline[2, 1] + cav.msm$Qmatrices$sex[2, 1] * 
                                         (1 - cav.msm$qcmodel$covmeans))
cav.q[2, 1]

# standard errors
sqrt(diag(cav.msm$covmat))

## @knitr transition_probability
# 1-year transition probability
pmat <- pmatrix.msm(cav.msm, t=1, covariates = list(sex = 1))
MatrixExp(qmat$estimates)
print(pmat)

# 10-year transition probability 
pmat10 <- pmatrix.msm(cav.msm, t = 10, covariates = list(sex = 1))
x <- c(1, 0, 0, 0)
for (t in 1:10){
  x <- x %*% phat
}
print(x)
print(pmat10[1, ])

## @knitr posterior_distribution
#boot.msm(cav.msm, stat = pmatrix.msm, B  = 10, file = NULL, cores = NULL)
n <- 1000
x <- cav.msm
base.ind <- which(names(x$estimates) == "qbase")
psamp <- mvrnorm(n, cav.msm$estimates, cav.msm$covmat)
pmatl <- vector(mode = "list", n)
for (i in 1:n){
  q <- matrix(0, nrow = nrow(x$qmodel$imatrix), ncol = ncol(x$qmodel$imatrix))
  q[t(x$qmodel$imatrix) == 1] <- exp(psamp[i, base.ind])
  q <- t(q)
  diag(q) <- -rowSums(q)
  pmatl[[i]] <- MatrixExp(q)
}


