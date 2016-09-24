# This is the R code for the R Markdown file multinomial.Rmd
library("nnet")
library("MCMCpack")
set.seed(101)

## ---- SIMULATE DATA ----------------------------------------------------------
## @knitr sim_mlogit
# set up data
n <- 100000
x <- data.frame(int = rep(1, n),
                age = round(runif(n, 45, 70)),
                male = rbinom(n, 1, .5))
beta <- matrix(c(-.12, .05, .02, -.0886, .01, .02), nrow = 3, ncol = 2)

# simulate multinomial logit probability
mlogit_prob <- function(x, beta) {
  log.odds <- as.matrix(x) %*% beta
  odds <- cbind(exp(0), exp(log.odds))
  odds.sum <- apply(odds, 1, sum) 
  p <- odds/odds.sum
  return(p)
}
p <- mlogit_prob(x, beta)
sim_mlogit <- function(x, beta){
  p <- mlogit_prob(x, beta)
  m = t(apply(p, 1, rmultinom, n = 1, size = 1))
  dat <- cbind(x, data.frame(trans = apply(m, 1, function(x) which(x==1))))
  return(dat)
}
dat <- sim_mlogit(x, beta)

## ---- ESTIMATE MODEL ---------------------------------------------------------
## @knitr ml_multionomial_logit
mlogit.ml <- multinom(trans ~ age + male, 
                   data = dat)

## @knitr bayesian_multionomial_logit
mlogit.bayes <- MCMCmnl(trans ~ age + male, 
                        data = dat,
                        mcmc.method = "IndMH")

## @knitr compare_coef
rbind(c(coef(mlogit.ml)),
      summary(mlogit.bayes)$stat[, "Mean"])
rbind(c(summary(mlogit.ml)$standard.errors),
      summary(mlogit.bayes)$stat[, "SD"])

## @knitr predic_prob
newdat <- matrix(c(1, 1, 55, 55, 0, 1), nrow = 2)
mlogit_prob(as.matrix(newdat), t(coef(mlogit.ml)))
mlogit_prob(newdat, beta)

## ---- SAVE RESULTS -----------------------------------------------------------
saveRDS(matrix(mlogit.bayes), "cea/")
