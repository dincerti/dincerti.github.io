# This is the R code for the R Markdown file multinomial.Rmd
library("nnet")
library("MCMCpack")
set.seed(101)

## ---- SIMULATE DATA ----------------------------------------------------------
## @knitr sim_setup
n <- 100000
x <- data.frame(int = rep(1, n),
                age = round(runif(n, 45, 70)),
                male = rbinom(n, 1, .5))
beta.good <- matrix(c(-.05, -.02, -.02, -.05, -.03, -.02), nrow = 3, ncol = 2)
beta.fair <- matrix(c(-.12, .05, .02, -.0886, .01, .02), nrow = 3, ncol = 2)
beta.poor <- matrix(c(0.02, .01, .05, 0.1, .05, .02), nrow = 3, ncol = 2)
colnames(beta.poor) <- c("Fair", "Poor")
rownames(beta.poor) <- colnames(x)
beta.poor

## @knitr mlogit_predict_prob
mlogit_prob <- function(x, beta) {
  log.odds <- as.matrix(x) %*% beta
  odds <- cbind(exp(0), exp(log.odds))
  odds.sum <- apply(odds, 1, sum) 
  p <- odds/odds.sum
  return(p)
}
p.good <- mlogit_prob(x, beta.good)
p.fair <- mlogit_prob(x, beta.fair)
p.poor <- mlogit_prob(x, beta.poor)

## @knitr sim_mlogit_data
sim_mlogit <- function(x, beta){
  p <- mlogit_prob(x, beta)
  m = t(apply(p, 1, rmultinom, n = 1, size = 1))
  dat <- cbind(x, data.frame(trans = apply(m, 1, function(x) which(x==1))))
  return(dat)
}
dat.good <- sim_mlogit(x, beta.good)
dat.fair <- sim_mlogit(x, beta.fair)
dat.poor <- sim_mlogit(x, beta.poor)

## ---- MAXIMUM LIKELIHOOD -----------------------------------------------------
## @knitr ml_mnl
ml.fair <- multinom(trans ~ age + male, 
                   data = dat.fair)

## @knitr predic_prob
newdat <- matrix(c(1, 1, 55, 55, 0, 1), nrow = 2)
colnames(newdat) <- colnames(model.matrix(ml.fair))
mlogit_prob(as.matrix(newdat), t(coef(ml.fair)))
predict(ml.fair, newdat, "probs")
mlogit_prob(newdat, beta)

## ---- BAYESIAN ---------------------------------------------------------------
## @knitr bayesian_mnl
bayes.good <- MCMCmnl(trans ~ age + male, 
                      data = dat.good,
                      mcmc.method = "IndMH")
bayes.fair <- MCMCmnl(trans ~ age + male, 
                      data = dat.fair,
                      mcmc.method = "IndMH")
bayes.poor <- MCMCmnl(trans ~ age + male, 
                      data = dat.poor,
                      mcmc.method = "IndMH")

## @knitr compare_coef
rbind(c(coef(ml.fair)),
      summary(bayes.fair)$stat[, "Mean"])
rbind(c(summary(ml.fair)$standard.errors),
      summary(bayes.fair)$stat[, "SD"])

## ---- SAVE RESULTS -----------------------------------------------------------
## @knitr save
bayes <- vector(3, mode = "list")
bayes[[1]] <- as.matrix(bayes.good)
bayes[[2]] <- as.matrix(bayes.fair)
bayes[[3]] <- as.matrix(bayes.poor)
names(bayes) <- c("good", "fair", "poor")
saveRDS(as.matrix(bayes), "cea/output/mlogit_bayes.rds")
