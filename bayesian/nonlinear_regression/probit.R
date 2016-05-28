## run a probit regression in R and C++
rm(list = ls())
library("Rcpp")
library("RcppArmadillo")
library("inline")
library("MCMCpack")
library("mvtnorm")
library("microbenchmark")
library("learnMCMC")

## ---- SIMULATE DATA ----------------------------------------------------------
set.seed(101)
n <- 100000
beta <- c(.2, .3, -.5)
x <- cbind(rep(1, n), rnorm(n, 10, 5), rnorm(n, 5, 2))
colnames(x) <- paste0("x", seq(1, ncol(x)))
y <- rbinom(n, 1, pnorm(x %*% beta))

## ---- PRIORS -----------------------------------------------------------------
b0 <- rep(0, length(beta))
B0inv <- solve(diag(10, length(beta)))
nsims <- 1000

## ---- MCMC -------------------------------------------------------------------
# MCMC pack
dat <- data.frame(y = y, x)
ptm <- proc.time()
bprobit.mcmcpack <- MCMCprobit(mcmc = nsims, y ~ x2 + x3, b0 = b0, B0 = B0inv,
                            data = dat)
proc.time() - ptm

# gibbs sampler with R
ptm <- proc.time()
bprobit.r <- bprobitR(sims = nsims, x = x, y = y, b0 = b0, B0inv = B0inv,
                      beta.start = beta)
proc.time() - ptm

# gibbs sampler with C++
ptm <- proc.time()
bprobit.c <- bprobitC(sims = nsims, X = x, y = y, b0 = b0, B0inv = B0inv,
                  beta_start = beta)
proc.time() - ptm

## ---- COMPARE PARAMETRS ------------------------------------------------------
summary(bprobit.mcmcpack)[[1]]
summary(mcmc(bprobit.r$beta))[[1]]
summary(mcmc(bprobit.c$beta))[[1]]

