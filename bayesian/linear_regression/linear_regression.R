## run a linear regression with informative priors by hand in R and C++
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
n <- 10000
beta <- c(1, 2, 1.5)
sigma <- 1.2
x <- cbind(rep(1, n), rnorm(n, 10, 5), rnorm(n, 5, 2))
colnames(x) <- paste0("x", seq(1, ncol(x)))
y <- x %*% beta + rnorm(n, 0, sigma)

## ---- PRIORS -----------------------------------------------------------------
b0 <- rep(0, length(beta))
B0inv <- solve(diag(10, length(beta)))
c0 <- .001
d0 <- .001
nsims <- 10000

## ---- MCMC -------------------------------------------------------------------
# MCMC pack
lm.dat <- data.frame(y = y, x)
ptm <- proc.time()
blm.mcmcpack <- MCMCregress(mcmc = nsims, y ~ x2 + x3, b0 = b0, B0 = B0inv,
                            c0 = c0, d0 = d0, data = lm.dat)
proc.time() - ptm

# gibbs sampler with R
ptm <- proc.time()
blm.r <- blmR(sims = nsims, x = x, y = y, b0 = b0, B0inv = B0inv, c0 = .001,
              d0 = .001, beta.start = beta, sigma2.start = sigma^2)
proc.time() - ptm

# gibbs sampler with C++
ptm <- proc.time()
blm.c <- blmC(sims = nsims, y = y, X = x, b0 = b0, B0inv = B0inv,
              sigma2 = sigma^2)
proc.time() - ptm

## ---- COMPARE PARAMETRS ------------------------------------------------------
summary(blm.mcmcpack)[[1]]
rbind(summary(mcmc(blm.r$beta))[[1]], summary(mcmc(blm.r$sigma2))[[1]])
rbind(summary(mcmc(blm.c$beta))[[1]], summary(mcmc(blm.c$sigma2))[[1]])

