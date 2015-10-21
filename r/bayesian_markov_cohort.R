
# This is the R code for the R Markdown file bayesian_markov_cohort.Rmd
setwd("C:/Users/Devin/Dropbox/Projects/dincerti.github.io")
rm(list=ls()) 

## ---- PARAMETERS -------------------------------------------------------------
## @knitr transition_table
tt <- matrix(c(1251, 350, 116, 17,
             0, 731, 512, 15,
             0, 0, 1312, 437),
             ncol = 4, nrow = 3, byrow = TRUE)

## @knitr relative_risk
rr.est <- .509
lrr.est <- log(.509)
lrr.upper <- log(.710)
lrr.se <- (lrr.upper - lrr.est)/qnorm(.975)

## ---- JAGS MODEL -------------------------------------------------------------
## @knitr jags_data
S <- 4
n <- apply(tt, 1, sum)
alpha <- rep(1, S)
mu.lrr <- lrr.est
tau.lrr <- 1/lrr.se^2
params <- c("rr", "p")
jags.data <- list("mu.lrr", "tau.lrr", "S", "n", "alpha", "tt") 

## @knitr run_jags
library("R2jags")
set.seed(100)
jagsfit <- jags(data = jags.data, parameters.to.save = params,
                 model.file = "jags/markov_cohort.txt", n.chains = 3,
                 n.iter = 10000, n.thin = 5, progress.bar = "none")


## @knitr jags_output
print(jagsfit)

## @knitr jags_combine_chains
jagsfit.mcmc <- do.call("rbind", as.mcmc(jagsfit))

## @knitr conjugate_prior
library(MCMCpack)
summary(rdirichlet(nrow(jagsfit.mcmc), tt[1, ] + alpha[1])[, 2])
summary(as.numeric(jagsfit.mcmc[, "p[1,2]"]))

## ---- MARKOV MODEL -----------------------------------------------------------
## @knitr costs_effects
c.zidovudine <- 2278
c.lamivudine <- 2086.50
c.0 <- c(2756 + c.zidovudine, 3052 + c.zidovudine, 9007 + c.zidovudine, 0)
c.1 <- c(c.0[1:3] + c.lamivudine, 0)
qolw <- c(1, 1, 1, 0)

## @knitr transition_matrix
TMatrix <- function(probs, rr){
  P0 <- matrix(c(probs, rep(0, 3), 1), nrow = 4, ncol = 4, byrow = TRUE)
  P1 <- P0
  P1[upper.tri(P1)] <- P1[upper.tri(P1)] * rr
  P1[1, 1] <- 1 - sum(P1[1, 2:4]) 
  P1[2, 2] <- 1 - sum(P1[2, 3:4])
  P1[3, 3] <- 1 - sum(P1[3, 4])
  return(list(P0 = P0, P1 = P1))
}

## @knitr simulation
source("r/markov.R")
ncycles <- 20
delta.c <- rep(NA, nrow(jagsfit.mcmc))
delta.e <- rep(NA, nrow(jagsfit.mcmc))
for (i in 1:nrow(jagsfit.mcmc)){
  tp <- TMatrix(probs = jagsfit.mcmc[i, 2:13], rr = jagsfit.mcmc[, "rr"][i])
  sim0 <- MarkovCohort(P = tp$P0,  z0 = c(1000, 0, 0, 0), ncycles = ncycles,
                       costs = c.0, qolw = qolw, 
                       discount = 0.06)  
  sim1 <- MarkovCohort(P = c(replicate(2, tp$P1, simplify = FALSE), 
                             replicate(ncycles - 2, tp$P0, simplify = FALSE)),
                       z0 = c(1000, 0, 0, 0), ncycles = ncycles,
                       costs = c(replicate(2, c.1, simplify = FALSE),
                                 replicate(ncycles - 2, c.0, simplify = FALSE)),
                       qolw = qolw, discount = 0.06)
  delta.c[i] <- (sum(sim1$c) - sum(sim0$c))
  delta.e[i] <- (sum(sim1$e) - sum(sim0$e))
}

## ---- DECISION ANALYSIS ------------------------------------------------------
## @knitr ce_plane
ce.dat <- data.frame(delta.c = delta.c, delta.e = delta.e)
ggplot(data = ce.dat, aes(x = delta.e, y = delta.c)) + geom_point() +
  xlab("Effectiveness difference") + ylab("Cost difference")

## @knitr ce_accept_curve
pce <- rep(NA, 10000)
for (k in 1:10000){
  pce[k] <- mean((k * delta.e - delta.c) > 0)
}
ceac.dat <- data.frame(k = seq(1, length(pce)), pce = pce )
ggplot(data = ceac.dat, aes(x = k, y = pce)) + geom_line() +
  xlab("Willigness to pay") + ylab("Probabiliy cost effective")


