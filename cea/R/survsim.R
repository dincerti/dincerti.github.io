library("ggplot2")
library("msm")
library("cea")
theme_set(theme_bw())

## ---- SAMPLE FROM CUMULATIVE HAZARD ------------------------------------------
## @knitr haz_sample_fun
cumhaz_sample <- function(time, Haz, size = 1, replace = TRUE){
  ### Function to sample from a distribution given by a cumulative hazard
  ### Input:
  ###     time: vector containing time points to sample from
  ###     Haz: the cumulative hazard at those time points
  ###     size: the size of the sample
  ###     replace: sample with or without replacement? (default=TRUE, with replacement)
  ### Output:
  ###     A random realisation from H
  p <- diff(c(0, 1 - exp(-Haz)))
  p <- c(p, exp(-Haz[length(Haz)])) # add probability of sampling time=Inf 
  return(sample(c(time, Inf), size = size, prob = p, replace = replace))
}

## @knitr weibull_hazard
cumhaz_weibull <- function(time, shape, scale){
  haz <- shape/scale * (time/scale)^(shape - 1)
  cumhaz <- cumsum(haz[-1] * diff(time))
  cumhaz <- c(haz[1], cumhaz)
  return(cumhaz)
}
cumhaz.weibull <- data.frame(t = seq(0, 50, .01))
cumhaz.weibull$Haz <- with(cumhaz.weibull, 
                           cumhaz_weibull(t, shape = 1.5, scale = 10))
ggplot(cumhaz.weibull, aes(x = t, y = Haz)) + geom_line() + xlab("Time") + 
  ylab("Cumulative hazard")

## @knitr haz_sample_weibull
n <- 10000
weicumhaz.sim <- cumhaz_sample(cumhaz.weibull$t, 
                            cumhaz.weibull$Haz, n)
weidist.sim <- rweibull(n, shape = 1.5, scale = 10)
sim <- data.frame(t = c(weicumhaz.sim, weidist.sim),
                  lab = rep(c("Sample from cumulative hazard",
                              "Sample from Weibull distribution"), each = n))
ggplot(sim, aes(x = t, col = lab)) + geom_density() + xlab("Simulated time") + 
  ylab("Density") + theme(legend.position="bottom",
                                          legend.title = element_blank()) 

## ---- SAMPLE FROM LIFETABLE --------------------------------------------------
## @knitr rpexp
lifetable_male$rate <- -log(1 - lifetable_male$qx) # assume constant rate
pexp.sim <- rpexp(n, rate = lifetable_male$rate,
                  t = lifetable_male$age)
summary(pexp.sim)

## @knitr sample_survdist
lifetable_male$surv <- lifetable_male$lx/lifetable_male$lx[1]
survdist.sim <- sample(lifetable_male$age, 
                        size = n, 
                        prob = diff(c(0, 1 - lifetable_male$surv)), 
                        replace = T)
summary(survdist.sim)

## @knitr sample_binomial
qx_sample <- function(ages, qx, nsims){
  age.death <- rep(NA, nsims)
  n.ages <- length(ages)
  for (i in 1:nsims){
    death <- rbinom(n.ages, 1, prob)
    age.death[i] <- ages[which(death == 1)[1]]
  }
  return(age.death)
}
qx.sim <- qx_sample(lifetable_male$age, lifetable_male$qx, n)
summary(qx.sim)

## @knitr lifetable_sim_compare
lifetable.sim <- data.frame(sim = c(pexp.sim, survdist.sim, qx.sim),
                            lab = rep(c("Sample from piecewise exponential", 
                                    "Sample from survival curve",
                                    "Sample from binomial distributions"), 
                                    each = n))
ggplot(lifetable.sim, aes(x = sim, col = lab)) + geom_density() + 
  xlab("Simulated age at death") + 
  ylab("Density") + theme(legend.position="bottom",
                          legend.title = element_blank()) +
  guides(col = guide_legend(nrow=2, byrow=TRUE))


