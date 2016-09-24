library("ggplot2")
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
  p <- diff(c(0, 1-exp(-Haz)))
  p <- c(p, exp(-Haz[length(Haz)])) # add probability of sampling time=Inf
  return(sample(c(time, Inf), size = size, prob = p, replace = replace))
}

## @knitr weibull_hazard
cumhaz_weibull <- function(time, shape, scale){
  #haz <- shape/scale * (time/shape)^(shape - 1)
  f <- dweibull(time, shape, scale)
  s <- 1 - pweibull(time, shape, scale)
  haz <- f/s
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