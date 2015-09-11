
# This is the R code for the R Markdown file dynamic_twopart.Rmd
setwd("C:/Users/Devin/Dropbox/Projects/dincerti.github.io")

## ---- MORTALITY --------------------------------------------------------------
## @knitr lifetable
library(XML)
lt <- readHTMLTable("http://www.ssa.gov/oact/STATS/table4c6.html", 
                    stringsAsFactors = FALSE)
lt <- lt[[2]][, c("V1", "V2")]
lt <- lt[-c(1:2), ]
lt$V1 <- as.numeric(lt$V1)
lt$V2 <- as.numeric(lt$V2)
colnames(lt) <- c("age", "mrate")

## @knitr mort_fit
CSage <- function(x) (x - 65)/10
lt$c_age <- CSage(lt$age)
lt$qnorm_mrate <- qnorm(lt$mrate)
mrate.lm <- lm(qnorm_mrate ~ c_age +  I(c_age^2) + I(c_age^3), data = lt)

## @knitr dynamic_twopart_mortplot
library(ggplot2)
theme_set(theme_bw())
lt$phat <- predict(mrate.lm, lt)
ggplot(lt, aes(x = age, y = qnorm_mrate)) + geom_point(size = 1) + 
  geom_line(aes(y = phat), col = "blue") +
  xlab("Age") + ylab(expression(Phi^-1*"(mortality rate)"))

## ---- SIMULATING DATA --------------------------------------------------------
## @knitr sim_setup
library(mvtnorm) # draw from multivariate normal
library(data.table)
set.seed(100)

## @knitr initdat_func
InitData <- function(n){
  d.0 <- rbinom(n, 1, .80)
  y.0 <- d.0 * rlnorm(n, meanlog = 7.5, sdlog = 1)
  usage <- data.frame(min = seq(0, 85, 5), max = c(seq(4, 85, 5), 120), 
                      perc = c(6.5, 6.6, 6.7, 7.1, 7, 6.8, 6.5, 6.5, 6.8, 7.4,
                               7.2, 6.4, 5.4, 4, 3, 2.4, 1.9, 1.8)) 
  agegrp <- sample(nrow(usage), n, replace=TRUE, prob = usage$perc/100)
  age.0 <- round(usage$min[agegrp] + runif(n) * (usage$max[agegrp] - usage$min[agegrp])) 
  return(data.frame(age = age.0, y = y.0))
}

## @knitr true_params
alpha <- c(.5, .05, .5)
beta <- c(6, .1, .25)
kappa <- coef(mrate.lm)
sigma2 <- 1
Sigma <- matrix(c(.5, .25, .25, .3), nrow = 2, ncol = 2)

## @knitr sim_func
# FUNCTION TO SIMULATE DATA
sim <- function(sim.T = 5, sim.n = 10000, SIGMA, Alpha = alpha, 
                Beta = beta, Kappa = kappa, Sigma2 = sigma2) {
  
  # RANDOM INTERCEPTS
  b <- rmvnorm(sim.n, mean = rep(0, 2), sigma = SIGMA) 
  
  # INITIALIZE LOOP
  init <- InitData(sim.n)
  N <- sim.n * sim.T
  age <- c(init$age, rep(NA, N))
  l_d <- l_ly <- rep(NA, N + sim.n)  
  m <- c(rep(0, sim.n), rep(NA, N))
  d <- c(1 * (init$y > 0), rep(NA, N))
  y <- c(init$y, rep(NA, N))
  l_obs <- 1:sim.n 
  
  # RECURSIVELY SIMULATE MODEL
  for (t in 1:sim.T){
    # update time varying data
    obs <- (l_obs[1] + sim.n):(l_obs[sim.n] + sim.n)
    age[obs] <- age[l_obs] + 1
    c_age <- CSage(age[obs])
    l_d[obs] <- d[l_obs]
    l_ly[obs] <- ifelse(y[l_obs] > 0, log(y[l_obs]), 0)
    x1 <- as.matrix(data.table(int = 1, c_age = c_age, l_d = l_d[obs]))
    x2 <- as.matrix(data.table(int = 1, c_age = c_age, l_ly = l_ly[obs]))
    
    # mortality 
    z <- as.matrix(data.table(int = 1, c_age = c_age, c_age2 = c_age^2, 
                              c_age3 = c_age^3))
    m.tmp <- rbinom(sim.n, 1, pnorm(z %*% Kappa))
    m[obs] <- ifelse(m[l_obs] == 0, m.tmp, 1) # set m = 1 in all periods after death
    
    # expenditures
    d[obs] <- rbinom(sim.n, 1, pnorm(x1 %*% Alpha + b[, 1]))
    y[obs] <- d[obs] * rlnorm(sim.n, meanlog = x2 %*% Beta + b[, 2], 
                              sdlog = sqrt(Sigma2))
    
    # counter
    l_obs <- obs 
  }
  
  # RESULTS
  dat <- data.table(id = rep(seq(1, sim.n), (sim.T + 1)), 
                    year = rep(seq(0, sim.T), each = sim.n),
                    y = y, d = d, m = m, b1 = rep(b[, 1], each = (sim.T + 1)), 
                    b2 = rep(b[, 2], each = (sim.T + 1)),
                    age = age, l_d = l_d, l_ly = l_ly)
  dat <- dat[order(id, year)]
  dat[, int := 1]
  dat[, c_age := CSage(age)]
  dat[, ly := ifelse(y > 0, log(y), NA)]
  dat <- dat[m == 0] # drop individuals after death
  return(dat)
}

## ---- A MODEL WITH NO RANDOM INTERCEPTS --------------------------------------
## @knitr sim_noint
library(mvtnorm) 
dat <- sim(SIGMA = Sigma * 0)
dat <- dat[year > 0]

## @knitr d_probit
d.probit <- glm(d ~ c_age + l_d, family=binomial(link=probit), dat)
cbind(coef(d.probit), confint.default(d.probit), alpha)

## @knitr ly_lm
ly.lm <- lm(ly ~ c_age + l_ly, dat)
cbind(coef(ly.lm), confint(ly.lm), beta)

## @knitr cov95
nsims <- 1000
cov.95 <- matrix(NA, length(beta), nsims)
for (i in 1:nsims){
  dat <- sim(sim.n = 1000, SIGMA = Sigma * 0) 
  lm.sim <- lm(ly ~ c_age + l_ly, dat[year > 0])
  beta.hat <- coef(lm.sim)
  beta.se <- summary(lm.sim)$coef[, 2]
  cov.95[, i] <- abs(beta - beta.hat) < abs(qnorm(.025)) * beta.se
}
apply(cov.95, 1, mean)

## ---- MCMC -------------------------------------------------------------------
## @knitr sim
dat <- sim(sim.T = 5, SIGMA = Sigma)
dat <- dat[year > 0]
ni <- dat[, .N, by = "id"]$N
x1 <- as.matrix(dat[, .(int, c_age, l_d)])
x2 <- as.matrix(dat[, .(int, c_age, l_ly)])
b1 <- dat[, .(b1 = mean(b1)), by = "id"] 
b2 <- dat[, .(b2 = mean(b2)), by = "id"] 

## @knitr mcmc_init
library(MCMCpack)   # iwish distribution used in updating Sigma
priors <- list(malpha = rep(0, ncol(x1)), Valpha = diag(10, ncol(x1)), 
               mbeta = rep(0, ncol(x2)), Vbeta = diag(10, ncol(x2)),
               sigma2.shape = 1, sigma2.rate = 1,
               Sigma.v = 3, Sigma.S = diag(2))
inits <- list()
inits$alpha <- rnorm(length(alpha), coef(d.probit), summary(d.probit)$coef[, 2])
inits$beta <- rnorm(length(beta), coef(ly.lm), summary(ly.lm)$coef[, 2])
inits$sigma2 <- rnorm(1, summary(ly.lm)$sigma^2, .2)
inits$Sigma <- riwish(v = 75, S = 79 * Sigma)
inits$b <- cbind(b1$b1, b2$b2)

## @knitr gibbs
source("r/dynamic_twopart_mcmc.R")
gibbs <- Gibbs(nsim = 10000, thin = 10, burn = 5000, y = dat$y, 
               x1 = x1, x2 = x2, ni = ni, id = dat$id, 
               priors = priors, init = inits)
save(gibbs, file = "output/gibbs.RData")

## @knitr convert_mcmc
load("output/gibbs.RData")
library(coda)
alpha.mcmc <- as.mcmc(gibbs$alpha)
beta.mcmc <- as.mcmc(gibbs$beta)
sigma2.mcmc <- as.mcmc(gibbs$sigma2)
Sigma.mcmc <- as.mcmc(gibbs$Sigma) 

## @knitr dynamic_twopart_traceplot
plot(sigma2.mcmc, density = FALSE, 
     main = expression("Traceplot of"~sigma[epsilon]^2))

## @knitr post_quantiles
cbind(summary(alpha.mcmc)$quant, alpha)
cbind(summary(beta.mcmc)$quant, beta)
c(summary(sigma2.mcmc)$quant, sigma2)
cbind(summary(Sigma.mcmc)$quant, c(Sigma))