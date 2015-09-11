---
layout: rmd
title: Dynamic Two-part Models
---
* TOC
{:toc}


### Overview
This page describes the dynamic two-part model that I have used in my research (see [here](papers\longterm_spending.pdf)). It extends the [cross-sectional two-part model](twopart.html) to a longitudinal setting and provides the R code necessary for estimation.  

With longitudinal data, it is necessary to model the persistence in spending from one period to the next. The model does this in two ways: first, it includes lagged dependent variables and second, it allows for individual specific random regression coefficients. I focus on the case in which only intercepts vary across individuals so the model can be referred to as a dynamic random-intercept two-part model.

### The Model
Let $d_{it}\equiv I(d_{it}^{\ast} > 0) = I(y_{it} > 0)$ so that the latent variable $d_{it}^{\ast}$ describes whether expenditures are positive or zero. The model can then be written as,

$$
\begin{aligned}
d_{it}^{\ast} &= x_{1it}^{T}\alpha + b_{1i} + \epsilon_{1it},\\
\log (y_{it}|d_{it}^{\ast} >0) &=  x_{2it}^{T}\beta + b_{2i} + \epsilon_{2it},
\end{aligned}
$$

where $\alpha$ and $\beta$ are the vectors of coefficients for the explanatory variables, and $b_{1i}$ and $b_{2i}$ are the random intercepts. We assume that $\epsilon_{1it}\sim N(0, 1)$ so that the binomial component is a probit model and $\epsilon_{2it}\sim N(0, \sigma^2_\epsilon)$ so that the continuous component is a lognormal model. Importantly, both $x_{1it}$ and $x_{2it}$ contain lagged values of (functions of) spending. In particular, $x_{1it}$ contains $d_{it-1}$ and $x_{2it}$ contains a variable equal to $\log(y_{it-1})$ if $y_{it-1}>0$ and $0$ otherwise. Other functions of lagged spending are possible but we will use those here for illustrative purposes.

The random intercepts are assumed to be jointly normal,

$$
\begin{aligned}
b_i = \begin{bmatrix}
b_{1i}\\
b_{2i}
\end{bmatrix} &\sim  N
\begin{bmatrix}
\begin{pmatrix}
0 \\
0
\end{pmatrix},
\Sigma_b = 
\begin{pmatrix}
\sigma_1^2 & \rho\sigma_{1}\sigma_{2}\\
\rho\sigma_{1}\sigma_{2} & \sigma_{2}^2 \\
\end{pmatrix}
\end{bmatrix},
\end{aligned}
$$

where $\rho$ is the correlation between $b_{1i}$ and $b_{2i}$. It is worth noting that if the design matrices, $x_{1it}$ and $x_{2it}$, contain an intercept, then $b_{1i}$ and $b_{2i}$ can be treated as error terms that are constant across individuals. In some disciplines, such as economics, these error terms are referred to as unobserved heterogeneity. 

To illustrate, consider the second part of the two-part model and let $\eta_{it} = \epsilon_{2it} + b_{2i}$. The variances and covariances of $\eta_{it}$ are then,

$$
\begin{aligned}
\rm{Var}(\eta_{it})&=\sigma^2_\epsilon + \sigma_2^2\\
\rm{Cov}(\eta_{it}, \eta_{is}) &= b_{2i} \; \rm{for}\; s\neq t\\
\rm{Cov}(\eta_{it}, \eta_{js}) &= 0\; \rm{for}\; s\neq t\; \rm{and}\; i \neq j.
\end{aligned}
$$

In other words the error terms are correlated within individuals but not over time across different individuals.

### Mortality
In order to model spending over time it is necessary to model mortality as well. Here, we ignore possible correlations between spending and mortality and assume that individuals die during during period $t$ at a rate, $m_{it}$, that is increasing with age. We model the death rate with a simple probit model,

$$
Pr(m_{it}=1)=\Phi(z_{it}^{T}\kappa),
$$

where $z_{it}$ contains an intercept and age covariates, and $\kappa$ is the corresponding vector of coefficients. To obtain a reasonable estimate for $\kappa$, we can download  estimates of the probability of dying within one year by age from an actuarial life table published by Social Security using the `XML` package. For simplicity, we will use the death rates for males.

The function `readHTMLTable` returns a list of HTML tables. We extract the first and second columns from the second table which contains the information we want.

{% highlight r %}
library(XML)
lt <- readHTMLTable("http://www.ssa.gov/oact/STATS/table4c6.html", 
                    stringsAsFactors = FALSE)
lt <- lt[[2]][, c("V1", "V2")]
lt <- lt[-c(1:2), ]
lt$V1 <- as.numeric(lt$V1)
lt$V2 <- as.numeric(lt$V2)
colnames(lt) <- c("age", "mrate")
{% endhighlight %}
With the data in hand, we can model the death rates as a function of age. To ensure that coefficient estimates are on a reasonable scale and that the intercept has an interesting interpretation, we center and scale all age variables using the function `CSage` below. Since we are assuming that death rates can be modeled with a probit model, it follows that $\Phi^{-1}(p_{it}) = z_{it}^{T}\kappa$ where $p_{it}= Pr(m_{it}=1)$. We can therefore estimate the regression coefficients with simple OLS.

{% highlight r %}
CSage <- function(x) (x - 65)/10
lt$c_age <- CSage(lt$age)
lt$qnorm_mrate <- qnorm(lt$mrate)
mrate.lm <- lm(qnorm_mrate ~ c_age +  I(c_age^2) + I(c_age^3), data = lt)
{% endhighlight %}
The model fits the data quite well after age 25 or so, which suggests that we can predict mortality very accurately as a function of age.

{% highlight r %}
library(ggplot2)
theme_set(theme_bw())
lt$phat <- predict(mrate.lm, lt)
ggplot(lt, aes(x = age, y = qnorm_mrate)) + geom_point(size = 1) + 
  geom_line(aes(y = phat), col = "blue") +
  xlab("Age") + ylab(expression(Phi^-1*"(mortality rate)"))
{% endhighlight %}

<img src="/figs/dynamic_twopart_mortplot-1.png" title="plot of chunk dynamic_twopart_mortplot" alt="plot of chunk dynamic_twopart_mortplot" style="display: block; margin: auto;" />

### Simulating Data
A great way to understand a model is to simulate data according to the assumed data generating process. That is, first set the parameters to their "true" values and simulate the model given these parameter values. Next, estimate parameters based on the simulated data using a chosen statistical method and compare the estimated values to the true values.

Before doing this, we will load a couple of R packages. The `mvtnorm` package can be used to draw the random intercepts from their bivariate normal distribution and the `data.table` package is useful for manipulating larger datasets. We should also set the seed to ensure that the results are reproducible.

{% highlight r %}
library(mvtnorm) # draw from multivariate normal
library(data.table)
set.seed(100)
{% endhighlight %}
In order to simulate the model, we will need period $0$ data, which is assumed to be known as baseline. For simplicity, we will use small design matrices, $x_{1}$ and $x_{2}$, that contain an intercept, a covariate for age, and a function of lagged spending. The model is therefore only dependent on initial spending levels and age. The function `InitData` creates the necessary initial values for a desired sample size. Values for $y$ and age are drawn from distributions so that they are resonably consistent with observed health spending and the age distribution in the United States.

{% highlight r %}
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
{% endhighlight %}
The true values for the parameters are set as follows.

{% highlight r %}
alpha <- c(.5, .05, .5)
beta <- c(6, .1, .25)
kappa <- coef(mrate.lm)
sigma2 <- 1
Sigma <- matrix(c(.5, .25, .25, .3), nrow = 2, ncol = 2)
{% endhighlight %}
We can now create a function that simulates longitudinal data using the model. The simulation begins with a set number of individuals alive during period $1$. Expenditures are predicted for a chosen number of simulation periods, say $T$, althogh some individuals die before reaching period $T$. 

{% highlight r %}
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
{% endhighlight %}
The function can be summarized as follows. First, it draws the random intercepts from their bivariate normal distribution. Second, it draws period $0$ data and initializes the vectors that we would like to store for analysis. Third, it recursively simulates the model for the desired number of periods by:

1. Updating the design matrices to reflect spending in the previous period and being 1 year older.
2. Drawing death indicators.
3. For survivors, simulating expenditures according to the two-part model. 
4. Incrementing the number of periods by 1 unless t = T.

The simulation ultimately returns a dataset with all relevant variables from the simulation.

### A Model with No Random Intercepts
Before moving on to the full model, we will investigate a model in which everyone has the same intercept. This is similar to a [cross-sectional two-part model](twopart.html) except that it includes lagged dependent variables. We will simulate data for 10, 000 individuals over 5 periods using the `sim` function.

{% highlight r %}
library(mvtnorm) 
dat <- sim(SIGMA = Sigma * 0)
dat <- dat[year > 0]
{% endhighlight %}
The parameters can be esimated using a probit model for the first part and OLS for the second part. The sample size is reasonably large so the estimated parameters should be close to their true values. 

{% highlight r %}
d.probit <- glm(d ~ c_age + l_d, family=binomial(link=probit), dat)
cbind(coef(d.probit), confint.default(d.probit), alpha)
{% endhighlight %}



{% highlight text %}
##                             2.5 %     97.5 % alpha
## (Intercept) 0.49885632 0.46938895 0.52832370  0.50
## c_age       0.04875333 0.04296868 0.05453798  0.05
## l_d         0.49120307 0.46285221 0.51955393  0.50
{% endhighlight %}
As expected, the estimated parameters are very close to the true values and that they fall with the 95\% confidence intervals. The OLS estimates are more precisely estimated and even closer to the true values.

{% highlight r %}
ly.lm <- lm(ly ~ c_age + l_ly, dat)
cbind(coef(ly.lm), confint(ly.lm), beta)
{% endhighlight %}



{% highlight text %}
##                             2.5 %    97.5 % beta
## (Intercept) 5.99886608 5.97283940 6.0248928 6.00
## c_age       0.09651274 0.09187113 0.1011544 0.10
## l_ly        0.24948909 0.24614595 0.2528322 0.25
{% endhighlight %}
it is important to note however that although the confidence intervals worked for this particular simulation, we can't be sure that the coverage probabilities are correct. To check coverage probabilities we would need to simulate the data multiple times, estimate the parameters for each simulation, and check whether the true values were contained within the 95\% confidence intervals for 95\% of the simulations. As an example, let's check coverage for the second part of the model by simulating data 1000 times,

{% highlight r %}
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
{% endhighlight %}



{% highlight text %}
## [1] 0.950 0.959 0.949
{% endhighlight %}
The confidence intervals are working as intended as the true coefficient values are contained within the 95\% confidence interval approximately 950 out of 1000 times.

### MCMC
Estimating a model with random intercepts that are correlated across both components of the two-part model is more difficult. As a result, this section is somewhat technical and requires Bayesian Markov Chain Monte Carlo (MCMC) methods. This means that we will need to write down the full joint posterior density for the model. Letting $\theta = (\alpha^T, \beta^T, \sigma^2_\epsilon)^T$ and $T_i$ represent the number of years that individual $i$ is observed before death, the conditional density of expenditures for indidivual $i$, $f(y_{i1}, y_{i2}, \ldots, y_{iT} |y_{i0}, \theta, b_i)$, can (under the model) be written as,

$$
\begin{aligned}
\prod_{t=1}^{T_i} f(y_{it}|y_{it-1}, \theta, b_i) &= \prod_{t=1}^{T_i}  (1-\Phi(\mu_{1it})^{1-d_{it}} \left[\Phi(\mu_{1it})\times \rm{LN}(y_{it};\mu_{2it}, \sigma^2_\epsilon)\right]^{d_{it}},\\
\mu_{1it}&= x^T_{1it} \alpha + b_{1i}, \\
\mu_{2it}&= x^T_{2it} \beta + b_{2i},
\end{aligned}
$$

where $\rm{LN}(\cdot)$ is the lognormal distribution. Given prior distributions, $p(\theta)$ and $p(\Sigma_b)$, the joint posterior density is then,

$$
\begin{aligned}
p(\theta, b_i, \Sigma_b|y) &\propto p(\theta)p(\Sigma_b)\prod_{i=1}^n \prod_{t=1}^{T_i}  f(y_{it}|y_{it-1}, \theta, b_i)p(b_i|\Sigma_b),
\end{aligned}
$$

where $y$ is the stacked vector of $y_{it}$ and there are $n$ individuals. 

A [Gbbs](https://en.wikipedia.org/wiki/Gibbs_sampling) sampling algorithm estimates the parameters by partitioning the joint posterior distribution into conditional distributions. For full details see Appendix C [here](..\papers\longterm_spending.pdf). The R function `gibbs` contained in the R file [dynamic_twopart_mcmc.R](../r/dynamic_twopart_mcmc.R) implements the Gibbs sampler. The function `gibbs` relies on five functions which sample $\alpha$, $\beta$, $\sigma^2_\epsilon$, $b_i$ and $\Sigma_b$. Conjugate priors were chosen for all of the parameters except $b_i$ so sampling straightforward. The conditional distribution of $b_i$ is nonstandard so it is sampled using a random-walk [Metropolis](https://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_algorithm) step.

Before implementing the Gibbs sampler we will simulate data assuming that there is unobserved heterogeneity. 

{% highlight r %}
dat <- sim(sim.T = 5, SIGMA = Sigma)
dat <- dat[year > 0]
ni <- dat[, .N, by = "id"]$N
x1 <- as.matrix(dat[, .(int, c_age, l_d)])
x2 <- as.matrix(dat[, .(int, c_age, l_ly)])
b1 <- dat[, .(b1 = mean(b1)), by = "id"] 
b2 <- dat[, .(b2 = mean(b2)), by = "id"] 
{% endhighlight %}
We must then specify priors and initial values. We need the `MCMCpack` package to sample from an inverse Wishart distribution. 

{% highlight r %}
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
{% endhighlight %}
The Gibbs sampler is run on 10,000 iterations. The sequence is thinned by keeping every 10th draw and the first 5,000 iterations are discarded as burn-in. Since the simulation takes around 17 minutes on my ThinkPad W530 Mobile Workstation, it is not a bad idea to save the output.

{% highlight r %}
source("r/dynamic_twopart_mcmc.R")
gibbs <- Gibbs(nsim = 10000, thin = 10, burn = 5000, y = dat$y, 
               x1 = x1, x2 = x2, ni = ni, id = dat$id, 
               priors = priors, init = inits)
save(gibbs, file = "output/gibbs.RData")
{% endhighlight %}
The simulated posterior densities are returned in list. It is useful to convert the parameter vectors and matrices to a Markov Chain Monte Carlo object using the `coda` package.

{% highlight r %}
library(coda)
alpha.mcmc <- as.mcmc(gibbs$alpha)
beta.mcmc <- as.mcmc(gibbs$beta)
sigma2.mcmc <- as.mcmc(gibbs$sigma2)
Sigma.mcmc <- as.mcmc(gibbs$Sigma) 
{% endhighlight %}
Before summarizing the posterior densities, we must test to see whether our chains have converged. Convergence can be inspected visually with a traceplot, which plots a draw from the posterior density at each iteration against the iteration number. A chain that has converged should have reached a stationary distribution with a relatively constant mean and variance. The parameter values should jump around the posterior density rather than getting stuck in certain regions. To illustrate, consider the traceplot for $\sigma^2_\epsilon$.

{% highlight r %}
plot(sigma2.mcmc, density = FALSE, 
     main = expression("Traceplot of"~sigma[epsilon]^2))
{% endhighlight %}

<img src="/figs/dynamic_twopart_traceplot-1.png" title="plot of chunk dynamic_twopart_traceplot" alt="plot of chunk dynamic_twopart_traceplot" style="display: block; margin: auto;" />
The variance parameter appears to have converged; however, this plot could could be misleading if $\sigma^2_\epsilon$ has only converged to a local region and has not explored the full posterior. It is in general a good idea to run multiple chains with dispersed starting values to ensure that this is not the case. The Gelman-Rubin convergence diagnostic can then be used to test convergence. We will not do that here to keep things simple, but it is good practice.

Now lets look at the posterior quantiles for some of the parameters and compare them to their true values.

{% highlight r %}
cbind(summary(alpha.mcmc)$quant, alpha)
{% endhighlight %}



{% highlight text %}
##             2.5%        25%        50%        75%      97.5% alpha
## int   0.43667047 0.46759930 0.48504669 0.49987207 0.53040524  0.50
## c_age 0.03701827 0.04240177 0.04554813 0.04879016 0.05514034  0.05
## l_d   0.47623616 0.50304348 0.51558136 0.53004923 0.55465324  0.50
{% endhighlight %}



{% highlight r %}
cbind(summary(beta.mcmc)$quant, beta)
{% endhighlight %}



{% highlight text %}
##            2.5%        25%        50%        75%     97.5% beta
## int   5.9433682 5.96312116 5.97560613 5.98680677 6.0071493 6.00
## c_age 0.0871133 0.09194349 0.09420716 0.09663533 0.1014231 0.10
## l_ly  0.2469678 0.24913068 0.25051688 0.25187327 0.2548241 0.25
{% endhighlight %}



{% highlight r %}
c(summary(sigma2.mcmc)$quant, sigma2)
{% endhighlight %}



{% highlight text %}
##      2.5%       25%       50%       75%     97.5%           
## 0.9850542 0.9966137 1.0027020 1.0084922 1.0213431 1.0000000
{% endhighlight %}



{% highlight r %}
cbind(summary(Sigma.mcmc)$quant, c(Sigma))
{% endhighlight %}



{% highlight text %}
##              2.5%       25%       50%       75%     97.5%     
## Sigma11 0.4509504 0.4725421 0.4884947 0.5036067 0.5391468 0.50
## Sigma12 0.2075306 0.2199094 0.2265067 0.2338062 0.2480678 0.25
## Sigma12 0.2075306 0.2199094 0.2265067 0.2338062 0.2480678 0.25
## Sigma22 0.2612393 0.2715018 0.2778633 0.2849412 0.2952004 0.30
{% endhighlight %}
The estimated parameters seem to be converging to their true values, which suggests that the Bayesian algorithm is working as intended. 

