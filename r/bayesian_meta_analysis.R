
# This is the R code for the R Markdown file bayesian_meta_analysis.Rmd
setwd("C:/Users/Devin/Dropbox/Projects/dincerti.github.io")

## ---- PREVIOUS RCTS ----------------------------------------------------------
## @knitr rcts
rct <- data.frame(study = c("New York", "Malamo I", "Kopparberg", "Ostergotland",
              "Canada I", "Canada II", "Stockholm", "Goteborg", "UK age trial"))
rct$year <- c(1963, 1976, 1977, 1978, 1980, 1980, 1981, 1982, 1991)
rct$d1 <- c(218, 87, 126, 135, 105, 107, 66, 88, 105)
rct$n1 <- c(31000, 20695, 38589, 38491, 25214, 19711, 40318, 21650, 53884)
rct$d0 <- c(262, 108, 104, 173, 108, 105, 45, 162, 251)
rct$n0 <- c(31000, 20783, 18582, 37403, 25216, 19694, 19943, 29961, 106956)

## @knitr relative_risk
rct$p1 <- rct$d1/rct$n1
rct$p0 <- rct$d0/rct$n0
rct$rr <- rct$p1/rct$p0
rct$lrr <- log(rct$rr)
rct$lse <- sqrt((1 - rct$p1)/(rct$p1 * rct$n1) + (1 - rct$p0)/(rct$p0 * rct$n0))
rct$lower <- exp(rct$lrr - qnorm(.975) * rct$lse)
rct$upper <- exp(rct$lrr + qnorm(.975) * rct$lse)

## @knitr forest_plot
library("metafor")
p <- forest(x = rct$rr, ci.lb = rct$lower, ci.ub = rct$upper, 
       slab = paste(rct$study, rct$year, sep = ", "), refline = 1)
text(min(p$xlim), .88 * max(p$ylim), "Study and Year", pos = 4, font = 2)
text(max(p$xlim), .88 * max(p$ylim), "Relative Risk [95% CI]", pos = 2, font = 2)

## ---- COMPLETE POOLING -------------------------------------------------------
## @knitr fixed_effects
me.fe <- rma(rct$lrr, rct$lse^2, method = "FE")
c(exp(me.fe$b), exp(me.fe$ci.lb), exp(me.fe$ci.ub))

## @knitr weighed_mean
exp(weighted.mean(rct$lrr, 1/(rct$lse^2)))

## ---- MAXIMUM LIKELIHOOD HIERARCHICAL MODEL ----------------------------------
## @knitr random_effects
me.re <- rma(rct$lrr, rct$lse^2)
c(exp(me.re$b), exp(me.re$ci.lb), exp(me.re$ci.ub))

## ---- BAYESIAN META ANALYSIS -------------------------------------------------
## @knitr stan_data
library("rstan")
J <- nrow(rct)
stan.dat <- list(J = J, y = rct$lrr, sigma = rct$lse)

## @knitr stan_fit
fit <- stan(file = "stan/bayesian_meta_analysis.stan",
            data = stan.dat, iter = 2000, chains = 4)
post <- extract(fit, permuted = TRUE)

## @knitr ci_mu
quantile(exp(post$mu), probs = c(.025, .5, .975))

## @knitr ci_tau
quantile(post$tau, probs = c(.025, .5, .975))

## @knitr theta_plot
library("ggplot2")
theme_set(theme_bw())
p.dat <- apply(exp(post$theta), 2, quantile, probs = c(.025, .5, .975))
p.dat <- data.frame(lower = p.dat[1, ], rr = p.dat[2, ], upper = p.dat[3, ])
p.dat <- rbind(p.dat, rct[, c("lower", "upper", "rr")])
p.dat$lab <- rep(c("Theta", "Y"), each = J)
p.dat$id <- rep(seq(9, 1), 2)
p.dat$idlab <- factor(p.dat$id, labels = rev(paste(rct$study, rct$year, sep = ", ")))
ggplot(p.dat, aes(x = idlab, y = rr, ymin = lower, ymax = upper, col = lab)) +  
  geom_pointrange(aes(col = lab), position = position_dodge(width = 0.50)) +
  coord_flip() + geom_hline(aes(yintercept = mean(exp(post$mu))), lty = 2) +  xlab("") + 
  ylab("")  + theme(legend.position="bottom") + 
  scale_colour_discrete(name="", 
                        labels = c("Theta" = bquote("Random effect:"~exp(theta[J])~" "),
                                    "Y"= bquote("Relative risk:"~exp(Y[J]))))

## @knitr prediction
n.sims <- nrow(post$mu)
theta.new <- rep(NA, n.sims)
for (i in 1:n.sims){ 
  theta.new[i]  <- rnorm(1,  post$mu[i],  post$tau[i]) 
}

## @knitr prediction_ci
quantile(exp(theta.new), probs = c(.025, .5, .975))




