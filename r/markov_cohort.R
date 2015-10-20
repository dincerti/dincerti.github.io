
# This is the R code for the R Markdown file markov_cohort.Rmd
setwd("C:/Users/Devin/Dropbox/Projects/dincerti.github.io")
rm(list=ls()) 

## ---- PARAMETERS -------------------------------------------------------------
## @knitr transition_matrices
P.0 <- matrix(c(.721, .202, .067, .01, 
              0, .581, .407, .012,
              0, 0, .75, .25,
              0, 0, 0, 1),
            ncol = 4, nrow = 4, byrow = TRUE)
P.1 <- matrix(c(.858, .103, .034, .005,
             0, .787, .207, .006,
             0, 0, .873, .127,
             0, 0, 0, 1),
             ncol = 4, nrow = 4, byrow = TRUE)

## @knitr costs
c.zidovudine <- 2278
c.lamivudine <- 2086.50
c.0 <- c(2756 + c.zidovudine, 3052 + c.zidovudine, 9007 + c.zidovudine, 0)
c.1 <- c(c.0[1:3] + c.lamivudine, 0)

## @knitr effects
qolw <- c(1, 1, 1, 0)

## ---- COHORT SIMULATION ------------------------------------------------------
## @knitr markov_cohort
source("r/markov.R")

## @knitr simulation
ncycles <- 20
sim0 <- MarkovCohort(P = P.0,  z0 = c(1000, 0, 0, 0), ncycles = ncycles,
                    costs = c.0, qolw = qolw, 
                    discount = 0.06)
sim1 <- MarkovCohort(P = c(replicate(2, P.1, simplify = FALSE), 
                              replicate(ncycles - 2, P.0, simplify = FALSE)),
                     z0 = c(1000, 0, 0, 0), ncycles = ncycles,
                     costs = c(replicate(2, c.1, simplify = FALSE),
                               replicate(ncycles - 2, c.0, simplify = FALSE)),
                     qolw = qolw, discount = 0.06)

## ---- DECISION ANALYSIS ------------------------------------------------------
## @knitr survival
library("ggplot2")
theme_set(theme_bw())
surv.df <- data.frame(surv = c(sim0$e, sim1$e),
                      cylce = rep(seq(1, ncycles), 2),
                       lab = rep(c("Monotherapy", "Combination therapy"), 
                       each = ncycles))
ggplot(dat = surv.df, aes(x = cylce, y = surv, col = lab)) + geom_line() + 
  xlab("Cycle") + ylab("Fraction surviving") +
  theme(legend.title=element_blank()) + 
  theme(legend.position = "bottom")


## @knitr icer
icer <- (sum(sim1$c) - sum(sim0$c))/(sum(sim1$e) - sum(sim0$e))
print(icer)


