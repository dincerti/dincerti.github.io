# This is the R code for the R Markdown file markov_cohort.Rmd
library("cea")
rm(list=ls()) 
mlogit.beta <- readRDS("cea/output/mlogit_bayes.rds")

## ---- TRANSITION PROBABILITIES -----------------------------------------------
## @knitr transprob_alive_setup
ages <- seq(45, 100, 1)
x <- matrix(c(rep(1, length(ages)), ages), ncol = 2)
x <- x[rep(seq_len(nrow(x)), times = 2), ]
x <- cbind(x, rep(c(0, 1), each = length(ages)))
colnames(x) <- c("int", "age", "male")
n.x <- nrow(x)

## @knitr transprob_alive
mlogit_prob <- function(x, beta) {
  log.odds <- as.matrix(x) %*% beta
  odds <- cbind(exp(0), exp(log.odds))
  odds.sum <- apply(odds, 1, sum) 
  p <- odds/odds.sum
  return(p)
}
nsims <- nrow(mlogit.beta[[1]])
p <- matrix(NA, nrow = nsims * n.x, ncol = 3 * 3)
id <- data.frame(sim = rep(seq(1, nsims), each = n.x),
                 age = rep(x[, "age"], times = nsims),
                 male = rep(x[, "male"], times = nsims))
id$cycle <- id$age - ages[1]
j1 <- 1; j2 <- 3
for (j in 1:3){
  for (i in 1:nsims){
    betaij <- matrix(mlogit.beta[[j]][i, ], ncol = 2, byrow = T)
    indices <- ((i - 1) * nrow(x) + 1):(i * nrow(x))
    p[indices, j1:j2] <- mlogit_prob(x, betaij)
  }
  j1 <- j1 + 3; j2 <- j2 + 3
}

## @knitr adjust_death
lt <- read.csv("data/period_lt_2013.csv")
lt$male <- ifelse(lt$sex == "male", 1, 0)
dat <- merge(dat, lt[, c("age", "male", "death_prob")], by = c("age", "male"))
dat$death <- rbinom(nrow(dat), 1, dat$death_prob)
dat$trans2 <- ifelse(dat$death == 1, 4, dat$trans)

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
z0 <- matrix(c(1000, 0, 0, 0), nrow = 1)
sim0 <- markov_trans(z0, ncycles, c(t(P.0)), nsims = 1)
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


