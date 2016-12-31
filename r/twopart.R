# This is the R code for the R Markdown file twopart.Rmd

## ---- HEALTH EXPENDITURE DATA FROM THE MEPS ----------------------------------
## @knitr meps_load
meps <- readRDS("data/meps-fyc-2012.rds")

## @knitr twopart_hist
library(data.table)
library(ggplot2)
theme_set(theme_bw())
qplot(totexp12, data = meps[totexp12 < quantile(meps$totexp12, probs = .95)])

## @knitr twopart_loghist
library(MASS)
meps <- meps[, logtotexp12 := ifelse(totexp12 > 0, log(totexp12), NA)]
ggplot(meps, aes(x = logtotexp12)) + 
  geom_histogram(aes(y = ..density..), colour = "black", fill = "White") +
  stat_function(fun = dnorm, args = fitdistr(meps[totexp12 > 0 , logtotexp12],
                                             "normal")$estimate)

## ---- PREDICTING MEAN EXPENDITURES -------------------------------------------
## @knitr meps_vars
meps <- meps[, age := ifelse(age12x < 0, NA, age12x)]
meps <- meps[, age2 := age^2]
meps <- meps[, srh := as.factor(ifelse(rthlth53 < 0, NA, rthlth53))] # health status
meps <- meps[, hisp := ifelse(racethx == 1, 1, 0)] # hispanic race
meps <- meps[, black := ifelse(racebx == 1 | racebx == 2, 1, 0)]
meps <- meps[, prvins := ifelse(inscov12 == 1, 1, 0)] # private insurance
meps <- meps[, pubins := ifelse(inscov12 == 2, 1, 0)] # public insurance
meps <- meps[, d_totexp12 := ifelse(totexp12 == 0, 0, 1)] # indicator for positive spending

## @knitr meps_subset
meps <- meps[, id := seq(1, nrow(meps))]
xvars <- c("age", "age2", "srh", "hisp", "black", "prvins", "pubins")
meps <- meps[, c("id", "totexp12", "d_totexp12", "logtotexp12", xvars), with = FALSE]
meps <- meps[complete.cases(meps[, xvars, with = FALSE])]

## @knitr meps_sample
set.seed(100)
meps <- meps[, sample := sample(c("train", "test"), nrow(meps), replace = TRUE)]

## @knitr fit
Fm <- function(y, xvars){
  return(as.formula(paste(y, "~", paste(xvars, collapse = "+"))))
}
# part 1
logistic.fit <- glm(Fm("d_totexp12", xvars), meps[sample == "train"], family = binomial)

# part 2
ols.fit <- lm(Fm("totexp12", xvars), meps[totexp12 > 0 & sample == "train"])
logols.fit <- lm(Fm("logtotexp12", xvars), meps[totexp12 > 0 & sample == "train"])
gamma.fit <- glm(Fm("totexp12", xvars), meps[totexp12 > 0 & sample == "train"], 
                 family = Gamma(link = log))

## @knitr pred
phat <- predict(logistic.fit, meps[sample == "test"], type = "response")
pred <- data.table(totexp12 = meps[sample == "test", totexp12])
pred$ols <- phat * predict(ols.fit, meps[sample == "test"])
pred$logols <- phat * exp(predict(logols.fit, meps[sample == "test"]) + summary(logols.fit)$sigma^2/2)
pred$gamma <- phat * predict(gamma.fit, meps[sample == "test"], type = "response")

## @knitr rmse
RMSE <- function(x, y)  sqrt(mean((y - x)^2, na.rm = TRUE))
rmse <- c(RMSE(pred$totexp12, pred$ols),
          RMSE(pred$totexp12, pred$logols),
          RMSE(pred$totexp12, pred$gamma))
names(rmse) <- c("OLS", "Log OLS", "Gamma")
print(rmse)

## @knitr smearing
meps <- meps[, agecat := cut(age, breaks = c(0, 1, seq(5, 90, 5)), 
                             right = FALSE)]
epsilon <- data.table(age = logols.fit$mode$age, res = logols.fit$res)
epsilon <- epsilon[, agecat := cut(age, breaks = c(0, 1, seq(5, 90, 5)), 
                                   right = FALSE)]
epsilon <- epsilon[, .(phihat = mean(exp(res))), by = "agecat"]
meps <- merge(meps, epsilon, by = "agecat", all.x = TRUE)
meps <- meps[order(id)]
pred$logols_smear <- phat * exp(predict(logols.fit, meps[sample == "test"])) * mean(exp(logols.fit$res))
pred$logols_hetsmear <- phat * exp(predict(logols.fit, meps[sample == "test"])) * meps[sample == "test", phihat]
rmse <- c(RMSE(pred$totexp12, pred$logols_smear),
          RMSE(pred$totexp12, pred$logols_hetsmear))
names(rmse) <- c("Log OLS Homoskedastic Smearing", "Log OLS Heteroskedastic Smearing")
print(rmse)

## ---- PREDICTIVE SIMULATION --------------------------------------------------
## @knitr logistic_normal
n <- nrow(meps[sample == "test"])
d <- rbinom(n, 1, phat)
y.norm <- d * rnorm(n, pred$ols, summary(ols.fit)$sigma)

## @knitr logistic_lognormal
y.lognorm <- d * rlnorm(n, predict(logols.fit, meps[sample == "test"]) , 
                        summary(logols.fit)$sigma)

## @knitr mom_dispersion
res <- (gamma.fit$model$totexp12 - gamma.fit$fit)/gamma.fit$fit # this is equivalent to gamma.fit$res
c(sum(res^2)/gamma.fit$df.res, summary(gamma.fit)$dispersion)

## @knitr logistic_gamma
a <- gamma.shape(gamma.fit)$alpha
b <- a/pred$gamma
y.gamma <- d * rgamma(n, shape = a , rate = b)

## @knitr twopart_yrepden
y <- meps[sample == "test", totexp12]
p.dat <- data.table(y = c(y, y.norm, y.lognorm, y.gamma),
                    lab = c(rep("Observed", n), rep("Normal", n), 
                            rep("Lognormal", n), rep("Gamma", n)))
p <- ggplot(p.dat[y > 0 & y < 10000], aes(x = y, col = lab)) + 
  geom_density(kernel = "gaussian") +
  xlab("Expenditures") + ylab("Density") +
  theme(legend.position="bottom") + labs(col = "") +
  scale_color_manual(values=c(Observed = "black", Normal = "red", 
                              Lognormal = "blue", Gamma = "green")) 
print(p)

## @knitr yrep_quantiles
MySum <- function(x){
  q <- c(0.30, 0.5, 0.75, .9, .95, .98)
  dat <- c(100 * mean(x == 0, na.rm = TRUE),
           min(x, na.rm = TRUE), quantile(x, probs = q, na.rm = TRUE), 
           max(x, na.rm = TRUE))
  names(dat) <- c("PercentZero", "Min", paste0("Q", 100 * q), "Max")
  return(round(dat, 0))
} 
sumstats <- rbind(MySum(y), MySum(y.norm), 
                  MySum(y.lognorm), MySum(y.gamma))
rownames(sumstats) <- c("Observed", "Normal", "Lognormal", "Gamma")
print(sumstats)
