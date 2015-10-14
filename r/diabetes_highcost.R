
# This is the R code for the R Markdown file diabetes_highcost.Rmd
setwd("C:/Users/Devin/Dropbox/Projects/dincerti.github.io")

## ---- DATA -------------------------------------------------------------------
## @knitr data
library("data.table")
meps <- fread("data/mepsdiab.csv")
summary(meps$totexp)

## @knitr lorenz_curve
library("ineq")
library("ggplot2")
theme_set(theme_bw())
lorenz <- Lc(meps$totexp)
lorenz <- data.table(p = lorenz$p, L = lorenz$L)
frac <- c(seq(.3, .9, .1), .95)
quants <-  data.table(quant = quantile(lorenz$L, probs = frac),
                      frac)
ggplot(lorenz, aes(x = p, y = L)) + geom_line() + 
  geom_abline() + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  xlab("Fraction of individuals ordered by expenditures") + 
  ylab("Cumulative share of expenditures") + theme(legend.title=element_blank()) + 
  scale_x_continuous(breaks = seq(0, 1, by = .1))  + 
  scale_y_continuous(breaks = seq(0, 1, by = .1))  + 
  theme(legend.position = "bottom") +
  geom_text(data = quants, aes(x = frac - .02, y = quant + .02,
                              label = round(quant, 2), size = 8),
            show_guide  = F)

## ---- TRAINING AND TEST DATA -------------------------------------------------
## @knitr variables
ccsnames <- names(meps)[grep("^ccs\\d\\d\\d", names(meps))]
ds <- c("dsdiet", "dsmed", "dsinsu", "dseypr", "dskidn")
xvars <- c("l_topexp", "l_totexp", "age", "srh", ds, ccsnames)
meps <- meps[, c("dupersid", "panel", "year", "topexp", "totexp",
                 xvars), with = FALSE]

## @knitr datasets
train <- meps[panel <= 15]
train <- train[complete.cases(train)]
test <- meps[panel > 15]
test <- test[complete.cases(test)]

## ---- LOGISTIC REGRESSION ----------------------------------------------------
## @knitr logistic_regression
Fm <- function(y, x){
  return(as.formula(paste(y, "~", paste(xvars, collapse = "+"))))
}
logistic <- glm(Fm("topexp", xvars), family = binomial, train)
p.logistic <- predict.glm(logistic, test, type = "response")

## ---- REGULARIZED LOGISTIC REGRESSION ----------------------------------------
## @knitr logistic_ridge_regression
library("glmnet")
x.train <- model.matrix(Fm("topexp", xvars), train)[, -1]
y.train <- train$topexp
x.test <- model.matrix(Fm("topexp", xvars), test)[, -1]
y.test <- test$topexp
cv.ridge <- cv.glmnet(x.train, y.train, family = "binomial", alpha = 0)
p.ridge <- predict(cv.ridge, x.test, s = "lambda.min", type = "response")

## ---- RANDOM FORESTS ---------------------------------------------------------
## @knitr rf
library("randomForest")
set.seed(101)
train$srh <- as.numeric(train$srh)
x.train <- model.matrix(Fm("topexp", xvars), train)[, -1]
y.train <- as.factor(y.train)
rf <- randomForest(x = x.train, y = y.train, ntree = 501,
                   importance = TRUE, do.trace= FALSE)
test$srh <- as.numeric(test$srh)
p.rf <- predict(rf, test, type = "prob")

## ---- CLASSIFIER PERFORMANCE -------------------------------------------------
## @knitr brier
Brier <- function(p, y) mean((p - y)^2)
c(Brier(p.logistic, y.test), Brier(p.ridge, y.test), 
  Brier(p.rf[, 2], y.test))

## @knitr roc
library("ROCR")

# predictions
p.mat <- cbind(p1 = p.logistic, p2 = p.ridge[, 1], p3 = p.rf[, 2])
labels.mat <- matrix(y.test, nrow = length(y.test), ncol = ncol(preds))
pred <- prediction(p.mat, labels.mat)

# roc curve
perf <- performance(pred, "tpr", "fpr")
roc.dat <- data.frame(fpr = do.call("c", perf@x.values), 
                       tpr = do.call("c", perf@y.values))
len <- do.call("c", lapply(perf@x.values, length))
roc.dat$lab <- rep(c("Logistic regression", "Logistic ridge regression", 
                     "Random forest"), times = len)
ggplot(roc.dat, aes(x = fpr, y = tpr, col = lab)) + geom_line() +
  xlab("False positive rate") +
  ylab("True positive rate") + theme(legend.title=element_blank()) + 
  theme(legend.position = "bottom")

## @knitr auc
auc <- as.numeric(performance(pred, "auc")@y.values)
print(auc)

## @knitr ppv
ppv.dat <- data.frame(pred = p.ridge[, 1], y = y.test, exp = test$totexp)
ppv.dat$ecdf.p <- ecdf(ppv.dat$pred)(ppv.dat$pred)
pctile <- seq(.7, 1, .01)
tp <- rep(NA, length(pctile))
for (i in 1:length(tp)){
  yhat <-  1 * (ppv.dat$ecdf.p > pctile[i])
  tp[i] <- mean(ppv.dat$y[which(yhat == 1)])
}
ggplot(data.frame(x = pctile, y = tp), aes(x = x, y = y)) + geom_line() +
  xlab("Percentile cutoff") + ylab("Positive predicted value") + 
  geom_vline(x = 0.9, lty = "dotted", col = "red")

## @knitr spending_ratio
fsaf <- sum(ppv.dat$exp[ppv.dat$ecdf.p >= .9])/sum(ppv.dat$exp[ppv.dat$y == 1])
print(fsaf)


