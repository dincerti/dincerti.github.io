# This is the R code for the R Markdown file mortality.Rmd
library("ggplot2")
theme_set(theme_bw())

## @knitr load_data
# life table can be found at: https://www.ssa.gov/oact/STATS/table4c6.html
lt <- read.csv("data/period_lt_2013.csv")

## @knitr probit
lt$qnorm_dprob <- qnorm(lt$death_prob)
lm.dprob <- lm(qnorm_dprob ~ age +  I(age^2) + I(age^3) + sex, data = lt)

## @knitr model_fit
lt$qnorm_dprobhat <- predict(lm.dprob, lt)
ggplot(lt, aes(x = age, y = qnorm_dprob)) + geom_point(size = .2) + 
  geom_line(aes(y = qnorm_dprobhat), col = "blue") +
  facet_wrap(~sex) + xlab("Age") + 
  ylab(expression(Phi^-1*"(death probability)"))

## @knitr predict_prob
lt$dprobhat <- pnorm(lt$qnorm_dprobhat)