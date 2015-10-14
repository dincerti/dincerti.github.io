
# This is the R code for the R Markdown file markov_contraception.Rmd
setwd("C:/Users/Devin/Dropbox/Projects/dincerti.github.io")

## ---- MODEL SET UP -----------------------------------------------------------
p1.1 <- c(.992, 1, 1, 1, 1)
p1.2 <- .008 * .03
p1.3 <- .008 * .4462
p1.4 <- .008 * .1649
p1.5 <- .008 * .3589
pi <- matrix(data