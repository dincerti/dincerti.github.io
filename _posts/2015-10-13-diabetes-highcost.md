---
layout: post
title: Predicting High Cost Diabetes Patients using Ridge Regression and Random Forests
---
* TOC
{:toc}


### Overview
Disease management programs provide organized, proactive care to patients with specific conditions in order to help control costs associated with avoidable complications. Since funds are limited, it is important to identify patients that will benefit the most from these programs.

This section uses a few machine learning techniques (logistic regression, ridge regression, and random forests) to predict future high cost diabetes patients.

R code for the analysis can be found [here](https://raw.githubusercontent.com/dincerti/dincerti.github.io/master/_rmd-posts/diabetes_highcost.R), which needs [this]({{site.url}}/data/mepsdiab.csv) dataset.

### Data
We first load the data and examine some summary statistics.

```r
library("data.table")
meps <- fread("data/mepsdiab.csv")
summary(meps$totexp)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##       0    2054    4958   11186   11795  616631
```
The data consists of all diabetes patients from the 10th - 17th panels of the Medical Expenditure Panel Survey (MEPS). The data combines information from the full-year consolidated data files and the medical conditions files.

Mean expenditures are high, although the expenditure distribution is quite skewed. We can examine the skewness of the distribution with a Lorenz curve, which plots the cumulative share of expenditures on the y-axis and the cumulative share of the population on the x-axis.  

```r
library("ineq")
```

```
## Error in library("ineq"): there is no package called 'ineq'
```

```r
library("ggplot2")
theme_set(theme_bw())
lorenz <- Lc(meps$totexp)
```

```
## Error in Lc(meps$totexp): could not find function "Lc"
```

```r
lorenz <- data.table(p = lorenz$p, L = lorenz$L)
```

```
## Error in data.table(p = lorenz$p, L = lorenz$L): object 'lorenz' not found
```

```r
frac <- c(seq(.3, .9, .1), .95)
quants <-  data.table(quant = quantile(lorenz$L, probs = frac),
                      frac)
```

```
## Error in quantile(lorenz$L, probs = frac): object 'lorenz' not found
```

```r
ggplot(lorenz, aes(x = p, y = L)) + geom_line() + 
  geom_abline() + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  xlab("Fraction of individuals ordered by expenditures") + 
  ylab("Cumulative share of expenditures") + theme(legend.title=element_blank()) + 
  scale_x_continuous(breaks = seq(0, 1, by = .1))  + 
  scale_y_continuous(breaks = seq(0, 1, by = .1))  + 
  theme(legend.position = "bottom") +
  geom_text(data = quants, aes(x = frac - .02, y = quant + .02,
                              label = round(quant, 2), size = 8),
            show.legend  = F)
```

```
## Error in ggplot(lorenz, aes(x = p, y = L)): object 'lorenz' not found
```





















