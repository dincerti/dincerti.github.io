---
layout: rmd
title: Markov Cohort Models
---
* TOC
{:toc}


### Overview
This page shows how to employ a Markov cohort model for cost-effectiveness analysis using R. The analysis replicates a paper ([Chanc97](references.html#Chanc97)) that is commonly used for teaching purposes ([BSC06](references.html#BSC06), [Drum05](references.html#Drum05)), which compares combination HIV therapy using zidovudine and lamivudine to zidovudine monotherapy.

To replicate this page you will need [this](r/markov_cohort.R) R script. 

### Setting Up a Markov Model
In a Markov model patients move from one mutually exclusive disease state to another over a series of discrete time periods. A fully specified model consists of 1) disease states, 2) the probability of transitioning from one stage to the next, 3) treatment costs in each stage, and 4) quality of life weights by stage. 

**States**

In our example, HIV patients can be in one of the following 4 states at a given point in time:

* State A: 200 < cd4 < 500 cells/mm$$^3$$,
* State B: cd4 < 200 cells/mm$$^3$$,
* State C: Aids,
* State D: Death.

**Transition Matrix**

At any given points in time, the row vector $$z_{kt} = [z_{1kt}\; z_{2kt}\; z_{3kt}\; z_{4kt}]$$ represents the total number of HIV patients in each of the 4 states using treatment $$k$$ at time $$t$$. This vector, $$z_{kt}$$, changes over time according to a transition matrix, $$P_{kt}$$,

$$
\begin{aligned}
\underbrace{z_{kt}}_{1 \times 4} &= \underbrace{z_{k,t-1}}_{1 \times 4} \underbrace{P_{kt}}_{4 \times 4},
\end{aligned}
$$

where

$$
\begin{aligned}
P_{kt} &=
\begin{bmatrix}
p_{11} & p_{12} & p_{13} & p_{14}\\
p_{21} & p_{22} & p_{23} & p_{24}\\
p_{31} & p_{32} & p_{33} & p_{34}\\
p_{41} & p_{42} & p_{43} & p_{44}
\end{bmatrix}.
\end{aligned}
$$

The probability of transitioning from state $$i$$ to state $$j$$ at time $$t$$ is denoted by $$p_{ij}$$. Under monotherapy (treatment $$0$$), the estimated transition matrix is constant and given by,

$$
\begin{aligned}
P_0 &=
\begin{bmatrix}
0.721 & 0.202 & 0.067 & 0.010\\
0 & 0.581 & 0.407 & 0.012\\
0 & 0 & 0.750 & 0.250\\
0 & 0 & 0 & 1
\end{bmatrix}.
\end{aligned}
$$

Transition probabilities for combination therapy (treatment $$1$$) are based on an estimated relative risk of disease progression of $$0.509$$ from a meta-analysis of 4 comparative trials. The relative risk is assumed to reduce the probability of transitioning to a worse state (i.e. $$p_{12}, p_{13}, p_{14}, p_{23}, p_{24}, p_{34}$$) by a factor of .509.

We set up the transition probabilities in R.


{% highlight r %}
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
{% endhighlight %}

**Costs**

Treatment costs are due to a) direct medical and community expenses and b) the prices of each drug. The drug costs of zidovudine and lamivudine are &pound;2278 and &pound;2086 respectively. Medical and community expenses in each state are as follows:

* State A: &pound;2256
* State B: &pound;3052
* State C: &pound;9007
* State D: &pound;0

In R, we have:


{% highlight r %}
c.zidovudine <- 2278
c.lamivudine <- 2086.50
c.0 <- c(2756 + c.zidovudine, 3052 + c.zidovudine, 9007 + c.zidovudine, 0)
c.1 <- c(c.0[1:3] + c.lamivudine, 0)
{% endhighlight %}


**Quality of Life Weights**

Treatment effects are typically measured by quality-adjusted life-years (QALYs). However, we follow the original paper which measured effectiveness with (unadjusted) life-years. In mathematical terms, this means that the 4 states are weighted with the row vector $$[1\; 1\; 1\; 0]$$. 

{% highlight r %}
qolw <- c(1, 1, 1, 0)
{% endhighlight %}
### Cohort Simulation in R
A cohort simulation uses the Makov model to measure the experiences of a hypothetical cohort, say 1000 patients, over a set period of time (i.e. 20 years), under each treatment option. 

We load an R function that simulates the costs and effects of an intervention given a transition matrix ```P```; an initial row vector ```z0``` containing the number of patients in each state at time $$0$$; the number of cycles ```ncycles```, the costs in each state; the quality of life weights; and a discount factor for future costs. (Note that, like the original paper, the discount factor is only applied to costs, although a discount could be applied to life-years as well.) 


{% highlight r %}
source("r/markov.R")
{% endhighlight %}


{% highlight r %}
MarkovCohort <- function(P, z0, ncycles, costs, qolw, discount){
  # Calculates quantities for cost-effectiveness analysis using a markov
  # cohort model.
  #
  # Args:
  #   P: State transition matrix.
  #   z0: Initial values of z
  #   ncycles: number of cycles
  #   costs: vector containing costs in each state. May be a time 
  #         constant vector or time varying list of vectors.
  #   qolw: vector containing quality of life weight in each state.
  #           May be a time constant vector or time varying list of
  #           vectors
  #   discount: cycle discount rate
  #
  # Returns:
  #   List containing proportion in each state by cycle, costs, and effects.
  if (is.list(P)){
    nc <- ncol(P[[1]])
  } else {
    nc <- ncol(P)
  }
  z <- rbind(z0, matrix(NA, nrow = ncycles, ncol = nc))
  c <- e <- rep(NA, ncycles)
  delta <- 1/(1 + discount)^(seq(0, ncycles))
  for (t in 2:(ncycles+1)){
    # z
    if (is.list(P) == TRUE) {
      z[t, ] <- z[t-1, , drop = FALSE] %*%  P[[t-1]] 
    } else{
      z[t, ] <- z[t-1, , drop = FALSE] %*%  P
    }
    
    # costs
    if (is.list(costs)){
      c[t] <- delta[t] * z[t, ] %*% costs[[t-1]]/sum(z0)
    } else{
      c[t] <- delta[t] * z[t, ] %*% costs/sum(z0)
    }
    
    # effectiveness
    if (is.list(effects)){
      e[t] <- z[t, ] %*% qolw[[t-1]]/sum(z0)
    } else{
      e[t] <- z[t, ] %*% qolw/sum(z0)
    }
  }
  return(list(z = z, c = c[-1], e = e[-1]))
}
{% endhighlight %}
At each cycle, the function calculates the number of patients in each state ```z[t, ]```; total costs per patient, ```c[t]```; and life-years ```e[t]```. The algorithm is very simple but appears more complicated because it allows for both time constant and time-varying transition matrices and costs. The time-varying transition matrices are necessary because lamivudine is only assumed to be given for the first two-years of treatment. 

Using the function we simulate outcomes under both treatments.

{% highlight r %}
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
{% endhighlight %}

### Decision Analysis
Armed with the simulation results, we can plot survival curves by treatment.

{% highlight r %}
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
{% endhighlight %}

<img src="/figs/survival-1.png" title="plot of chunk survival" alt="plot of chunk survival" style="display: block; margin: auto;" />
Due to the estimated relative risk of disease progression, the survival curve for combination therapy lies above the curve for monotherapy. That said, costs are higher for combination therapy because of 1) additional drug costs for lamivudine and 2) patients being treated longer due to increased survival.

To summarize the cost-effectiveness of the intervention we can use the incremental cost-effectiveness ratio (ICER), or,

$$
\begin{aligned}
\frac{c_1 - c_0}{e_1 - e_0},
\end{aligned}
$$

where $$c_k$$ and $$e_k$$ refer to the costs and effects under treatment $$k$$ respectively. In R,

{% highlight r %}
icer <- (sum(sim1$c) - sum(sim0$c))/(sum(sim1$e) - sum(sim0$e))
print(icer)
{% endhighlight %}



{% highlight text %}
## [1] 6276.083
{% endhighlight %}
which suggest that a decision-maker should only approve combination therapy if he or she is willing to pay over ``6276`` pounds for an additional life-year gained.
