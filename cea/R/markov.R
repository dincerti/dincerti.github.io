
# this file contains functions for markov cost-effectiveness models

## ---- MARKOV COHORT MODEL ----------------------------------------------------
## @knitr markov_cohort_func
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