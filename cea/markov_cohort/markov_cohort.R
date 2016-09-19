ncycles <- 500
z0 <- c(1000, 0)

# 1000 simulations
ptm <- proc.time()
sims <- 10000
tmp <- markovCohortC(z0, ncycles, discount = .03, nsims = sims,
                             P = array(c(.95, 0, .05, 1), c(4, sims, 1)),
                             costs = array(c(25, 50), c(2, sims, 1)),
                             qol = array(c(1, 0), c(2, sims, 1)),
                             P_indx = rep(1, ncycles), cost_indx = rep(1, ncycles),
                             qol_indx = rep(1, ncycles))
proc.time() - ptm

# single simulation
tmp2 <- markovCohort(z0 = z0, ncycles = ncycles, discount = .03, nsims = 1,
             P = c(.95, 0, .05, 1), costs = c(25, 50), qol = c(1, 0))

# tornado plot
p.dat <- data.frame(par = c("Cost in state 1", "Utility in state 1",
                            "Cost in state 1", "Utility in state 1"),
                    par.values = c(.7, .6, .3, .2),
                    labels = c("High", "High", "Low", "Low"),
                    output.values = c(150, 100, -100, -120))
ggplotTornado(data = p.dat, par = "par", par.values = "par.values",
              labels = "labels", output.values = "output.values",
                          output.primary = 50, output.label = "ICER")

