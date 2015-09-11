
# UPDATE ALPHA -----------------------------------------------------------------
UpdateAlpha <- function(Y = y, X = x, alpha, V = Valpha,
                       b1 = b[, 1], m0, V0inv) {
  obs <- length(Y)
  dstar <- rep(0, obs)
  dstar.mean <- X %*% alpha + rep(b1, times = ni)		# data augmentation with current alpha
  dstar[Y == 0] <- qnorm(runif(obs, 0, pnorm(0, dstar.mean, 1)), dstar.mean, 1)[Y == 0] # left truncated normal
  dstar[Y > 0] <- qnorm(runif(obs, pnorm(0, dstar.mean, 1), 1), dstar.mean, 1)[Y > 0]   #  right truncated normal
	# Valpha was calculated outside of this function to avoid repitition
  malpha <- V %*% (V0inv %*% m0 + crossprod(X, dstar - rep(b1, times = ni)))
  return(c(rmvnorm(1, malpha, V)))
}

# UPDATE BETA ------------------------------------------------------------------
UpdateBeta <-function(Y = y, logY = ly, X, b2 = b[, 2], sigma2inv = sigma2inv, 
                      m0, V0inv) {
  Vbeta <- solve(V0inv + sigma2inv * crossprod(X[Y > 0,], X[Y > 0,]))
  mbeta <- Vbeta %*% (V0inv %*% m0 + sigma2inv * t(X[Y > 0, ]) %*% (logY - rep(b2, times = ni))[Y>0]) 
  return(c(rmvnorm(1, mbeta, Vbeta)))
}

# UPDATE SIGMA2inv -------------------------------------------------------------
UpdateSigma2inv <- function(Y = y, X, logY = yl, prior.shape, prior.rate, b, beta){
  shape <- prior.shape + length(Y[Y>0])/2      
  b2 <- rep(b[, 2], times = ni) 
  rate <- prior.rate +  sum((logY- X %*% beta - b2)[Y > 0, ]^2)/2       	 
  return(rgamma(1, shape = shape, rate = rate))
}

# UPDATE SIGMA (COVARIANCE MATRIX FOR b_i) -------------------------------------
UpdateSigma <- function(prior.v, prior.S, b){
  return(riwish(prior.v + nrow(b), prior.S + crossprod(b, b)))
}

# UPDATE b_i (CORRELATED RANDOM EFFECTS) ---------------------------------------
Updateb <- function(Y = y, X1, X2, b, ni, id, Sigma.prop = diag(2), Sigma,
                    alpha, beta, sigma2 = 1/sigma2inv) {
  l <- l.prop <- rep(NA, length(Y))
  n <- nrow(b)
  
  # current log likelihood
  dstar.mean <- X1 %*% alpha + rep(b[, 1], times = ni)
  p <- pnorm(dstar.mean)
  ly.mean <- (X2 %*% beta + rep(b[, 2], times = ni))[Y > 0]
  l[Y == 0] <- log(1- p[Y == 0])		
  l[Y > 0] <- log(p[Y > 0] * dlnorm(Y[Y > 0], ly.mean, sqrt(sigma2)))
  
  # random walk proposal 
  b.prop <- b + rmvnorm(n, sigma = Sigma.prop)  		# or for student t: b + rmvt(sigma = Sigma.prop, 3)
  dstar.mean.prop <- X1 %*% alpha + rep(b.prop[, 1], times = ni)
  p.prop <- pnorm(dstar.mean.prop)
  ly.mean.prop <- (X2 %*% beta + rep(b.prop[, 2], times = ni))[Y > 0]  
  l.prop[Y == 0] <- log(1 - p.prop[Y == 0])			
  l.prop[Y > 0]<- log(p.prop[Y > 0] * dlnorm(Y[Y > 0], ly.mean.prop, sqrt(sigma2)))
  
  # compare likelihoods
  L.prop <- data.table(l = l.prop, id = id)
  L.prop <- L.prop[, .(l = sum(l)), by = "id"][, l]
  L <- data.table(l = l, id = id)
  L <- L[, .(l = sum(l)), by = "id"][, l]
  ratio <- L.prop + dmvnorm(b.prop, c(0,0), Sigma, log = TRUE) - # note: adding pr(b_i) to likelihood
         (L + dmvnorm(b, c(0, 0), Sigma, log = TRUE))
  
  # metropolis jumping rule
  u <- 1 * (log(runif(n)) < ratio)
  b[u == 1, ] <- b.prop[u==1, ]
  return(list(b = b, acc = u))
}

# GIBBS SAMPLER ----------------------------------------------------------------
Gibbs <- function(nsim = 10000, thin = 5, burn, y, x1, x2, ni, id, 
                  priors = list(malpha, Valpha, mbeta, Vbeta, sigma2.shape, sigma2.rate),
                  init = list(alpha, beta, sigma2, b, Sigma)){
    
  # simulation parameters
  if(missing(burn)) burn <- nsim * 0.5
  lastit <- (nsim-burn)/thin  # Last stored value
  
  # initial values
  alpha <- init$alpha
  beta <- init$beta
  sigma2 <- init$sigma2
  b <- init$b
  Sigma <- init$Sigma
  
  # data
  N <- length(y)
  n <- nrow(b)
  k1 <- ncol(x1)
  k2 <- ncol(x2)
  ly <- rep(0, N)
  ly[y>0] <- log(y[y > 0]) # logy for y > 0

  # variances and precisions 
  prior.Vinvalpha <- solve(priors$Valpha)
  prior.Vinvbeta <- solve(priors$Vbeta)
  sigma2inv <- 1/sigma2
  Valpha <- solve(prior.Vinvalpha+ crossprod(x1, x1)) # variance for probit regression
  
  # storage objects
  Alpha <- matrix(NA, lastit, k1)  
  colnames(Alpha) <- colnames(x1)
  Beta <- matrix(NA, lastit, k2) 
  colnames(Beta) <- colnames(x2)
  Sigma2 <- rep(NA, lastit)    # Lognormal SDs
  b.acc <- rep(0, n) # set initial acceptance rates at zero 
  B.acc <- matrix(NA, lastit, n)
  B1 <- B2 <- matrix(NA, lastit, n)
  SIGMA <- matrix(NA, lastit, 4)
  colnames(SIGMA) <- c("Sigma11", "Sigma12", "Sigma12", "Sigma22") 
   
  # sampler
  for (i in 1:nsim) {
    # UPDATE alpha beta and sigma2
    alpha <- UpdateAlpha(Y = y, X = x1, alpha = alpha, b1 = b[, 1], V = Valpha,
                         m0 = priors$malpha, V0inv = prior.Vinvalpha) 
    beta <- UpdateBeta(Y = y, logY = ly, X = x2, b2 = b[, 2], sigma2inv = sigma2inv,
                       m0 = priors$mbeta, V0inv = prior.Vinvbeta)
    sigma2inv <- UpdateSigma2inv(Y = y, X = x2, logY = ly, prior.shape = priors$sigma2.shape, 
                                 prior.rate = priors$sigma2.rate, b = b, beta = beta)
    
    # UPDATE b_i (with Metropolis) and Sigma 
    b.results <- Updateb(Y = y, X1 = x1, X2 = x2, b = b, ni = ni, id = id,
                         Sigma.prop = 1.3 * Sigma, Sigma = Sigma,
                 alpha = alpha, beta = beta, sigma2 = 1/sigma2inv)
    b <- b.results$b
    Sigma <- UpdateSigma(prior.v = priors$Sigma.v, prior.S = priors$Sigma.S, b = b)
    
    # store parameters
    print(i)
    if (i > burn & i %% thin==0) {
      j <- (i-burn)/thin
      Alpha[j, ] <- alpha
      Beta[j, ] <- beta
      Sigma2[j] <- 1/sigma2inv
      B1[j, ] <- b[, 1]
      B2[j, ] <- b[, 2]
      SIGMA[j, ] <- c(Sigma)
      b.acc <- b.acc + b.results$acc
      B.acc[j, ] <- b.acc/j
    }
    #if (i %% 100==0) print(i)
  }
  return(list(alpha = Alpha, beta = Beta, sigma2 = Sigma2, 
              b1 = B1, b2 = B2, Sigma = SIGMA, b.acc = B.acc))
}
