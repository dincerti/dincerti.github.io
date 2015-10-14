
# EVALUATE CLASSIFIER ----------------------------------------------------------
Evaluate <- function(prob, y){
  # Uses ggplot to plot an ROC curve.
  #
  # Args:
  #   prob: matrix with each column equal to the predicted probability that 
  #         y = 1 for model j.
  #   y: test data with observed 0 or 1 values.
  #
  # Returns:
  #   ggplot object
  w <- which(y == 1)
  cutoff <- seq(0, .99,length = 1000) 
  perf <- matrix(NA, nrow = length(cutoff), 4) 
  colnames(perf) <- c("cutoff", "tpr", "tnr", "acc")
  for (i in 1:length(cutoff)){
    yhat <- 1 * (prob >= cutoff[i])
    perf[i, 1] <- cutoff[i]
    perf[i, 2] <- mean(yhat[w] == 1) 
    perf[i, 3] <- mean(yhat[-w] == 0) 
    perf[i, 4] <- mean(y == yhat) 
  }
  return(perf)
}

# AREA UNDER ROC CURVE ---------------------------------------------------------
AUC <- function(x, y){
  h <- (y[-1] + y[-length(y)])/2
  w <- -diff(x)
  return(sum(h * w))
}