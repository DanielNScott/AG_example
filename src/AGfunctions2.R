## Functions to fit ancestral graphs from fmri data
require(ggm)

# make object "brain"
as.brain <- function(acts, rois = NA, cond = NA, ...) {
  if(!is.na(rois[1])) acts <- acts[,rois]
  if(!is.na(cond[1])) acts <-acts[cond,]
  attr(acts,"n") <- length(cond)
  class(acts) <- c("brain","data.frame")
  return(acts)
}


getAIC <- function(fit, S, n) {
  p <- length(fit$Shat[1,])
  q <- p*(p - 1) /2 - fit$df
  npars <- p + q
  fitS  <- -(n/2)*log(2*pi) -(n/2)*log(det(S)) -(n*p)/2
  aic   <- fit$dev + fitS + 2*npars
  return(list(aic = aic, npars = npars))
}


fitAG <- function(acts, adjacency, rois = NA, cond = NA, detail = "AIC", ...) {
  # Fits all subjects' data using an AG model
  #
  # Args:
  #   acts     : Array-like activation data. Trials as rows and rois as columns.
  #   adjacency: Adjacency matrix
  #   rois     : Regions of interest
  #   cond     : conditions
  #   detail   : Level of detail required: 'AIC', 'LR', or 'both'
  #
  # Returns:
  #   

  N <- length(acts)
  n <- length(cond)
  p <- length(rois)
  results <- numeric(0)
  
  for(i in 1:N) {

    data  <- as.brain(acts[[i]], rois, cond)
    covar <- makeCov(data)
    
    covar.fit <- fitAncestralGraph(adjacency, covar, attr(covar,"n"))

    chi2  <- T2(s.fit$dev, n)
    p     <- pchisq(chi2, s.fit$df, lower.tail = FALSE)    
    fit   <- list(chi2, p)

    aic   <- getAIC(s.fit, s, n)
  
    results <- c(results, fit, aic$aic, aic$npars)

  }
  return(results)
}

makeCov <- function(acts) {
    S <- cov(acts)
    attr(S,"n") <- attr(acts,"n")
    return(S)
}

T2 <- function(fit, n) {
    fit/(1+(fit/n))
}
