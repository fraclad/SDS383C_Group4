
# Finds the L2 norm of point x in R2
vnorm <- function(x) {
  return(sqrt(sum(x^2)))
}
# Evaluates log of vMF constant
log_vMF_c <- function(tau, d) {
  return((d/2-1)*log(tau)-log(besselI(tau,nu=(d/2-1),expon.scaled = TRUE)))
}

# Evaluate Fisher-Gaussian pdf at certain point
# Inputs: 
  # - x: point in Rd
  # - parameters:
  #     - ctr (center)
  #     - r (radius)
  #     - sigma.sq (variance for normal)
  #     - mu
  #     - tau a scalar (specifies where on sphere points are concentrated)
dfg <- function(x, ctr, r, sigma.sq, mu, tau) {
  d = length(x)
  z = x - ctr
  denom = vnorm(tau*mu + r*z/sigma.sq)
  term <- log_vMF_c(tau, d)-(d/2)*log(2*pi*sigma.sq)-log_vMF_c(denom,d) - (1/(2*sigma.sq))*(vnorm(z)^2+r^2-2*sigma.sq*denom)-tau
  return(exp(term))
  }

# Testing to see if dfg gives me a number
# x <- c(1,1)
# center <- c(0, 0)
# r <- 1
# variance <- 0.001
# mu <- c(1/sqrt(2), 1/sqrt(2))
# tau <- 3
# dfg(x, center, r, variance, mu, tau)

# Given results of mixture model, evaluate density at test points
# Stuff inside results: "Sp_wgts","Kern_wgts","Centers","Radius","mu","tau","trace_sigma","likelihood_chain","inclusion_matrix"
mixture_density <- function(X, results) {
  # Get parameters
  ctr <- results$Centers
  r <- results$Radius
  mu <- results$mu
  tau <- results$tau
  lambda <- results$Sp_wgts
  pi <- results$Kern_wgts
  f_xs <- NULL
  K <- nrow(tau)
  M <- ncol(tau)
  D <- ncol(X)
  ss <- results$trace_sigma[length(results$trace_sigma)]
  
  # Iterate through points in X matrix and evaluate density
  for (i in 1:nrow(X)) {
    x0 <- X[i,] # grab row, this is one x
    fvmf <- NULL
    for (l in 1:K) { # Each vMF kernel
      fvmf[l] <- sum(sapply(1:M, function(m){return(pi[l,m]*dfg(x0,ctr[,l],r[l],ss,mu[,l,m],tau[l,m]))}))
    }
    f_xs[i] <- t(fvmf)%*%lambda
  }
  return(f_xs)
}

# NOT FINISHED! I AM NOT SURE HOW TO DO THIS
# Returns 2 things: pts (predictive without errors), w (with normal errors)
plot_mixture_density <- function(N, results) {
  # Get parameters
  ctr <- results$Centers
  r <- results$Radius
  mu <- results$mu
  tau <- results$tau
  lambda <- results$Sp_wgts
  pi <- results$Kern_wgts
  f_xs <- NULL
  K <- nrow(tau)
  M <- ncol(tau)
  D <- ncol(X)
  ss <- results$trace_sigma[length(results$trace_sigma)]
  
  # Initialize matrix to hold the points
  pts <- matrix(data=NA, nrow=N, ncol=D)
  w <- matrix(data=NA, nrow=N, ncol=D)
  
  
  
}







