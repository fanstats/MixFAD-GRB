# CEM algorithm for Gaussian Mixture of Factor Analyzers

# GMM: Mixtures of K Gaussian distributions
# Dimension reduction: For each mixture: q.k-factor model 
#                      -- assume common # of factors: q.1 <-,..,<- q.K <- q
# The (k,q)-factor model is selected by BIC

# Inputs:
# X: n by p data matrix
# init: a list of initial parameters M, L, D, prob (optional)

source("initial.R")
source("cem_cycle.R")
source("logf.R")
source("utilis.R")

mixfad <- compiler::cmpfun(function(X, K, q, init = NULL,
                                    maxiter = 500, epsi = 1e-4){
  p <- ncol(X) # number of coordinates
  n <- nrow(X) # number of obs.
  iter <- 0    # iteration
  
  # Initialize K Mean & Covariance Components
  if (is.null(init)) {
    init <- random.init(p,q,K)
  }
  M <- init$M
  L <- init$L
  D <- init$D
  prob <- init$prob
  
  pf <- logf(X,M,L,D,prob) + rep(log(prob), each = n)
  llk0 <- sum(log(rowSums(exp(pf))))
  
  for (iter in 0:maxiter){
  
    out <- cem_cycle(X,M,L,D,prob)
    M <- out$M
    L <- out$L
    D <- out$D
    prob <- out$prob
    llk <- out$loglik
    iter <- iter+1
    
    if (abs((llk-llk0)/llk) < epsi){
      break
    }
    
    llk0 <- llk
  }
  
  return(list(iter = iter, 
              mmu = out$M, mD = out$D, mL = out$L, prob = out$prob,
              loglik = out$loglik, 
              bic = -2*out$loglik + log(n)*K*(2*p + p*q + 1)
  ))
  
})
    