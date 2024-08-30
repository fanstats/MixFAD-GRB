# A CEM algorithm for Mixture of Factor Analyzers

# Latent variables: Mixtures of K Gaussian
# Dimension reduction: For each Gaussian, assume a q_k-factor model
#                      with common # of factors: q_1 = q_2 =,..,= q_K = q
# Model selection: The (k,q)-factor model is selected by BIC

# Inputs:
# X: A n by p data matrix
# init: A list of initial parameters M, L, D, prob (optional)

# Application: 
# init = random.init(p = p, q = q, K = k)
# res = mixfad(X = data_matrix, K = k, q = q, init = init)

source("initial.R")
source("cem_cycle.R")
source("logf.R")
source("utilis.R")

mixfad <- compiler::cmpfun(function(X, K, q, init = NULL,
                                    maxiter = 500, epsi = 1e-4){
  X = as.matrix(X)
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
  
  f <- exp(logf(X,M,L,D,prob) + rep(log(prob), each = n))
  mix.pf <- log(rowSums(f*is.finite(f), na.rm=TRUE))
  llk0 <- sum(mix.pf*is.finite(mix.pf), na.rm=TRUE)
  
  for (iter in 0:maxiter){
  
    out <- cem_cycle(X,M,L,D,prob)
    M <- out$M
    L <- out$L
    D <- out$D
    llk <- out$loglik
    iter <- iter+1
    
    if(abs(llk) < 10^-6){
      llk.dif = epsi
    }else{
      llk.dif = abs((llk-llk0)/llk)
      }
    
    if (llk.dif < epsi){
      break
    }
    
    llk0 <- llk
  }
  
  return(list(iter = iter, init = init,
              mmu = out$M, mD = out$D, mL = out$L, prob = out$prob,
              loglik = out$loglik, 
              bic = -2*out$loglik + log(n)*K*(2*p + p*q + 1)
              #ebic = -2*out$loglik + (log(n) + 2*max(0,1-0.5/log(p,base = n))*log(p))*k*(2*p+p*q+1))
  ))
  
})
    