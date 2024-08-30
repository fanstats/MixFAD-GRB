cem_cycle <- compiler::cmpfun(function(X,M,L,D,prob){
  
  # Purpose: one CEM cycle for ML method
  # X: n by p data matrix
  n <- nrow(X)
  p <- ncol(X)
  q <- dim(L)[2]
  K <- length(prob)
  
  ### E-step: obtain n x K matrix of gamma_{i,k}: P(Z_i = k|X_i) ###
  Ga <- matrix(NA,n,K)
  pf <- logf(X,M,L,D,prob) + rep(log(prob), each = n)
  for (k in 1:K) {
    Ga[,k] <- 1/rowSums(exp(pf - pf[,k]))
  }
  
  pf <- matrix(NA,n,K)
  for (k in 1:K) {
    
    ### CM-step: prob_k ###
    nk <- sum(Ga[,k])
    pk <- nk/n
    
    ### CM-step: mu_k ###
    muk <- colSums(Ga[,k]*X)/nk
    
    ### CM-step: D_k, L_k ###
    X.c <- sweep(X, 2, muk)
    X.wc <- sqrt(Ga[,k]/pk)*X.c
    sd <- sqrt(colMeans(X.wc^2))
    
    fad.k <- try(fad:::fad.fit.X(X = X.wc, q = q, 
                           iSD = 1/sd/sqrt(n), mu = rep(0,p)),
                 silent = T)
    if(class(fad.k) == "try-error"){ 
      break
    }
    
    Lk <- sd*fad.k$loadings
    Dk <- sd^2*fad.k$uniquenesses
    
    prob[k] <- pk
    M[,k] <- muk
    L[,,k] <- Lk
    D[,k] <- Dk
    
  }
  
  prob <- prob/sum(prob)
  
  f <- exp(logf(X,M,L,D,prob) + rep(log(prob), each = n))
  mix.pf <- log(rowSums(f*is.finite(f), na.rm=TRUE))
  loglik <- sum(mix.pf*is.finite(mix.pf), na.rm=TRUE)

  
  return(list(prob = prob, M = M, L = L, D = D, loglik = loglik))
})