random.init = function(p,q,K){
  # Purpose: randomly initialize M,L,D and set pk = 1/K
  M = D = SD <- matrix(NA,p,K)
  L <- array(NA,c(p,q,K))
  
  for (k in 1:K) {
    M[,k] <- rnorm(p)
    Lk <- matrix(rnorm(p*q),p,q)
    D[,k] <- runif(p,0.2,0.8)
    Gamma <- crossprod(Lk,Lk/D[,k])
    eig <- eigen(Gamma)$vectors
    L[,,k] <- Lk%*%eig
  }
  
  prob = rep(1/K,K)
  
  return(list(M = M,L = L,D = D,prob = prob))
}

random.init.error = function(p,q,K){
  # Purpose: randomly initialize M,L,D and set pk = 1/K
  M = U <- matrix(NA,p,K)
  L <- array(NA,c(p,q,K))
  
  tau <- runif(p, 0.002,0.008)
  alpha <- rep(15, p)
  
  for (k in 1:K) {
    M[,k] <- rnorm(p)
    Lk <- matrix(rnorm(p*q),p,q)
    U[,k] <- runif(p,0.2,0.8)
    #D[,k] <- U[,k]+tau
    Gamma <- crossprod(Lk,Lk/U[,k])
    eig <- eigen(Gamma)$vectors
    L[,,k] <- Lk%*%eig
  }
  
  prob = rep(1/K,K)
  
  return(list(M = M,L = L,U = U,tau = tau,alpha = alpha,prob = prob))
}
