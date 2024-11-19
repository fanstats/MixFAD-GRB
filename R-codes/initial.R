random.init = function(p,q,K){
  # Purpose: randomly initialize M,L,D and set pk = 1/K
  if(q>0){
    M = D <- matrix(NA,p,K)
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
    
  }else{
    M <- matrix(NA,p,K)
    S <- array(NA,c(p,p,K))
    for (k in 1:K) {
      M[,k] <- rnorm(p)
      A <- matrix(rnorm(p^2), p, p) 
      S[,,k] <- t(A) %*% A
    }
    prob = rep(1/K,K)
    
    return(list(M = M, S = S, prob = prob))
  }
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

# data-driven initial
data.init = function(X, prop, p, q, K,seed){
  options(warn = -1)
  
  n = nrow(X)
  size = round(prop*n)
  
  set.seed(seed)
  x = X[sample.int(n,size = size,replace = F),]
  cl.kmeans = kmeans(x,centers = K,iter.max = 1,nstart = 1)$cluster
  
  M = D = matrix(NA,p,K)
  L = array(NA,c(p,q,K))
  for (k in 1:K) {
    x.k = x[cl.kmeans==k,]
    res.k = fad::fad(x.k,q)
    M[,k] = colMeans(x.k)
    D[,k] = res.k$sd^2*res.k$uniquenesses
    L[,,k] = res.k$sd*res.k$loadings
  }
  prob = rep(1/K,K)
  
  return(list(M = M,L = L,D = D,prob = prob))
}

