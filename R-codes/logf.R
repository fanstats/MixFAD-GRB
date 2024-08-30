logf <- compiler::cmpfun(function(X,M,L,D,prob){
  # Purpose: compute the n x K matrix of logfk(Xi)
  # X: n by p data matrix
  n <- nrow(X)
  p <- ncol(X)
  K <- length(prob)
  
  logf <- matrix(NA,n,K)
  
  for (k in 1:K) {
    muk <- as.numeric(M[,k])
    Dk <- as.numeric(D[,k])
    Lk <- as.matrix(L[,,k])
    
    cmatdg <- fad:::cmdg(Lk,Dk)
    logdet <- sum(log(Dk)) + 2*sum(log(1/sqrt(cmatdg)))
    xm <- t(X) - muk
    comp <- colSums(xm^2/Dk) - colSums(Lk%*%(cmatdg*crossprod(Lk,xm/Dk))*xm/Dk)
    logf[,k] <- {-p}*log(2*pi)*0.5 - logdet*0.5 - comp*0.5 
    
  }
  
  return(logf)
})