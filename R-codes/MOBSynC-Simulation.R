source("logf.R")
source("utilis.R")
source("initial.R")
library(mclust)
library(MixSim)
library(RColorBrewer)

# Load the 3D shape dataset
load("mobsync-shapedata.rda")
X.sim = as.matrix(df[,-4])
id = df[,4]

# Fit with GMM
bics = mclustBIC(X.sim,modelName = "VVV")
which.max(bics) # optimal 5 groups
res.sim <- Mclust(X.sim, modelName = "VVV", G = 5)

# Initial clustering results

# estimated clusters
cl = res.sim$classification
# match the estimated cluster labels and the true labels
id = EMCluster::recolor(cl, as.numeric(id)) 
# confusion matrix
table(id,cl)
# classification error
classError(cl, id)

paras = res.sim$parameters
ppi = paras$pro
mmu = paras$mean
s = paras$variance$sigma
K = 5
# generalized overlap
overlapGOM(Pi = ppi, Mu = t(mmu), S = s)

# map of pairwise overlaps 
colpal.5 = c(brewer.pal(8,"Dark2")[3:5], brewer.pal(8,"Set2")[1:2])

overlap_mat <- overlap(ppi, t(mmu), s)$OmegaMap
overlap.map((2*overlap_mat - diag(K)), lab.col = colpal.5, 
            map.cex = 2.75, lab.cex = 4,
            legend.cex = 1.5, font = 2, scale.pos = 0.35, legend.width = 0.1)

# Merge initial clusters 1,2,3 and clusters 4,5 
# kappa = 1

# relabel for composite clusters
cl.m = cl
ind = which(cl.m==1)
cl.m[ind] = 1
ind = which(cl.m==2)
cl.m[ind] = 1
ind = which(cl.m==3)
cl.m[ind] = 1
ind = which(cl.m==4)
cl.m[ind] = 2
ind = which(cl.m==5)
cl.m[ind] = 2
cl.m = factor(cl.m)
table(cl.m)

## Compute overlaps of composite clusters by Monte Carlo (MC) method
size = 1e6
K = nlevels(cl.m)
W = diag(K)
prob = paras$pro
p = ncol(X.sim)

for (i in 1:K) {
  for (j in 1:K) {
    # omega_j|i: sample Xi from group i 
    # check if (log density for group j - log density group i) > 0
    
    # density for i
    if(i == 1){
      ppi = prob[1:3]/sum(prob[1:3])
      mm = mmu[,1:3]
      ss = s[,,1:3]
      
      A = MixSim::simdataset(n = size, Pi = ppi, Mu = t(mm), S = ss)
      xi = A$X
      
      f <- exp(logf.nfa(xi,mm,ss,prob[1:3]) + rep(log(prob[1:3]), each = nrow(xi)))
      lf.i <- log(rowSums(f))
      
    }
    
    if(i == 2){
      ppi = prob[4:5]/sum(prob[4:5])
      mm = mmu[,4:5]
      ss = s[,,4:5]
      
      A = MixSim::simdataset(n = size, Pi = ppi, Mu = t(mm), S = ss)
      xi = A$X
      
      f <- exp(logf.nfa(xi,mm,ss,prob[4:5]) + rep(log(prob[4:5]), each = nrow(xi)))
      lf.i <- log(rowSums(f))
      
    }
    
    # density for j
    if(j == 1){
      ppi = prob[1:3]/sum(prob[1:3])
      mm = mmu[,1:3]
      ss = s[,,1:3]
      
      A = MixSim::simdataset(n = size, Pi = ppi, Mu = t(mm), S = ss)
      xi = A$X
      
      f <- exp(logf.nfa(xi,mm,ss,prob[1:3]) + rep(log(prob[1:3]), each = nrow(xi)))
      lf.j <- log(rowSums(f))
    }
    
    if(j == 2){
      ppi = prob[4:5]/sum(prob[4:5])
      mm = mmu[,4:5]
      ss = s[,,4:5]
      
      A = MixSim::simdataset(n = size, Pi = ppi, Mu = t(mm), S = ss)
      xi = A$X
      
      f <- exp(logf.nfa(xi,mm,ss,prob[4:5]) + rep(log(prob[4:5]), each = nrow(xi)))
      lf.j <- log(rowSums(f))
    }
    
    W[j,i] = sum(lf.j>lf.i)/size
  }
}

# matrix of pairwise overlaps
O = W + t(W)
diag(O) = rep(1,K)
O

# generalized overlap
gomega = (max(eigen(O)$values)-1)/(K-1)

# map of pairwise overlaps
overlap.map(O, lab.col = colpal[1:4], 
            map.cex = 2.75, lab.cex = 4,
            legend.cex = 1.5, font = 2, scale.pos = 0.35, legend.width = 0.1)

## gomega < 0.001, stop merging.

