Tmisc()
library(mvtnorm)
setMKLthreads(8)
filename <- Rfilename("overfitting_test2d", seed=1000)
n <- 1000
m <- 2000
overlap <- 0
p <- 20000
a <- 10
pi <- overlap/max(n,m)
b <- 1


a2 <- a^2

for(j in 1:50) {
  M <- rMat(a2, n+m-overlap)
  corM <- cor(M)
  
  nsamples <- 50
  j <- sample(m, nsamples)
  k <- sample(m, nsamples)
  discard <- j == k
  j <- j[!discard]
  k <- k[!discard]
  nsamples <- length(j)
  r <- sample(n, nsamples)
  rho_jk <- corM[n-overlap+j, n-overlap+k]
  rho_jr <- corM[n-overlap+j, r]
  rho_kr <- corM[n-overlap+k, r]
  M <- array(M,dim = c(nsamples, nsamples, nsamples))
  for(i in 1:nsamples) {
    M[,,i] <- outer(rho_jr[,i], rho_kr[,i]) * rho_jk
  }
  fill.in.results.table(mean3=mean(as.vector(M)), 
                        mean1=mean(offdiags(corM)), 
                        mean2=mean(as.vector(rho_jk^2)))
}

