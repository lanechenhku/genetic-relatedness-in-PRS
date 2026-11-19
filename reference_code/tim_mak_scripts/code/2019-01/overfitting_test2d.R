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
M <- rMat(a2, n+m-overlap)
corM <- cor(M)

X_Z <- t(rmvnorm(p, sigma=corM))
X <- X_Z[1:n, ]
Z <- X_Z[(n-overlap+1):(m+n-overlap), ]

XZT <- tcrossprod(X, Z)

Tr1 <- sum(XZT^2)
print(Tr1/(p^2*n*m))
components1 <- c(pi/max(n,m), 1/a2, 1/p)
print(c(components1, sum(components1)))

W <- crossprod(XZT, XZT) 
tcpZ <- tcrossprod(Z)
Tr2 <- sum(W * tcpZ)
print(Tr2/(n*p^2*m*(3*m+p)))

Tr1 / sqrt(Tr2) / sqrt(p*n)

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
hist(as.vector(M), 200)
  