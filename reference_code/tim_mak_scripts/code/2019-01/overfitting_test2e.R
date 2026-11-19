Tmisc()
library(mvtnorm)
setMKLthreads(8)
filename <- Rfilename("overfitting_test2d", seed=1000)
n <- 200
m <- 2000
overlap <- 0
p <- 200
# a <- 20
pi <- overlap/max(n,m)
b <- 1


# a2 <- a^2
# M <- rMat(a2, n+m-overlap)
# corM <- cor(M)

X <- rMat(n, p)
Z <- rMat(m, p)

XZT <- tcrossprod(X, Z)

Tr1 <- sum(XZT^2)
sqrt(Tr1) / sqrt(p^2*n)

cpZ <- crossprod(Z)
cpX <- crossprod(X)

W <- cpZ %*% cpX
Tr2 <- sum(cpZ * W)
Tr2b <- sum(diag(cpZ %*% cpX %*% cpZ))
print(Tr2/(n*p^2*m))

# For m < 5000 or so only
V <- crossprod(XZT)
tcpZ <- tcrossprod(Z)
Tr2c <- sum(V*tcpZ)
sum(diag(V*tcpZ))/Tr2