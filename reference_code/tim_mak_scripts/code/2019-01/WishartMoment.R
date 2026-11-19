Tmisc()
library(numbers)
Sigma <- c(2,3,5,
           3,7,11,
           5,11,13)
Sigma <- matrix(Sigma, 3,3)
index <- c(18,2)

W <- Sigma %x% Sigma %x% Sigma
vecSigma <- as.vector(Sigma)
V <- outer(vecSigma, vecSigma) %x% Sigma

commMat <- function(m,n) {
  MMM <- matrix(integer(0), 0, m*n)
  
  for(i in 1:n) {
    MM <- matrix(integer(0), m, 0)
    for(j in 1:m) {
      M <- matrix(0, m, n)
      M[j,i] <- 1
      MM <- cbind(MM, M)
    }
    MMM <- rbind(MMM, MM)
  }
  return(MMM)
}
P <- diag(3) %x% commMat(3,3)
Q <- commMat(3,3) %x% diag(3)

p <- function(n) paste(primeFactors(n), collapse = ",")
pp <- function(M) return(p(M[index[1], index[2]]))
ppp <- function(x) {
  MM <- NULL
  for(i in nchar(x):1) {
    mat <- substr(x, i, i)
    M <- get(mat, envir=.GlobalEnv,inherits = F)
    if(is.null(MM)) MM <- M else {
      MM <- M %*% MM
    }
  }
  return(pp(MM))
}
#### Terms ####
n3terms <- c("W")
n2terms <- c("V", "QW", "QPVPQ", "QPW", "PW", "PVP", "QPQW")
n1terms <- c("QPVP", "PQW", "QPV", "VP", "VPQ", "PV", "PVPQ")

print(Sigma)
sapply(n3terms, ppp)
sapply(n2terms, ppp)
sapply(n1terms, ppp)
