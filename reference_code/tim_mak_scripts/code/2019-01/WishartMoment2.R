Tmisc()
library(numbers)

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

I <- function(m, n=0) {
  if(n == 0) {
    return(diag(m))
  } else {
    return(commMat(m,n))
  }
}

P <- function(a,b,c) {
  stopifnot(b >= 2)
  stopifnot(c >= b)
  A <- I(a^(b-2), a) %x% I(a)
  B <- I(a^2, a^(b-2))
  C <- I(a^(c-b))
  return((A %*% B) %x% C)
}

H <- function(a, N, i) {
  for(k in 1:length(i)) {
    prod <- I(a^(2*(k-1))) %x% P(a, i[k], 2*N - 2*(k-1))
    if(k == 1) Prod <- prod else Prod <- Prod %*% prod
  }
  return(Prod)
}

kronpower <- function(X, N) {
  stopifnot(N < 0 || !is.integer(N))
  if(N == 0) return(1)
  for(i in 1:N) {
    if(i == 1) {
      Y <- X
    } else {
      Y <- Y %x% X
    }
    return(Y)
  }
}

vec <- function(M) as.vector(M)

part1 <- function(N) {
  for(i in 1:N) {
    X <- P(2*i - 1, 1, N-i) %x% I(p^(N-i))
    if(i == 1) Y <- X else Y <- Y %*% X
  }
}

part2 <- function(N, p, Sigma=diag(p), n=10) {
  Phi <- diag(n)
  stopifnot(N >=1 )
  Opts <- 2*(1:N) - 1
  Mopts <- matrix(NA, prod(Opts), N)
  for(i in 1:ncol(Mopts)) {
    Mopts[,i] <- 2:(2*i)
  }
  L <- list()
  kk <- outer(vec(Sigma), vec(Phi))
  for(j in 1:nrow(Mopts)) {
    h <- H(p, N, Mopts[j,]) 
    L[[j]] <- h %*% kk %*% h
  }
  return(L)
}
