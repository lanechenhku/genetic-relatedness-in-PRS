Tmisc()
setMKLthreads(8)
Tim.load(Rplink)
Tim.load(lassosum)
filename <- Rfilename("overfitting_test2b", seed=1000)
testn <- 700
validn <- 700
p <- 4000

n.causal <- p
herit <- 0.8
reps <- 40
Trainn <- c(700, 1400)

for(j in 1:length(Trainn)) {
  
  trainn <- Trainn[j]
  n <- testn + trainn + validn
  
  group <- sample(c(rep(0,trainn), rep(1, testn), rep(2,validn)))
  trainsample <- group == 0
  testsample <- group == 1
  validsample <- group == 2
  
  X <- rMat(n,p)
  X <- scale(X)
  
  for(i in 1:reps) {
    
    # i <- 1
    beta <- rep(0, p)
    beta1 <- scale(rnorm(n.causal))
    beta[sample(p, n.causal)] <- beta1
    Xb <- as.vector(X %*% beta)
    error <- scale(rnorm(n, sd = sqrt(var(Xb)/herit*(1-herit))), scale = F)
    pheno <- Xb + error
    # sd.pheno <- sd(pheno)
    sd.pheno <- 1
    pheno <- pheno / sd.pheno
    beta <- beta /sd.pheno
    Xb <- Xb / sd.pheno
    error <- error / sd.pheno
    
    # bhat <- as.vector(cov(pheno[trainsample], X[trainsample,])) * (trainn - 1)
    # bhat1 <- as.vector(cov(Xb[trainsample], X[trainsample,])) * (trainn - 1)
    # bhat2 <- as.vector(cov(error[trainsample], X[trainsample,])) * (trainn - 1)
    bhat <- as.vector(cov(pheno[trainsample], X[trainsample,])) 
    bhat1 <- as.vector(cov(Xb[trainsample], X[trainsample,])) 
    bhat2 <- as.vector(cov(error[trainsample], X[trainsample,])) 
    Xbhat <- as.vector(X[testsample, ] %*% bhat)
    Xbhat1 <- as.vector(X[testsample, ] %*% bhat1)
    Xbhat2 <- as.vector(X[testsample, ] %*% bhat2)
    Xb <- Xb[testsample]
    error <- error[testsample]
    
    cov_Xbhat_Xb <- as.vector(cov(Xb, Xbhat)) * (testn - 1)
    cov_Xbhat1_Xb <- as.vector(cov(Xb, Xbhat1)) * (testn - 1)
    cov_Xbhat2_Xb <- as.vector(cov(Xb, Xbhat2)) * (testn - 1)
    var_Xbhat <- as.vector(var(Xbhat)) * (testn - 1)
    var_Xbhat1 <- as.vector(var(Xbhat1)) * (testn - 1)
    var_Xbhat2 <- as.vector(var(Xbhat2)) * (testn - 1)
    var_Xbhat12 <- as.vector(cov(Xbhat1, Xbhat2)) * (testn - 1)
    var_Xb <- as.vector(var(Xb)) * (testn - 1)
    
    cor <- cor(Xbhat, Xb)
    
    print.counter(i)
    results.table.vars <- c("cov_Xbhat_Xb", "cov_Xbhat1_Xb", "cov_Xbhat2_Xb", 
                            "var_Xbhat", "var_Xbhat1", "var_Xbhat2", 
                            "var_Xbhat12", "var_Xb", "cor", "var_b", "var_e", 
                            "trainn")
    var_e <- var(error)
    var_b <- var(beta)
    fill.in.results.table()
    
  }
}

d <- subset(results.table, trainn==700)
d2 <- subset(results.table, trainn==1400)
d$cor2 <- with(d, (cov_Xbhat1_Xb + cov_Xbhat2_Xb) / 
                 sqrt(var_Xbhat1 + 2*var_Xbhat12 + var_Xbhat2) /
                 sqrt(var_Xb))
summary(d$cor - d$cor2) # test

# Checks
summary(d$cov_Xbhat1_Xb)
summary(d$cov_Xbhat2_Xb)
summary(d$var_Xbhat1)
summary(d$var_Xbhat2)
summary(d$var_Xbhat12)
summary(d$var_Xb)
summary(d$cor)

Z <- X[trainsample,]
W <- X[testsample,]
cp <- tcrossprod(W, Z)
Tr <- sum(cp^2) 
ZWT <- tcrossprod(Z, W)
ZZT <- tcrossprod(Z)
Tr2 <- sum(tcrossprod(ZWT) * ZZT)

# Checks2
summary(d2$cov_Xbhat1_Xb)
summary(d$cov_Xbhat1_Xb)
summary(d2$var_Xbhat1)
summary(d$var_Xbhat1)



