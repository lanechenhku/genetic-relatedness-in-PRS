Tmisc()
setMKLthreads(8)
Tim.load(Rplink)
Tim.load(lassosum)
filename <- Rfilename("overfitting_test3", seed=1000)
testn <- 700
trainn <- 700
validn <- 700
p <- 100
n <- testn + trainn + validn

group <- sample(c(rep(0,trainn), rep(1, testn), rep(2,validn)))
trainsample <- group == 0
testsample <- group == 1
validsample <- group == 2

n.causal <- 100
herit <- 0.8
reps <- 100

for(i in 1:reps) {
  # i <- 1
  beta <- rep(0, p)
  beta[sample(p, n.causal)] <- rnorm(n.causal)
  X <- rMat(n,p)
  Xb <- as.vector(X %*% beta)
  error <- rnorm(n, sd = var(Xb)/herit*(1-herit))
  pheno <- Xb + error
  
  cor <- coef(lm(pheno[trainsample] ~ . -1, data=as.data.frame(X[trainsample,])))
  # cor <- as.vector(cor(pheno[trainsample], X[trainsample,]))
  est <- as.vector(X %*% cor)
  
  #### Tests ####
  # lapply(0:2, function(i) cor(est[group == i], Xb[group == i]))
  # lapply(0:2, function(i) cor(est[group == i], error[group == i]))
  #### Cor with X ####
  train.cor <- as.vector(cor(est[trainsample], X[trainsample,]))
  # train.cor0 <- as.vector(cor(Xb[trainsample], X[trainsample,]))
  test.cor <- as.vector(cor(est[testsample], X[testsample,]))
  # test.cor0 <- as.vector(cor(Xb[testsample], X[testsample,]))
  print.counter(i)
  fill.in.results.table(i = i, 
                        train.cor=cor(est[trainsample], Xb[trainsample]), 
                        test.cor=cor(est[testsample], Xb[testsample]), 
                        train.cov=cov(est[trainsample], Xb[trainsample]), 
                        test.cov=cov(est[testsample], Xb[testsample]), 
                        train.sd.est=sd(est[trainsample]), 
                        test.sd.est=sd(est[testsample]), 
                        train.sd.Xb=sd(Xb[trainsample]), 
                        test.sd.Xb=sd(Xb[testsample]), 
                        sd.train.cor = sd(train.cor), 
                        sd.test.cor = sd(test.cor))

}

boxplot(results.table[,c("sd.train.cor", "sd.test.cor")])
boxplot(results.table[,c("train.cor", "test.cor")])
boxplot(results.table[,c("train.cov", "test.cov")], 
        ylim=c(0, max(results.table[,c("train.cov", "test.cov")])))
boxplot(results.table[,c("train.sd.est", "test.sd.est")], 
        ylim=c(0, max(results.table[,c("train.sd.est", "test.sd.est")])))
boxplot(results.table[,c("train.sd.Xb", "test.sd.Xb")])
