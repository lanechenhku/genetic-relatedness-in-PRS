# snp wise correlation
#packages used
if(!("devtools" %in% rownames(installed.packages()))) install.packages("devtools"); .rootdir <- "/home/unfated/Tim"; source("/home/unfated/Tim/Tmisc/inst/Startup.R")
#17
## fixed.mem = 40gb
# a package
Tmisc()
Tim.load("Rplink",lib.dir="/home/unfated/Tim")
filename <- Rfilename("sim.10000.1000", seed = "prompt")
1
#0
# test.n <- 2000
# validation.n <- 2000
# m <- as.integer(nameoption(filename, 2))
# n.causal <- as.integer(nameoption(filename, 3))
library("RcppArmadillo")
library("data.table")
library("Matrix")

Tim.load("UKBBscripts", "/psychipc01/disk3/data/UKBB")
Tim.load("lassosum","/home/unfated/Tim")
Tim.load("crosspred", "/home/unfated/Tim")
#cl <- startParallel(10)

fixed.mem <- 40000
#bfile <- attachroot("~/DATA/UKBB2/data/genotype/whiteBritish1_train/whiteBritish1_train")
#bfile <- attachroot("/psychipc01/disk3/data/UKBB/data/genotype/whiteBritish1_train/whiteBritish1_train")
# bfile<-'/psychipc01/disk3/data/UKBB/data/genotype/imputed_0.01/maf_0.01_1'
bfile<-"/psychipc01/disk3/data/UKBB/data/genotype/whiteBritish1_train/whiteBritish1_train1"
#bfile<-'/home/unfated/chr1_block/block_1'
#bfile <- useWTCCC(.bfile = F)

N <- nrow.bfile(bfile)
m <- 100000
set.seed(1)
discovery.sample <- sample(N, m)
# others.sample <- setdiff(1:N, discovery.sample)
# validation.sample <- sample(others.sample, validation.n)
# others.sample <- setdiff(1:N, c(discovery.sample, validation.sample))
# 
discovery <- logical.vector(discovery.sample, N)
# validation <- logical.vector(validation.sample, N)

p <- ncol.bfile(bfile)
prop = 0.01
n.causal = ceiling(p*prop)
# beta <- rep(0, p)
library(STAT)
mu=rep(0,n.causal)
set.seed(1)
V=rWishart(1,n.causal,diag(n.causal))
V2=rWishart(1,n.causal,V[,,1])
corr=cov2cor(as.matrix(V2[,,1]))
library(MASS)
#cov=matrix(c(1,.5,.15,.15,0,.5,1,0.15,0.15,0,.15,.15,1,.25,0,.15,.15,.25,1,0,0,0,0,0,1),nrow=5,ncol=5)
set.seed(1)
effect=mvrnorm(n = 100, mu=mu, Sigma=corr, tol = 1e-6, empirical = FALSE)
Zscore=list()
Zscorefull=list()
for(i in 1:nrow(effect)){
beta <- rep(0, p)
d=n.causal+1
beta[c(2:d)]<- effect[i,]
# write.table(beta,'effectsize_chr1.txt')
#Xb <- pgs(bfile=bfile, weights=beta, cluster=cl, trace=1, keep=discovery)
# Xb <- pgs(bfile=bfile, weights=beta, chr=c(1,2), trace=1, keep=discovery)
Xb<- pgs(bfile=bfile, weights=beta, chr=c(1), trace=1, keep=discovery)
# write.table(Xb,'pgs_chr1.txt')
herit <- 0.8
error <- rnorm(m, sd=sqrt(var(Xb)/herit*(1-herit)))
y <- Xb + error
# write.table(y,'phenotype_chr1.txt')
Y <- rep(NA, N)
Y[discovery] <- y

#### Train ####
sum <- plink(bfile=bfile, cmd="--linear", plink2=TRUE, 
             pheno=Y, chr=c(1), silent=F, keep=discovery, 
             mem=fixed.mem)
summ <- read.table.plink(sum, ".PHENO1.glm.linear", header=T)

# est.beta <- summ$BETA
# est.beta[is.na(est.beta)] <- 0
# est.beta <- ifelse(summ$A1 == summ$ALT, est.beta, -est.beta)
# write.table(est.beta,"est.beta.chr1.txt")
# pval <- summ$P
# 
# pval[is.na(pval)] <- 1
# write.table(pval,'pval.chr1.txt')
# 
# out = summ[,c(1,2,3,8,11)]
# write.table(out,'summary.chr1.txt')
est.z <- summ$T_STAT
est.z[is.na(est.z)] <- 0
# est.z <- ifelse(summ$A1 == summ$ALT, est.z, -est.z)
Zscore[[i]]=est.z[2:9]
Zscorefull[[i]]=est.z
print(i)
}
write.table(Zscorefull,'zscore.txt')
ZZZ=as.data.frame(Zscore)
ZZ2=t(ZZZ)
head(ZZ2)
colnames(ZZ1)<-c('SNP1','SNP2','SNP3','SNP4','SNP5')
colnames(ZZ2)<-c('SNP2','SNP3','SNP4','SNP5','SNP6','SNP7','SNP8','SNP9')
corr1=cor(ZZ1)
corr2=cor(ZZ2)
library('corrplot')
corrplot(corr1, method="circle")
corrplot(corr2, method="circle")
corr=cor(data3[,1:8])
corr3=cor(cov)
cov
corrplot(cov, method="circle")
jpeg("rplot.jpg", width = 350, height = "350")
# 2. Create the plot
corrplot(cov, method="circle")
# 3. Close the file
dev.off()
