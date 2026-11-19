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
bfile<-'/psychipc01/disk3/data/UKBB/data/genotype/imputed_0.01/maf_0.01_1'
bfile<-"/psychipc01/disk3/data/UKBB/data/genotype/whiteBritish1_train/whiteBritish1_train1"
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
beta <- rep(0, p)
beta[sample(p, n.causal)] <- rnorm(n.causal)
write.table(beta,'effectsize_chr1.txt')
#Xb <- pgs(bfile=bfile, weights=beta, cluster=cl, trace=1, keep=discovery)
# Xb <- pgs(bfile=bfile, weights=beta, chr=c(1,2), trace=1, keep=discovery)
Xb<- pgs(bfile=bfile, weights=beta, chr=c(1), trace=1, keep=discovery)
write.table(Xb,'pgs_chr1.txt')
herit <- 0.8
error <- rnorm(m, sd=sqrt(var(Xb)/herit*(1-herit)))
y <- Xb + error
write.table(y,'phenotype_chr1.txt')
Y <- rep(NA, N)

Y[discovery] <- y

#### Train ####
sum <- plink(bfile=bfile, cmd="--linear", plink2=TRUE, 
             pheno=Y, chr=c(1), silent=F, keep=discovery, 
            mem=fixed.mem)
summ <- read.table.plink(sum, ".PHENO1.glm.linear", header=T)

est.beta <- summ$BETA
est.beta[is.na(est.beta)] <- 0
est.beta <- ifelse(summ$A1 == summ$ALT, est.beta, -est.beta)
write.table(est.beta,"est.beta.chr1.txt")
pval <- summ$P

pval[is.na(pval)] <- 1
write.table(pval,'pval.chr1.txt')

out = summ[,c(1,2,3,8,11)]
write.table(out,'summary.chr1.txt')

#### extract
# sum <- plink(bfile=bfile, cmd="--snps 10583-1892607 \ --make-bed \
#  --out block_1", plink2=TRUE, 
#              chr=c(1), silent=F, keep=discovery, 
#              mem=fixed.mem)
# y1=read.table('phenotype_chr1.txt')
# y=y1$V1
library('KRIS')
aa=c('/home/unfated/chr1_block/block_')
i=1
bed=paste0(aa,'.bed')
bim=paste0(aa,'.bim')
fam=paste0(aa,'.fam')
# bed=paste0(aa,i,'.bed')
# bim=paste0(aa,i,'.bim')
# fam=paste0(aa,i,'.fam')
# block1=read.bed(bed,bim,fam)
block=read.bed(bed,bim,fam)
z_1=block1$snp[discovery,]
data1=data.frame(y,z_1)
id=c()
for(i in 1:ncol(data1)){
  if(sum(is.na(data1[,i]))>50000)
  id=c(id,i)
}
data2=data1[,-id]

# id2=c()
# for(i in 1:nrow(data2)){
#   c=sum(is.na(data2[i,]))
#     id2[i]=c
# }
# max(id2)
#### regression
sum(is.na(data2))

# give 1 to NA
data3=data2
for(j in 1:ncol(data3)){
data3[is.na(data3[,j]),j]=1
}
#data3[1,1:50]
fit1<-lm(y~.,data=scale(data3))
#give 2 to NA

data4=data2
for(j in 1:ncol(data4)){
  data4[is.na(data4[,j]),j]=2
}
#data4[1,1:50]
fit2<-lm(y~.,data=data4)

#summary(fit)
beta.fit1=fit1$coefficients
beta.fit2=fit2$coefficients
beta.fit1[1:10]
beta.fit2[1:10]
write.table(beta.fit1,'glm.beta.imput_1.txt')
write.table(beta.fit2,'glm.beta.imput_2.txt')

### blasso
library('monomvn')
sum(is.na(data4))


### lasso
# library('lars')
# ls_fit2<-lars(as.matrix(data4[,2:5]), y, type = c("lasso", "lar", "forward.stagewise", "stepwise"), 
#      trace = FALSE, normalize = TRUE, intercept = TRUE, eps = .Machine$double.eps, max.steps=100000)
# summary(ls_fit2)
# ls_fit2$beta
# ls(ls_fit2)
library('glmnet')
x2=as.matrix(data4[,-1])
cv2 = cv.glmnet(x2, y)
beta.ls_fit2=coef(cv2, s = "lambda.1se")

x1=as.matrix(data3[,-1])
cv1 = cv.glmnet(x1, y)
beta.ls_fit1=coef(cv1, s = "lambda.1se")

write.table(as.matrix(beta.ls_fit2),'beta.ls_fit2.txt')
write.table(as.matrix(beta.ls_fit1),'beta.ls_fit1.txt')

#model = glmnet(x, y, type.gaussian="covariance", lambda=cv$lambda.1se)
#predict(model, type="coefficients")

save.image()
####

# overlap <- (0:10)/10
# 
# pl <- pthresh.pipeline(beta=est.beta, pvals=pval, 
#                        snp=summ$ID, keep.test = validation,
#                        A1=summ$ALT, exclude.ambiguous = F, 
#                        test.bfile = bfile, clump = TRUE, 
#                        cluster=cl, trace=1, 
#                        clump.options = list(r2=0.1))
# 
# Xb.validation <- pgs(bfile=bfile, weights=beta, cluster=cl, trace=1, 
#                      keep=validation)
# error.validation <- rnorm(validation.n, sd=sqrt(var(Xb)/herit*(1-herit)))
# y.validation <- Xb.validation + error.validation
# v <- validate(pl, pheno=y.validation)
# best.thres <- which(v$lambda == v$best.lambda)
# best.beta <- v$best.beta
# 
# for(i in 1:length(overlap)) {
#   # i <- 1
#   test.sample <- integer(0)
#   overlap.n <- test.n * overlap[i]
#   others.n <- test.n  - overlap.n
#   test.sample <- c(sample(discovery.sample, overlap.n), 
#                    sample(others.sample, others.n))
#   test <- logical.vector(test.sample, N)
#   
#   zero <- pl$beta[[1]]$beta == 0
#   pval2 <- pval
#   pval2[zero] <- 1
#   
#   o <- rank(pval2)
#   est.beta100 <- est.beta
#   est.beta100[o > 100] <- 0
#   BETA <- cbind(beta, best.beta, est.beta100, est.beta)
#   
#   PGS.test <- pgs(bfile=bfile, weights=BETA, cluster=cl, keep=test)
#   
#   cor.test <- as.vector(cor(PGS.test)[1,-1])
#   
#   fill.in.results.table(n=test.n, seed=attr(filename, "seed"),
#                         m=m, n.causal=n.causal, 
#                         cor = cor.test, 
#                         type2 = c("best", "100", "full"), 
#                         p = colSums(BETA[,-1] != 0), 
#                         P = nrow(BETA[,-1]), 
#                         overlap=overlap[i])
# 
# }
# 
# 
# save(results.table, file=paste0(filename, ".RData"))
