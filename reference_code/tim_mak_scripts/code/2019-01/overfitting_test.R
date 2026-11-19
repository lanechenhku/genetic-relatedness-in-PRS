Tmisc()
Tim.load(Rplink)
Tim.load(lassosum)
filename <- Rfilename("overfitting_test", seed=1000)
testn <- 7000
trainn <- 7000
usep <- 30000
useWTCCC()
bfile(subsample(p=usep, contiguous = F))
cl <- startParallel(10)
p <- ncol.bfile(.bfile)
n <- nrow.bfile(.bfile)

group <- sample(c(rep(0,trainn), rep(1, testn), rep(2,n - testn - trainn)))
testsample <- group == 1
validsample <- group == 2
trainsample <- group == 0

n.causal <- 1000
herit <- 0.8
beta <- rep(0, p)
beta[sample(p, n.causal)] <- rnorm(n.causal)
Xb <- pgs(.bfile, beta)
error <- rnorm(n, sd = var(Xb)/herit*(1-herit))
pheno <- Xb + error

stats <- plink.table(cmd="--linear", pheno=pheno, silent=F, keep=trainsample)
cor <- p2cor(stats$P,n = sum(trainsample),sign = stats$BETA)
lp <- lassosum.pipeline(cor=cor, chr=stats$CHR, snp=stats$SNP, A1=stats$A1, 
                        test.bfile=.bfile, keep.test=validsample, 
                        exclude.ambiguous = FALSE,
                        LDblocks = "EUR.hg19", cluster=cl)
v <- validate(lp, pheno=pheno[validsample])
est <- pgs(.bfile, v$best.beta)

#### Tests ####
lapply(0:2, function(i) cor(est[group == i], Xb[group == i]))
lapply(0:2, function(i) cor(est[group == i], error[group == i]))
#### Cor with X ####
train.cor <- plink.table(cmd="--linear standard-beta", keep=trainsample, pheno=est)
train.cor0 <- plink.table(cmd="--linear standard-beta", keep=trainsample, pheno=Xb)
test.cor <- plink.table(cmd="--linear standard-beta", keep=testsample, pheno=est)
test.cor0 <- plink.table(cmd="--linear standard-beta", keep=testsample, pheno=Xb)

