## fixed.mem = 40gb
Tmisc()
Tim.load("Rplink")
filename <- Rfilename("UKBBsim.1000", seed = "prompt")
0

Tim.load("UKBBscripts", "~/DATA/UKBB2/")
Tim.load("lassosum")
Tim.load("crosspred", "~/WORK/Projects/validation/")
cl <- startParallel(8)

fixed.mem <- 40000
bfile <- attachroot("~/DATA/UKBB2/data/genotype/whiteBritish1_train/whiteBritish1_train")

N <- nrow.bfile(bfile)
group <- rep("", N)
groups <- c(rep("train.100000", 100000), 
            rep("test.100000", 100000), 
            rep("train.10000", 10000), 
            rep("test.10000", 10000), 
            rep("train.1000", 1000), 
            rep("test.1000", 1000))
group[1:length(groups)] <- groups
group <- sample(group)
n <- as.integer(nameoption(filename, 2))
n.causal <- ceiling(n/10)
p <- ncol.bfile(bfile)

beta <- rep(0, p)
beta[sample(p, n.causal)] <- rnorm(n.causal)
tokeep <- group == paste0("train.", n)
Xb <- pgs(bfile=bfile, weights=beta, cluster=cl, trace=1, 
           keep=tokeep)
herit <- 0.5
error <- rnorm(n, sd=sqrt(var(Xb)/herit*(1-herit)))
y <- Xb + error
Y <- rep(NA, N)
Y[tokeep] <- y
#### Train ####
sum <- plink(bfile=bfile, cmd="--linear", plink2=TRUE, 
             pheno=Y, silent=F, keep=tokeep, threads = length(cl), 
             mem=fixed.mem)
summ <- read.table.plink(sum, ".PHENO1.glm.linear", header=T)
est.beta <- summ$BETA
est.beta[is.na(est.beta)] <- 0
pval <- summ$P
pval[is.na(pval)] <- 1
testkeep <- group == paste0("test.", n)
pl <- pthresh.pipeline(beta=est.beta, pvals=pval, 
                       snp=summ$ID, keep.test = testkeep,
                       A1=summ$ALT, exclude.ambiguous = F, test.bfile = bfile)
Xb.test <- pgs(bfile=bfile, weights=beta, cluster=cl, trace=1, 
          keep=testkeep)
error.test <- rnorm(n, sd=sqrt(var(Xb)/herit*(1-herit)))
y.test <- Xb.test + error.test
v <- validate(pl, pheno=y.test)
best.thres <- which(v$lambda == v$best.lambda)
best.beta <- v$best.beta

o <- order(pval)
est.beta100 <- est.beta
est.beta100[o > 100] <- 0
BETA <- cbind(best.beta, est.beta100, est.beta)
PGS.test <- pgs(bfile=bfile, weights=BETA, cluster=cl, trace=1, 
               keep=testkeep)
PGS.train <- pgs(bfile=bfile, weights=BETA, cluster=cl, trace=1, 
                keep=tokeep)

cor.test <- cor(PGS.test, Xb.test)
cor.train <- cor(PGS.train, Xb)     

fill.in.results.table(n=n, seed=attr(filename, "seed"),
                      cor = cor.train, 
                      type2 = c("best", "0.1", "full"), 
                      type="train",
                      p = colSums(BETA != 0), 
                      P = nrow(BETA))

fill.in.results.table(n=n, seed=attr(filename, "seed"),
                      cor = cor.test, 
                      type2 = c("best", "0.1", "full"), 
                      type="test",
                      p = colSums(BETA != 0), 
                      P = nrow(BETA))

save(results.table, file=paste0(filename, ".RData"))
