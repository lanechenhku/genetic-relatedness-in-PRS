## fixed.mem = 40gb
Tmisc()
Tim.load("Rplink")
filename <- Rfilename("sim.10000.1000", seed = "prompt")
0

test.n <- 2000
validation.n <- 2000
m <- as.integer(nameoption(filename, 2))
n.causal <- as.integer(nameoption(filename, 3))

Tim.load("UKBBscripts", "~/DATA/UKBB2/")
Tim.load("lassosum")
Tim.load("crosspred", "~/WORK/Projects/validation/")
cl <- startParallel(10)

fixed.mem <- 40000
bfile <- attachroot("~/DATA/UKBB2/data/genotype/whiteBritish1_train/whiteBritish1_train")
# bfile <- useWTCCC(.bfile = F)

N <- nrow.bfile(bfile)
discovery.sample <- sample(N, m)
others.sample <- setdiff(1:N, discovery.sample)
validation.sample <- sample(others.sample, validation.n)
others.sample <- setdiff(1:N, c(discovery.sample, validation.sample))

discovery <- logical.vector(discovery.sample, N)
validation <- logical.vector(validation.sample, N)

p <- ncol.bfile(bfile)

beta <- rep(0, p)
beta[sample(p, n.causal)] <- rnorm(n.causal)
Xb <- pgs(bfile=bfile, weights=beta, cluster=cl, trace=1, keep=discovery)
herit <- 0.8
error <- rnorm(m, sd=sqrt(var(Xb)/herit*(1-herit)))
y <- Xb + error
Y <- rep(NA, N)
Y[discovery] <- y

#### Train ####
sum <- plink(bfile=bfile, cmd="--linear", plink2=TRUE, 
             pheno=Y, silent=F, keep=discovery, 
             threads = length(cl), mem=fixed.mem)
summ <- read.table.plink(sum, ".PHENO1.glm.linear", header=T)

est.beta <- summ$BETA
est.beta[is.na(est.beta)] <- 0
est.beta <- ifelse(summ$A1 == summ$ALT, est.beta, -est.beta)

pval <- summ$P
pval[is.na(pval)] <- 1

overlap <- (0:10)/10

pl <- pthresh.pipeline(beta=est.beta, pvals=pval, 
                       snp=summ$ID, keep.test = validation,
                       A1=summ$ALT, exclude.ambiguous = F, 
                       test.bfile = bfile, clump = TRUE, 
                       cluster=cl, trace=1, 
                       clump.options = list(r2=0.1))

Xb.validation <- pgs(bfile=bfile, weights=beta, cluster=cl, trace=1, 
                     keep=validation)
error.validation <- rnorm(validation.n, sd=sqrt(var(Xb)/herit*(1-herit)))
y.validation <- Xb.validation + error.validation
v <- validate(pl, pheno=y.validation)
best.thres <- which(v$lambda == v$best.lambda)
best.beta <- v$best.beta

for(i in 1:length(overlap)) {
  # i <- 1
  test.sample <- integer(0)
  overlap.n <- test.n * overlap[i]
  others.n <- test.n  - overlap.n
  test.sample <- c(sample(discovery.sample, overlap.n), 
                   sample(others.sample, others.n))
  test <- logical.vector(test.sample, N)
  
  zero <- pl$beta[[1]]$beta == 0
  pval2 <- pval
  pval2[zero] <- 1
  
  o <- rank(pval2)
  est.beta100 <- est.beta
  est.beta100[o > 100] <- 0
  BETA <- cbind(beta, best.beta, est.beta100, est.beta)
  
  PGS.test <- pgs(bfile=bfile, weights=BETA, cluster=cl, keep=test)
  
  cor.test <- as.vector(cor(PGS.test)[1,-1])
  
  fill.in.results.table(n=test.n, seed=attr(filename, "seed"),
                        m=m, n.causal=n.causal, 
                        cor = cor.test, 
                        type2 = c("best", "100", "full"), 
                        p = colSums(BETA[,-1] != 0), 
                        P = nrow(BETA[,-1]), 
                        overlap=overlap[i])

}


save(results.table, file=paste0(filename, ".RData"))
