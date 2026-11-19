Tmisc()
Tim.load("Rplink")
filename <- Rfilename("UKBBrel", seed = 1000)
Tim.load("UKBBscripts", "~/DATA/UKBB2/")
bfile.fun <- function(i) {
  paste0(attachroot("~/DATA/UKBB2/data/genotype/whiteBritish1_train/whiteBritish1_train"), i)
}

N <- nrow.bfile(bfile.fun(1))
tokeep <- logical.vector(sample(N, 100), N)
p <- rep(NA, 22)
for(i in 1:22) {
  # i <- 22
  bfile <- bfile.fun(i)
  bfile2 <- plink(bfile, keep=tokeep, cmd="--make-bed", maf=0.01, silent=F, plink2 = T)
  print.counter(i, text="Chr")
  p[i] <- ncol.bfile(bfile2)
  X <- plink(bfile2, cmd="--make-rel square") # This is probably not exactly XX'/p but is close enough
  M <- read.table.plink(X, ".rel")
  if(i == 1) MM <- as.matrix(M) * p[i] else 
    MM <- MM + as.matrix(M) * p[i]
}
MM_UKBB <- MM/sum(p)
saveRDS(MM_UKBB, file = paste0(filename, ".rds"))