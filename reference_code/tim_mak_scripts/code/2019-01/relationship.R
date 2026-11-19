Tmisc()
Tim.load(lassosum)
Tim.load(Rplink)

# WTCCC
bfile <- useWTCCC(F)
tokeep <- logical.vector(sample(nrow.bfile(bfile), 100), nrow.bfile(bfile))
bfile2 <- plink(bfile, keep=tokeep, fill.missing.a2=T, cmd="--make-bed", 
                maf=0.01)
X <- plink(bfile2, cmd="--make-rel square") # This is probably not exactly XX'/p but is close enough
M <- read.table.plink(X, ".rel")
MM_WTCCC <- as.matrix(M)
summary(offdiags(M))

# 1000G GBR
samp <- get1000Gsample(FIDIID = F)
tokeep <- samp$pop=="GBR"
p <- rep(NA, 22)
for(i in 1:22) {
  # i <- 1
  bfile <- use1000G(i,.bfile = F)
  bfile2 <- plink(bfile, keep=tokeep, fill.missing.a2=T, cmd="--make-bed", 
                  maf=0.01)
  print.counter(i, text="Chr")
  p[i] <- ncol.bfile(bfile2)
  X <- plink(bfile2, cmd="--make-rel square") # This is probably not exactly XX'/p but is close enough
  M <- read.table.plink(X, ".rel")
  if(i == 1) MM <- as.matrix(M) * p[i] else 
    MM <- MM + as.matrix(M) * p[i]
}
MM_GBR <- MM/sum(p)

# 1000G EUR
samp <- get1000Gsample(FIDIID = F)
tokeep <- samp$super_pop == "EUR"
tokeep[tokeep] <- logical.vector(sample(sum(tokeep), 100), sum(tokeep))

p <- rep(NA, 22)
for(i in 1:22) {
  # i <- 1
  bfile <- use1000G(i,.bfile = F)
  bfile2 <- plink(bfile, keep=tokeep, fill.missing.a2=T, cmd="--make-bed", 
                  maf=0.01)
  print.counter(i, text="Chr")
  p[i] <- ncol.bfile(bfile2)
  X <- plink(bfile2, cmd="--make-rel square") # This is probably not exactly XX'/p but is close enough
  M <- read.table.plink(X, ".rel")
  if(i == 1) MM <- as.matrix(M) * p[i] else 
    MM <- MM + as.matrix(M) * p[i]
}
MM_EUR <- MM/sum(p)

# 1000G 
samp <- get1000Gsample(FIDIID = F)
tokeep <- logical.vector(sample(nrow(samp), 100), nrow(samp))

p <- rep(NA, 22)
for(i in 1:22) {
  # i <- 1
  bfile <- use1000G(i,.bfile = F)
  bfile2 <- plink(bfile, keep=tokeep, fill.missing.a2=T, cmd="--make-bed", 
                  maf=0.01)
  print.counter(i, text="Chr")
  p[i] <- ncol.bfile(bfile2)
  X <- plink(bfile2, cmd="--make-rel square") # This is probably not exactly XX'/p but is close enough
  M <- read.table.plink(X, ".rel")
  if(i == 1) MM <- as.matrix(M) * p[i] else 
    MM <- MM + as.matrix(M) * p[i]
}
MM_1000G <- MM/sum(p)

1/sd(offdiags(MM_WTCCC))
1/sd(offdiags(MM_GBR))
1/sd(offdiags(MM_EUR))
1/sd(offdiags(MM_1000G))





