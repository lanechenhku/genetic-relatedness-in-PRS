## fixed.mem= 40gb
Tmisc()
Tim.load("WishartMoment")
setwd0(attachroot("~/WORK/Projects/relatedness/moments/"))
source("complexMoments.R")
source("possible.arrangements.R")
filename <- Rfilename("VarWZZT", seed="prompt")
1

which <- attr(filename, "seed")
input <- c("jr", "hr", "hk", "js", "ls", "lk")
Input2 <- possible.arrangements(input, "rs", p.symbol="n")

Result <- list(arrangments=list(), coefTables=list())
for(i in 1:length(Input2$arrangements)) {
  # i <- 1
  arr <- possible.arrangements(Input2$arrangements[[i]], "jkhl", p.symbol="m")
  Result$arrangments <- c(Result$arrangments, arr$arrangements)
  cT <- lapply(arr$coefTables, function(x) 
    flatten.coefTable(outer(Input2$coefTables[[i]], x)))
  
  Result$coefTables <- c(Result$coefTables, cT)
}
# print(complexMoments(Result$arrangments, Result$coefTables))

result <- WishartMoment(Result$arrangments[[which]], IdentityMatrix = TRUE)
multiplier <- Result$coefTables[[which]]
save(result, multiplier, file=paste0(filename, ".RData"))

