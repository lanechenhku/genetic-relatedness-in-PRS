Tmisc()
Tim.load("WishartMoment")
setwd0(attachroot("~/WORK/Projects/relatedness/moments/"))
source("complexMoments.R")
source("possible.arrangements.R")

input <- c("jr", "hr", "hj")
Input2 <- possible.arrangements(input, "r", p.symbol="n")

Result <- list(arrangments=list(), coefTables=list())
for(i in 1:length(Input2$arrangements)) {
  # i <- 1
  arr <- possible.arrangements(Input2$arrangements[[i]], "jh", p.symbol="m")
  Result$arrangments <- c(Result$arrangments, arr$arrangements)
  cT <- lapply(arr$coefTables, function(x) 
    flatten.coefTable(outer(Input2$coefTables[[i]], x)))
  
  Result$coefTables <- c(Result$coefTables, cT)
}
print(complexMoments(Result$arrangments, Result$coefTables))

#' Calculate moments for Tr(ZX'XZ'ZZ' %*% ZX'XZ'ZZ')
#' Let W = ZX'XZ' 
#' W_{jh} = p^2 sum_r r_jr r_hr
#' ZZ'_{hk} = p r_hk
#' WZZ'_jk = p^3 sum_r sum_h r_jr r_hr r_hk  
#' Tr(ZX'XZ'ZZ' %*% ZX'XZ'ZZ') = 
#'   sum_j sum_k sum_h sum_l sum_r sum_s r_jr*r_hr*r_hk*r_js*r_ls*r_lk
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


