## fixed.mem=40gb
## total.walltime=40:00:00

Tmisc()
Tim.load("WishartMoment")
setwd0(attachroot("~/WORK/Projects/relatedness/moments/"))
filename <- Rfilename("Denominator3")
source("complexMoments.R")
source("possible.arrangements.R")

#' Calculate moments for Tr(WZZ')=Tr(ZX'XZ'ZZ')
#' W_{jh} = p^2 sum_r r_jr r_hr
#' ZZ'_{hk} = p r_hk
#' WZZ'_jj = p^3 sum_r sum_h r_jr r_hr r_hj  
input <- c("jr", "hr", "hj")
A <- possible.arrangements.wrapper(input, groups=c(m="jh", n="r"))
mean <- complexMoments(A$arrangements, A$coefTables)

#' Calculate moments for Tr(ZX'XZ'ZZ' %*% ZX'XZ'ZZ')
#' Let W = ZX'XZ' 
#' W_{jh} = p^2 sum_r r_jr r_hr
#' ZZ'_{hk} = p r_hk
#' WZZ'_jk = p^3 sum_r sum_h r_jr r_hr r_hk  
#' Tr(ZX'XZ'ZZ' %*% ZX'XZ'ZZ') = 
#'   sum_j sum_k sum_h sum_l sum_r sum_s r_jr*r_hr*r_hk*r_js*r_ls*r_lk

input <- c("jr", "hr", "hk", "js", "ls", "lk")
A <- possible.arrangements.wrapper(input, groups=c(m="jkhl", n="rs"))
var <- complexMoments(A$arrangements, A$coefTables)

save(mean, var, file=paste0(filename, ".RData"))
