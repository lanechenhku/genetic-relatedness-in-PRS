Tmisc()
Tim.load("WishartMoment")
setwd0(attachroot("~/WORK/Projects/relatedness/moments/"))
source("complexMoments.R")

#' Calculate moments for Tr(ZX'XZ'ZZ')
#' Let W = ZX'XZ' 
#' W_{jh} = p^2 sum_r r_jr r_hr
#' ZZ'_{hk} = p r_hk
#' WZZ'_jj = p^3 sum_r sum_h r_jr r_hr r_hj  
#'         = p^3 [n(m-1) E(r_jr r_hr r_hj) + n E(r_jr r_jr r_jj) ]

input <- c("jr", "hr", "hj")
Input <- list(input, 
              gsub("h", "j", input))

multipliers <- list(c(mn=1, n=-1), 
              c(n=1))
print(complexMoments(Input, multipliers))

#' Calculate moments for Tr(ZX'XZ'ZZ' %*% ZX'XZ'ZZ')
#' Let W = ZX'XZ' 
#' W_{jh} = p^2 sum_r r_jr r_hr
#' ZZ'_{hk} = p r_hk
#' WZZ'_jk = p^3 sum_r sum_h r_jr r_hr r_hk  
#' Tr(ZX'XZ'ZZ' %*% ZX'XZ'ZZ') = 
#'   sum_j sum_k sum_h sum_l sum_r sum_s r_jr*r_hr*r_hk*r_js*r_ls*r_lk
input <- c("jr", "hr", "hk", "js", "ls", "lk")

input_hj <- gsub("h", "j", input)
input_kj <- gsub("k", "j", input)
input_lj <- gsub("l", "j", input) 
input_hk <- gsub("h", "k", input)
input_lh <- gsub("l", "h", input)
input_lk <- gsub("l", "k", input) 

input_khj <- gsub("k", "j", input_hj)
input_lkj <- gsub("l", "j", input_kj)
input_lhj <- gsub("l", "j", input_hj)
input_lhk <- gsub("l", "k", input_hk)

input_lkhj <- gsub("l", "j", input_khj)

input_sr <- gsub("s", "r", input)


gsub("k", "j", input), 
gsub("l", "j", input), )
Input <- c(input, 

