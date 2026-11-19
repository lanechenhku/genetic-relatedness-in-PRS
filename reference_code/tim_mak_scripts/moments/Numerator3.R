Tmisc()
Tim.load("WishartMoment")
setwd0(attachroot("~/WORK/Projects/relatedness/moments/"))
filename <- Rfilename("Numerator3")
source("complexMoments.R")
source("possible.arrangements.R")

#' Calculate moments for Tr(W)=Tr(ZX'XZ')
#' Let W = ZX'XZ' 
#' W_{jj} = p^2 sum_r r_{jr}^2
#' E(Tr(W)) = p^2 sum_r sum_j E(r_{jr}^2) = m*n*p^2*E(r_{jr}^2)
input <- c("jr", "jr")
A <- possible.arrangements.wrapper(input, groups=c(m="j", n="r"))
mean <- complexMoments(A$arrangements, A$coefTables)


#' Calculate E(Tr(W^2))=E(Tr(ZX'XZ'ZX'XZ'))
#' Tr(W^2) = sum_j sum_k W_{jk}^2
#' W_{jk} = p^2 sum_r r_{jr}r_{kr}
#' (W_{jk})^2 = p^4 sum_r sum_s r_{jr}r_{kr}r_{js}r_{ks}

input <- c("jr", "kr", "js", "ks")
A <- possible.arrangements.wrapper(input, groups=c(m="jk", n="rs"))
var <- complexMoments(A$arrangements, A$coefTables)

save(mean, var, file=paste0(filename, ".RData"))