Tmisc()
Tim.load("WishartMoment")

#' Calculate moments for Tr(W)=Tr(ZX'XZ')
#' Let W = ZX'XZ' 
#' W_{jj} = p^2 sum_r r_{jr}^2
#' E(Tr(W)) = p^2 sum_r sum_j E(r_{jr}^2) = m*n*p^2*E(r_{jr}^2)
A <- WishartMoment(c("jr", "jr"), IdentityMatrix = TRUE)
#' = 1/p
#' Therefore, E(Tr(W)) = m*n*p

#' Calculate E(Tr(W^2))=E(Tr(ZX'XZ'ZX'XZ'))
#' Tr(W^2) = sum_j sum_k W_{jk}^2
#' W_{jk} = p^2 sum_r r_{jr}r_{kr}
#' E(Tr(W^2)) = m(m-1)n(n-1) E(r_jr r_kr r_js r_ks) + 
#'              mn(n-1) E(r_jr r_jr r_js r_js) + 
#'              m(m-1)n E(r_jr r_kr r_jr r_kr) + 
#'              mn      E(r_jr r_jr r_jr r_jr)

input <- c("jr", "kr", "js", "ks")
Input <- list(input, 
              input1 <- gsub("k", "j", input), 
              gsub("s", "r", input), 
              gsub("s", "r", input1))
W <- lapply(Input, WishartMoment, IdentityMatrix=TRUE)

terms <- list(c(m2n2=1, m2n=-1, mn2=-1, mn=1), 
              c(mn2=1, mn=-1), 
              c(m2n=1, mn=-1), 
              c(mn=1))

W2 <- lapply(1:length(Input), 
             function(i) {
               result <- kronecker(W[[i]], matrix(terms[[i]], 1))
               rownames(result) <- rownames(W[[i]])
               colnames(result) <- names(terms[[i]])
               return(result)
             })

print(do.call(mergeTables, W2))

colnames(A) <- r

