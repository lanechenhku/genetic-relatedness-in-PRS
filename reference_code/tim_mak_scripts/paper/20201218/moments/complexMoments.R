complexMoments <- function(list.of.input, list.of.multipliers) {
  Input <- list.of.input
  terms <- list.of.multipliers
  
  W <- lapply(Input, WishartMoment, IdentityMatrix=TRUE)
  W2 <- lapply(1:length(Input), 
               function(i) {
                 result <- kronecker(W[[i]], matrix(terms[[i]], 1))
                 rownames(result) <- rownames(W[[i]])
                 colnames(result) <- names(terms[[i]])
                 return(result)
               })
  
  return(do.call(mergeTables, W2))
  
}

