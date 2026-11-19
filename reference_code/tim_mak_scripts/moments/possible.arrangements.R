gettable <- function(group) {
  s <- strsplit(group, "")[[1]]
  n <- length(s)
  index <- 1:n
  names(index) <- s
  
  combs <- fact.comb(n)
  res <- do.call(rbind, combs)
  colnames(res) <- s
  return(res)
  
}
  
fact.comb <- function(n) {
  if(n == 1) R <- list(1) else {
    R <- list()
    for(i in 1:n) {
      r <- fact.comb(n-1)
      r <- lapply(r, function(x) return(c(x, i)))
      R <- c(R, r)
    }
  }
  return(R)
}


# group <- "jkhl"


possible.arrangements <- function(input, group, p.symbol) {
  # input <- c("jr", "hr", "hk", "js", "ls", "lk"); split <- "jkhl"
  
  Arr <- list()
  table <- gettable(group)
  coefTables <- list()
  index <- colnames(table)
  
  for(j in 1:nrow(table)) {
    ar <- input
    for(l in 1:ncol(table)) {
      ar <- gsub(index[l], index[table[j,l]], ar)
    }
    Arr[[j]] <- ar
    coefTables[[j]] <- coeffun(length(unique(table[j,])), p.symbol=p.symbol)
  }

  return(list(arrangements=Arr, coefTables=coefTables))
  
}

flatten.coefTable <- function(tab) {
  tab2 <- as.data.frame(as.table(tab))
  res <- tab2$Freq
  names(res) <- paste0(tab2$Var1, tab2$Var2)
  return(res)
}

possible.arrangements.wrapper <- function(input, groups) {
  # input <- c("jr", "hr", "hk", "js", "ls", "lk"); groups <- c("jkhl","rs")
  
  # Checks 
  s <- strsplit(paste0(groups, collapse = ""), "")[[1]]
  if(length(s) > length(unique(s))) {
    stop("groups must be mutually exclusive.")
  }

  if(!is.list(input)) {
    input <- list(input)
  }
  
  Result <- list(arrangements=list(), coefTables=list())
  names <- names(groups)
  if(is.null(names)) names <- letters[length(groups):1]
  
  if(length(groups) == 1) {
    for(i in 1:length(input)) {
      arr <- possible.arrangements(input[[i]], groups, p.symbol=names[1])
      Result$arrangements <- c(Result$arrangements, arr$arrangements)
      cT <- arr$coefTables
      Result$coefTables <- c(Result$coefTables, cT)
    }
  } else if(length(groups) >= 2) {
    res <- possible.arrangements.wrapper(input, groups[-1])
    for(i in 1:length(res$arrangements)) {
      arr <- possible.arrangements(res$arrangements[[i]], groups[1], 
                                   p.symbol = names[1])
      Result$arrangements <- c(Result$arrangements, arr$arrangements)
      cT <- lapply(arr$coefTables, function(x)
        flatten.coefTable(outer(res$coefTables[[i]], x)))
      Result$coefTables <- c(Result$coefTables, cT)
    }
  } 
  return(Result)
}
