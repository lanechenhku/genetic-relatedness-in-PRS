rearrange <- function(string) {
  # string <- "j2a1bcd"
  if(grepl("[0-9][0-9]", string)) {
    stop("powers cannot be greater than 9.")
  }
  if(grepl("[01]", string)) {
    stop("powers must be at least 2.")
  }
  
  s <- strsplit(string, "")[[1]]
  numbers <- grepl("[0-9]", s)
  vars <- sort(unique(s[!numbers]))
  tally <- rep(0, length(vars))
  names(tally) <- vars

  char <- ""  
  for(i in 1:length(s)) {
    # i <- 1
    if(numbers[i]) tally[char] <- tally[char] + as.integer(s[i]) - 1 else {
      char <- s[i]
      tally[char] <- tally[char] + 1
    }
  }
  # Not yet finished... 

}

