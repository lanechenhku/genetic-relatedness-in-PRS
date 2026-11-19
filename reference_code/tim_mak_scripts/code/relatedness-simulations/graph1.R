Tmisc()
setwd0(attachroot("~/WORK/Projects/relatedness/simulations/"))
load.results("sim", prompt = F)
filename <- Rfilename("graph1")

library(ggplot2)
t <- results.table
n.causal <- sort(unique(t$n.causal))
t$n.causal <- factor(t$n.causal, n.causal, paste0("n.causal=\n", n.causal))
m <- sort(unique(t$m))
t$m <- factor(t$m, m,paste0("m=", m))
t$type2 <- recode(t$type2, bestclumped100="100", copy = T)
# t$type2 <- factor(t$type2, labels=c(`100`="bestclumped100", `full`="full", 
#                                     `best`="best"))
touse <- subset(t, type2 != "best")
g <- ggplot(touse, aes(x=overlap, y=cor, color=type2)) + 
  facet_grid(as.factor(m) ~ as.factor(n.causal))
g + geom_line()
ggsave(paste0(filename, ".1.pdf"))

touse <- subset(t, T)
g <- ggplot(touse, aes(x=overlap, y=cor, color=type2)) + 
  facet_grid(as.factor(m) ~ as.factor(n.causal))
g + geom_line()
ggsave(paste0(filename, ".2.pdf"))
