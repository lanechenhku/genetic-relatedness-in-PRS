uniLasso.h1 <- read.csv("C:/Users/22367/Desktop/uniLasso-h1.txt", sep="")

rownames( uniLasso.h1 )<-c('FPR','FNR','FDR','power','avg.h2.m','n.causal','prop.causal')
colnames( uniLasso.h1 )<-c(1:10)
uniLasso.h1=data.frame(t( uniLasso.h1 ))
uniLasso.h1$Fscore=2*uniLasso.h1$power*(1-uniLasso.h1$FPR)/(uniLasso.h1$power+(1-uniLasso.h1$FPR))
library(ggplot2)
library(reshape2)
test_data_long1 <- melt(uniLasso.h1[,c(1,3,4,6,8)], id="n.causal") 
p2<-ggplot(data=test_data_long1, aes(x=n.causal, y=value, colour=variable)) + geom_line()+geom_point()# convert to long format ggplot(data=test_data_long, aes(x=date, y=value, colour=variable)) + geom_line()
# fills=c('Fscore'='grey','power'='red','FDR'='blue')
# p1<-ggplot(data=uniLasso.h1)+geom_line(aes(x=avg.h2.m,y=Fscore),color='grey',size=1)+ 
#   geom_line(aes(x=avg.h2.m,y=power,fill='Fscore'),color='red',size=1)+
#   geom_line(aes(x=avg.h2.m,y=FDR),color='blue',size=1)+
#   geom_point(aes(x=avg.h2.m,y=Fscore),color='grey',size=2)+ 
#   geom_point(aes(x=avg.h2.m,y=power),color='red',size=2)+
#   geom_point(aes(x=avg.h2.m,y=FDR),color='blue',size=2)+
#   scale_fill_brewer(limits=c('Fscore','power','FDR')) +
#   theme(plot.title = element_text(hjust = 0.5)) 
# p1
test_data_long <- melt(uniLasso.h1[,c(1,3,4,5,8)], id="avg.h2.m") 
p1<-ggplot(data=test_data_long, aes(x=avg.h2.m, y=value, colour=variable)) + geom_line()+geom_point()+labs(title='h^2 in [0.001,0.25]')
# p2<-ggplot(data=uniLasso.h1, aes(x=n.causal,y=Fscore)) + 
#   geom_line(color='red',size=1)+
#   geom_point()
# p3<-ggplot(data=uniLasso.h1, aes(x=prop.causal,y=Fscore)) + 
#   geom_line(color='blue',size=1)+
#   geom_point()
# p4<-ggplot(data=uniLasso.h1, aes(x=avg.h2.m,y=power)) + 
#   geom_line(color='green',size=1)+
#   geom_point()
# library(cowplot)
# p5 <- cowplot::plot_grid(p1, p2, nrow = 1, labels = LETTERS[1:2])
p1
p2
