if(!("devtools" %in% rownames(installed.packages()))) install.packages("devtools"); .rootdir <- "/home/unfated/Tim"; source("/home/unfated/Tim/Tmisc/inst/Startup.R")
#17
## fixed.mem = 40gb
# a package
Tmisc()
Tim.load("Rplink",lib.dir="/home/unfated/Tim")
library("RcppArmadillo")
library("data.table")
library("Matrix")
set.seed(1)
N=353465
discovery.sample <-sample(N,100000)
discovery <- logical.vector(discovery.sample, N)
itr=c(7:10)
option=c('lasso')
holder=c('/home/unfated/Tim/')
input=c('/home/unfated/chr1_block/imputed_block_')
input=c('/home/unfated/chr1_block/block_')
impute=c('1')
y=read.table('/home/unfated/Tim/phenotype_chr1.txt')
y=y$V1
length(y)
# Lasso(itr=itr,y=y,option=option,discovery=discovery,holder=holder)
Lasso(itr=itr,y=y,option=option,discovery=discovery,holder=holder,impute='0')
Lasso(itr=itr,y=y,option=option,discovery=discovery,holder=holder,impute='mean')
Lasso(itr=itr,y=y,option=option,discovery=discovery,holder=holder,impute='median')
Lasso(itr=itr,y=y,option=option,discovery=discovery,holder=holder,impute='KNN_med')
Lasso(itr=itr,y=y,option=option,discovery=discovery,holder=holder,impute='KNN_avg')
Lasso <- function(holder,itr,option,discovery,y,impute) {
  library('KRIS')
  library('glmnet')
  library('monomvn')
  library('DMwR')
  IID=list()
  for(i in itr){
  data1<-NA
  data3<-NA
  bed=paste0(input,i,'.bed')
  bim=paste0(input,i,'.bim')
  fam=paste0(input,i,'.fam')
  block=read.bed(bed,bim,fam)
  z_1=block$snp[discovery,]
  data1=data.frame(z_1)
  id=c()
  for(k in 1:ncol(data1)){
    if(sum(is.na(data1[,k]))>50000)
      id=c(id,k)
  }
  IID[[i]]=id
  data3=data1
  if(length(id>0))
    data3=data1[,-id]
  if(impute=='1'){
  for(j in 1:ncol(data3)){
    data3[is.na(data3[,j]),j]=1
  }}
  if(impute=='2'){
    for(j in 1:ncol(data3)){
      data3[is.na(data3[,j]),j]=2
    }}
  if(impute=='mean'){
    for(j in 1:ncol(data3)){
      me=mean(data3[,j],na.rm=T)
      data3[is.na(data3[,j]),j]=me
    }}
  if(impute=='median'){
    for(j in 1:ncol(data3)){
      md=median(data3[,j],na.rm=T)
      data3[is.na(data3[,j]),j]=md
    }}
  if(impute=='0'){
    for(j in 1:ncol(data3)){
      data3[is.na(data3[,j]),j]=0
    }}
  if(impute=='KNN_avg'){
    data2=knnImputation(data3, k = 5, scale = T, meth = "weighAvg", distData = NULL)
    data3=data2}
  if(impute=='KNN_med'){
    data2=knnImputation(data3, k = 5, scale = T, meth = "median", distData = NULL)
    data3=data2}
  x1=as.matrix(data3)
  cv1 = cv.glmnet(x1, y)
  if(option=='lasso'||option=='both'){
          beta.ls_fit1=coef(cv1, s = "lambda.1se")
         write.table(as.matrix(beta.ls_fit1),paste0(holder,'beta_',impute,'.ls_fit',i,'.txt'))}
  if(option=='blasso'||option=='both' ){
          bls<-blasso(x1, y, T = 500, thin = NULL, RJ = TRUE, M = NULL,
         beta = NULL, lambda2 = 1, s2 = var(y-mean(y)),
         case = c("default"), mprior = 0, rd = NULL,
         ab = NULL, theta = 0, rao.s2 = TRUE, icept = TRUE, 
         normalize = TRUE, verb = 1)
          
  write.table(as.matrix(beta.ls_fit1),paste0(holder,'beta_',impute,'.bls_fit',i,'.txt'))}
  }}

  