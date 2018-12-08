rm(list = ls()) 


library(HardyWeinberg)
library(ggplot2)
library("GenABEL")
library(LaplacesDemon)
library(kernlab)
library(MASS)

##################################
# function to compute observed MAF
# 
# genoCol: a vector of additive coding
##################################

MAF<- function(genoCol){
  count0= sum(genoCol==0)
  count1= sum(genoCol==1)
  count2= sum(genoCol==2)
  MAFtable= c(count0, count1, count2)
  p= maf(MAFtable)
  return(p)
}

############################
# standardize the genotype
############################

Gstand<- function(genoCol){
  count0= sum(genoCol==0)
  count1= sum(genoCol==1)
  count2= sum(genoCol==2)
  MAFtable= c(count0, count1, count2)
  p= maf(MAFtable)
  Gs=(genoCol-2*p)/(sqrt(2*p*(1-p)))
  return(Gs)
}

############################
# true genotype simulation
############################

set.seed(8899)

# simulate the haplotype in msdir
command= paste("/Users/miaozhang/Desktop/msdir/ms 1500 1 -t 935.68 -I 3 500 500 500 -r 935.68 1000 -n 1 1.980002 -n 2 4.914199 -n 3 6.696391 -eg 0 2 110.450397 -eg 0 3 139.447680 -ma x 0.743706 0.227223 0.743706 x 0.911111 0.227223 0.911111 x -ej 0.032140 3 2 -en 0.032140 2 0.254609 -ema 0.032140 3 x 4.457636 x 4.457636 x x x x x -ej 0.070325 2 1 -en 0.202422 1 1  >/Users/miaozhang/Desktop/msdir/ms.out")

system(command)

# read the simulated result by using the author's R function (the R script is in msdir package)
msout <- read.ms.output(file.ms.output="/Users/miaozhang/Desktop/msdir/ms.out" )

haplotype<- msout$gametes[[1]]

# subpopulation samples
n1= 300 # white
n2= 300 # Asian
n3= 200 # black

idx1= sample(1:500, n3, replace=T)
idx2= sample(1:500, n3, replace=T)
G1= haplotype[idx1,]
G2= haplotype[idx2,]
G_Black= G1+G2

idx1= sample(501:1000, n1, replace=T)
idx2= sample(501:1000, n1, replace=T)
G1= haplotype[idx1,]
G2= haplotype[idx2,]
G_White= G1+G2

idx1= sample(1001:1500, n2, replace=T)
idx2= sample(1001:1500, n2, replace=T)
G1= haplotype[idx1,]
G2= haplotype[idx2,]
G_Asian= G1+G2

G= rbind(G_Asian, G_White, G_Black)

########################################################################

library(HardyWeinberg)
#library(dplyr)
#library(tidyr)

##############################################
MAF<- function(genoCol){
  count0= sum(genoCol==0)
  count1= sum(genoCol==1)
  count2= sum(genoCol==2)
  MAFtable= c(count0, count1, count2)
  p= maf(MAFtable)
  return(p)
}
##############################################
########################
# add noise 
########################
# read 1000 Genome quality file
xq = read.table("/home/u11/miaozhang/dosage/xqual.txt", header=TRUE)

# function to add noise based on coverage depth
convert_byxq<-function(g, dp=30, xq, seed=NULL){
  if (! is.null(seed)){
    set.seed(seed)
  }
  n<-nrow(g); p<-ncol(g);
  
  m<-round(rgamma(n=n*p, shape=6.3, scale=dp/6.3))
  M<-array(m, c(n,p))
  
  q<-sapply(m, function(x) {dist<-abs(xq[,1]-x); sample(xq[dist==min(dist),2],1)})
  Q<-array(q, c(n,p))
  er<-10^(-Q/10)
  
  g1<-as.numeric(g>0)
  g2<-as.numeric(g>1)
  er1<-array(rbinom(n=n*p, size=1, prob=c(er)), c(n,p))
  er2<-array(rbinom(n=n*p, size=1, prob=c(er)), c(n,p))
  G<-abs(g1 - er1) + abs(g2 - er2)
  
  list(G=G, M=M, Q=Q, E=er)
}
##############################################
##############################################

noise= convert_byxq(G, dp=5, xq)
E= noise$E # prob for incorrect calling for ALT
P= 1-E  # prob for correct calling for ALT
G_obs= noise$G
DS= G_obs*P # genotype dosage value
Q= noise$Q


idx= which(DS==0) # find those genotype 0 
DS[idx]= 2*E[idx]
write.table(G_obs, file="/home/u11/miaozhang/dosage/Genotype_cov5.txt", col.names=F, row.names=F, quote=F)
write.table(DS, file="/home/u11/miaozhang/dosage/DosageGenotype_cov5.txt", col.names=F, row.names=F, quote=F)
write.table(Q, file="/home/u11/miaozhang/dosage/GQ_cov5.txt", col.names=F, row.names=F, quote=F)
write.table(E, file="/home/u11/miaozhang/dosage/Error_cov5.txt", col.names=F, row.names=F, quote=F)




