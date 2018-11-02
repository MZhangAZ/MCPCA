#######################################################################
# this script shows the example of using MCPCA to do dimensionality
# reduction on dosage value data
#
# In order to run MCPCA, need to dwonload matlab script from Soheil Feizi
# https://github.com/SoheilFeizi/MCPCA.git
########################################################################

library(HardyWeinberg)
library(ggplot2)
library(arules)
library(R.matlab)
library(grDevices)
library(optparse)
library(stringr)
library(dplyr)
library(tidyr)

############################################################################
# MCPCA_PopGen
# Input
# DS: data matrix of dosage value
# dis_method: discretize methods, "interval", "frequency", "cluster", "fixed"
# categoris: discretize categories
#
# Output
# G_Phi: transformed data matrix
# PC: top 10 projected features from MCPCA
############################################################################

option_list = list(
  make_option("--method", action="store", default=NA, type="character",
              help="discretize method", metavar="character"),
  make_option("--category", action="store", default=NA, type="double",
              help="category for non-FD method", metavar="character"),
  make_option("--q", action="store", default=NA, type="double",
              help="KyFan norm sum", metavar="character"),
  make_option("--iter", action="store", default=NA, type="double",
              help="iteration time", metavar="character"),
  make_option("--subsample", action="store", default=NA, type="character",
              help="sub sample size for each population", metavar="character"),
  make_option("--DSdata", type="character", default=NULL, 
              help="dosage data matrix", metavar="character"),
  make_option("--DS_Dis", action="store", default=NA, type="character",
              help="discretized dosage data matrix", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

set.seed(19880909)
n1= as.numeric(str_split(opt$subsample, "-")[[1]][1])
n2= as.numeric(str_split(opt$subsample, "-")[[1]][2])
n3= as.numeric(str_split(opt$subsample, "-")[[1]][3])

idx1= sample(1:500, n1)
idx2= sample(501:1000, n2)
idx3= sample(1001:1500, n3)

DS= data.matrix(read.table(file=opt$DSdata, header=F))

n= dim(DS)[1]
p= dim(DS)[2]

#############################
# discretize dosage 
#############################
DS_Dis= matrix(0, n, p)

# FD bin method
# equal interval/freq and k-means
if (opt$method=="FD"){
  bin= round( apply(DS, 2, function(x) nclass.FD(x)))
  
  for (i in 1:p){
    DS_Dis[,i]<-discretize(DS[,i], method="interval", categories = bin[i], labels = c(1:bin[i]))
  }} else {
    for (i in 1:p){
      DS_Dis[,i]<-discretize(DS[,i], method=opt$method, categories = as.numeric(opt$category), labels = c(1:as.numeric(opt$category)))
    }
  }



DS_Dis= DS_Dis[c(idx1, idx2, idx3),]


write.table(DS_Dis, file="/home/u11/miaozhang/dosage/DS_temp.txt",  quote = F, row.names = F, col.names = F)

#make a vector where each element is a line of MATLAB code
#matlab code reads in our variable x, creates two variables y and z,
#and write z in a csv file
matlab.lines <- c(
  
  "C= load('/home/u11/miaozhang/dosage/DS_temp.txt');",
  "[n,p] = size(C);",
  paste("q=", opt$q, ";", sep=""),
  "num_init=0;",
  paste("num_iter=", opt$iter, ";", sep=""),
  
  "oldpath=path;",
  "path(oldpath, '/home/u11/miaozhang/dosage/matlab');",
  "[phi_mat,fun_cell]=MCPCA_sample_disc_wrapper(C,q,num_iter,num_init);",
  "X_mc=normalize_matrix(phi_mat);",
  
  "save('/home/u11/miaozhang/dosage/DS_Phi_temp.mat', 'X_mc')"
)

#create a MATLAB script containing all the commands in matlab.lines
writeLines(matlab.lines, con="/home/u11/miaozhang/dosage/myscript.m")

#run our MATLAB script
system("/cm/shared/uaapps/matlab/r2017b/bin/matlab -nodisplay -r \"run('/home/u11/miaozhang/dosage/myscript.m'); exit\"")

# read MATLAB output in R
data <- readMat("/home/u11/miaozhang/dosage/DS_Phi_temp.mat")
G_Phi= data$X.mc

write.table(G_Phi, file=opt$DS_Dis,  quote = F, row.names = F, col.names = F)


############################################################################

Rscript --vanilla /home/u11/miaozhang/dosage/Rscript/DS_DisScript.R \
--q 20 \
--iter 10 \
--method interval \
--subsample 100-100-100 \
--category 10 \
--DSdata /home/u11/miaozhang/dosage/DosageGenotype_cov5.txt \
--DS_Dis /home/u11/miaozhang/dosage/FD-test.txt
##################################################################################

GT= read.table(file="/Users/miaozhang/Documents/MATLAB/Chapter3/Data/Genotype_T_n=1500.txt", header=F)
set.seed(19880909)
idx1= sample(1:500, 100)
idx2= sample(501:1000, 100)
idx3= sample(1001:1500, 100)

GT= GT[c(idx1, idx2, idx3),]

n= dim(GT)[1]
p= dim(GT)[2]

S1= read.table(file="/Users/miaozhang/Desktop/dissertationcomments/Simulation/K8.txt", header=F)
S2= read.table(file="/Users/miaozhang/Desktop/Bin/data/FD-test.txt", header=F)
S3= read.table(file="/Users/miaozhang/Desktop/Bin/data/Genotype_cov5.txt", header=F)
S3= S3[c(idx1, idx2, idx3),]


pca= prcomp(GT, center = T, scale. = T) 
pca_S1= prcomp(S1, center = T, scale. = T) 
pca_S2= prcomp(S2, center = T, scale. = T) 
pca_S3= prcomp(S3, center = T, scale. = T) 


eigenpca= (pca$sdev)^2
eigenSum_pca= sum(eigenpca)
eigenS1= (pca_S1$sdev)^2
eigenSum_S1= sum(eigenS1)
eigenS2= (pca_S2$sdev)^2
eigenSum_S2= sum(eigenS2)
eigenS3= (pca_S3$sdev)^2
eigenSum_S3= sum(eigenS3)

#######################################################################################
dat= data.frame(cbind(seq(1:min(n,p)), eigenpca, eigenS1, eigenS2, eigenS3))
colnames(dat)= c("q", "PCA", "MCPCA_S1", "MCPCA_S2", "MCPCA_S3")
dat= mutate(dat, ratioPCA = cumsum(PCA)/eigenSum_pca, ratioMCPCAs1 = cumsum(MCPCA_S1)/eigenSum_S1,
            ratioMCPCAs2 = cumsum(MCPCA_S2)/eigenSum_S2,
            ratioMCPCAs3 = cumsum(MCPCA_S3)/eigenSum_S3)
dat= select(dat, q, ratioPCA:ratioMCPCAs3)
colnames(dat)= c("q", "PCA", "MCPCA_S1", "MCPCA_S2", "MCPCA_S3")
dat= gather(dat, method, ratio, PCA:MCPCA_S3 )

g <- ggplot(dat, aes(x=q, y=ratio)) +
  geom_line(aes(linetype=method, col= method)) + ylab("Fraction of explained variance") + xlab("Eigenvalue")
g<- g + theme(axis.text=element_text(size=12,face="bold"),
              axis.title=element_text(size=14,face="bold"))
g<- g + scale_colour_manual(values=c(PCA="#339999",MCPCA_S1="#FF6699",MCPCA_S2="#CC0033",MCPCA_S3="#FF6666"),name="Method") + 
  scale_linetype_manual(values = c(rep("dashed", 3),"solid"), name= "Method")  
g<- g + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
g<- g + theme(legend.title = element_text(size=12, face="bold")) # Title appearance
g<- g + theme(legend.text = element_text(size = 10, face = "bold")) # Label appearance
g<- g + theme(legend.background = element_rect(fill="gray90", size=.5, linetype="dotted"))
g<- g + geom_vline(xintercept = 50, colour="green", linetype = "longdash")
g


#######################################################################################
PCAplot<- function(genotype){
  pca= prcomp(genotype, center = T, scale. = T) 
  top_PC= pca$x[,1:2]
  Pop= rep(c("S1","S2","S3"), each=100)
  df= data.frame(top_PC, Pop)
  g= ggplot(df, aes(PC1, PC2, color = as.factor(Pop)))+ geom_point()
  g
}
