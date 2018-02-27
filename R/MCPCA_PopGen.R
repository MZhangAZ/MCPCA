
###########################################################################################################################################
# this script is to show how apply MCPCA (Soheil Feizi and David Tse, Maximally Correlated Principal Component Analysis, arXiv:1702.05471)
# on genotype dosage value data to analyze the population structure
# In order to run the R script, you also need to download the Matlab script from Soheil Feizi
# https://github.com/SoheilFeizi/MCPCA.git
#
###########################################################################################################################################


library(arules)
library(R.matlab)
library(ggplot2)
# read dosage data

DS= read.table("/Users/miaozhang/Documents/MATLAB/R/DS.txt", header=F)

dis_method= "interval"
categories= 10

######################################################################
# MCPCA_PopGen
# function to apply MCPCA on dosage value data for population inference
# Input
# DS: genotype dosage value data matrix
# dis_method: methods to discretize dosage value, "interval", "frequency",
# "cluster" and "fix"
# categories: discretized categories
#
# Output
# G_Phi: transformed dosage data
######################################################################

MCPCA_PopGen<- function(DS, dis_method= "interval", categories= 10){

discretization<- function(xCol){
  cat<-discretize(xCol, method= dis_method, categories = categories)
  return(as.numeric(cat))
}

DS= data.frame(DS)
DS_clean= DS[vapply(DS, function(x) length(unique(x)) > 1, logical(1L))]
DS_Dis= apply(DS_clean, 2, discretization)

DS_Dis[which(is.na(DS_Dis))]= DS_Dis[which(is.na(DS_Dis))+1]

write.table(DS_Dis, file="/Users/miaozhang/Documents/MATLAB/R/DS_Dis.txt",  quote = F, row.names = F, col.names = F)

#write matlab command line
matlab.lines <- c(

  "C= load('/Users/miaozhang/Documents/MATLAB/R/DS_Dis.txt');",
  "[n,p] = size(C);",
  "q=5;",
  "num_init=0;",
  "num_iter=10;",

  "oldpath=path;",
  "path(oldpath, '/Users/miaozhang/Documents/MATLAB/MCPCA');",
  "[phi_mat,fun_cell]=MCPCA_sample_disc_wrapper(C,q,num_iter,num_init);",
  "X_mc=normalize_matrix(phi_mat);",

  "save('/Users/miaozhang/Documents/MATLAB/R/DS_Phi.mat', 'X_mc')"
)

#create a MATLAB script containing all the commands in matlab.lines
writeLines(matlab.lines, con="/Users/miaozhang/Documents/MATLAB/R/myscript.m")

#run our MATLAB script
system("/Applications/MATLAB_R2016a.app/bin/matlab -nodisplay -r \"run('/Users/miaozhang/Documents/MATLAB/R/myscript.m'); exit\"")

# read MATLAB output in R
data <- readMat("/Users/miaozhang/Documents/MATLAB/R/DS_Phi.mat")
G_Phi= data$X.mc

return(G_Phi)
}

G_Phi= MCPCA_PopGen(DS, dis_method= "interval", categories= 10)

############################################
# plot top 2 PCs for clustering
#
############################################

mcpca= prcomp(G_Phi, center = F, scale. = F)
percentage= (mcpca$sdev^2)/sum(mcpca$sdev^2)

title <- paste("MCPCA: PC1"," (",signif(percentage[1], digits=3)*100,"%)"," / PC2"," (",signif(percentage[2], digits=3)*100,"%)",sep="",collapse="")

group= c(rep("White",50), rep("Asian", 50), rep("Black", 50))

df= data.frame(scale(mcpca$x[,1:2]), group)

g <- ggplot(df, aes(PC1, PC2))
g<- g + geom_point(aes(colour = group, shape = group), size = 2)
g<- g + ggtitle(title)
g<- g + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
g

