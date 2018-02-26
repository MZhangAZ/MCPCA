
library(arules)
library(R.matlab)
library(ggplot2)
# read dosage data

DS= read.table("/Users/miaozhang/Documents/MATLAB/R/DS.txt", header=F)

dis_method= "interval"
categories= 10

discretization<- function(xCol){
  cat<-discretize(xCol, method= dis_method, categories = categories)
  return(as.numeric(cat))
}

DS_Dis= apply(DS, 2, discretization)

write.table(DS_Dis, file="/Users/miaozhang/Documents/MATLAB/R/DS_Dis.txt",  quote = F, row.names = F, col.names = F)

#make a vector where each element is a line of MATLAB code
#matlab code reads in our variable x, creates two variables y and z, 
#and write z in a csv file
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
#system("/cm/shared/uaapps/matlab/r2017b/bin/matlab -nodisplay -r \"run('/Users/miaozhang/Documents/MATLAB/R/myscript.m'); exit\"")

# read MATLAB output in R
data <- readMat("/Users/miaozhang/Documents/MATLAB/R/DS_Phi.mat")
G_Phi= data$X.mc

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

