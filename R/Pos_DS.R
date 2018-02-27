
########################################################
# read posterior prob from ANGSD output in R and calculate dosage
#
########################################################
dat= read.table(gzfile("/Users/miaozhang/Desktop/Pos_ALL.geno.gz")) 

snpName= dat[,2]
dat= dat[,-c(1,2)]

n= 42
p= dim(dat)[1]
ds= matrix(0, p, n)

for (i in 1:p){
  for(j in 1:n){
    ds[i,j]= dat[i,3*(j-1)+2] + dat[i,3*(j-1)+3]*2
  }
}

ds= t(ds)

colnames(ds)= snpName