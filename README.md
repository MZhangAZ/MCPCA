# MCPCA_PopGen
Program to do dimensionality reduction for NGS data

In order to run the R script, users need to install Matlab and download the MCPCA matlab scripts from 
https://github.com/SoheilFeizi/MCPCA.git

Reference: Soheil Feizi and David Tse, Maximally Correlated Principal Component Analysis, arXiv:1702.05471

# Data
## Example data
Attached txt file, DS.txt is an example input data file. Each row repersents a individual and each column repersents a SNP. The dosage value is a continuous variable between 0 and 2, which is the predicted dosage of the non reference allele given the data available.

## Bam file data input
To call the posterior probability of the genotype likelihood from bam file, user can use the program ANGSD

Korneliussen, T. S., Albrechtsen, A., and Nielsen, R. (2014). Angsd: analysis of next generation sequencing data. BMC bioinformatics, 15(1):356.

http://www.popgen.dk/angsd/index.php/Genotype_calling

write the posterior probability of all possible genotypes 

```{angsd}
./angsd -bam bam.filelist -GL 1 -out outfile -doMaf 2 -doMajorMinor 1 -SNP_pval 0.000001 -doGeno 8 -doPost 1 -postCutoff 0.95 
```
gives a output like this:
```{angsd}
chr1  10180 0.844336 0.145660 0.010004 0.865104 0.133812 0.001084 0.072660 0.790897 0.136443 
chr1  69511 0.948844 0.051152 0.000004 0.900312 0.097072 0.002617 0.900312 0.097072 0.002617 
chr1 404664 0.904511 0.093094 0.002395 0.904511 0.093094 0.002395 0.904511 0.093094 0.002395 
chr1 762174 1.000000 0.000000 0.000000 1.000000 0.000000 0.000000 1.000000 0.000000 0.000000 
chr1 762238 1.000000 0.000000 0.000000 0.999851 0.000149 0.000000 0.999995 0.000005 0.000000 
```
then use R script, Pos_DS to calculate dosage value from ANGSD output and gives a output like this:

```{r}
10180    69511    404664   762174   762238
0.165668 0.051160 0.097884      0 0.000000
0.135980 0.102306 0.097884      0 0.000149
1.063783 0.102306 0.097884      0 0.000005
```

# MCPCA R 
If the user doesn't want to call Matlab, there is an testing R script MCPCA_algorithm.R, which can run MCPCA without loading Matlab code. But since it has two for loops, the speed is very slow once the sample size is big
