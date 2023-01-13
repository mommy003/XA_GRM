# Phenotype simulation
Phenotypes can be simulated following the R script (Simu_phen_combined_pop.r) for combined populations (e.g. South Asian and African) using ancestry-specific allele frequency and ancestry-specific genetic architecture.
- To simulate phenotype in the working directory
```
R CMD BATCH --no-save Simu_phen_combined_pop.r
```
Please ensure that mtg2 and plink1.9 are already installed in the working directory 
- mtg2 software can be downloaded from (https://sites.google.com/view/s-hong-lee-homepage/mtg2)
- plink1.9 software can be downloaded from (https://www.cog-genomics.org/plink/)

# Estimating GRM
Genomic relationship matrix (GRM) can be estimated using -rtmx2 and/or -rtmx3 function using mtg2 software. Please see instruction below or refer MTG2 manual available in https://sites.google.com/view/s-hong-lee-homepage/mtg2 (see section 6.3).

## Genomic relationship matrix for different populations and scale factors
```
./mtg2 -plink toy -rtmx2 rtmx_parameters -out any_file_name
```

### In case of a single population
```
<'rtmx_parameter'> should have
1                                       #number of populations 
anyname.pop                 #file having the pop information (1 and 2)
-0.75                                #scale factor for pop 1 (same as 2 x alpha)
2                                       #when scaling, var(x) or 2p(1-p) can be selected
```
‘anyname.pop’ should have 3 columns (FID, IID, pop_info (1 always)) where FID/IID order should be the same as in .fam file.  
```
FID IID pop_info
fam1 id1 1
fam2 id2 1
fam3 id3 1
……
```
Observed SNP variance or expected SNP variance (2p(1-p)) can be selected. In the last line, 1 = 2p(1-p) and 2 = var(x_i).

### In case of 2 populations
```
<'rtmx_parameter'> should have
2                                       #number of populations 
anyname.pop                 #file having the pop information (1 and 2)
-0.75                                #scale factor for pop 1 (same as 2 x alpha)
-1.5                                  #scale factor for pop 2 (same as 2 x alpha)
2                                       #when scaling var(x) or 2p(1-p) can be selected
```
‘anyname.pop’ should have 3 columns (FID, IID, pop_info (1 or 2)) where FID/IID order should be the same as in .fam file. 

```
FID IID pop_info
fam1 id1 1
fam2 id2 1
fam3 id3 2
……
``` 

Observed SNP variance or expected SNP variance (2p(1-p)) can be selected. In the last line, 1 = 2p(1-p) and 2 = var(x_i).

### In case of N populations
```
<'rtmx_parameter'> should have
N                                       #number of populations 
anyname.pop                 #file having the pop information (1 – N)
-0.75                                #scale factor for pop 1 (same as 2 x alpha)
-1.5                                  #scale factor for pop 2 (same as 2 x alpha)
-1.0                                  #scale factor for pop N (same as 2 x alpha)
2                                       #when scaling var(x) or 2p(1-p) can be selected
```
‘anyname.pop’ should have 3 columns (FID, IID, pop_info (1 – N)) where FID/IID order should be the same as in .fam file. 

```
FID IID pop_info
fam1 id1 1
fam2 id2 1
fam3 id3 2
fam4 id4 4
fam5 id5 10
……
```

Observed SNP variance or expected SNP variance (2p(1-p)) can be selected. In the last line, 1 = 2p(1-p) and 2 = var(x_i).
When there are many missing genotypes, e.g. a substantial proportion of SNPs is totally missing in one ancestry, -rtmx3 can be used with the same input files.

```
./mtg2 -plink toy -rtmx3 rtmx_parameters -out any_file_name
```

-rtmx3 accounts for missing genotypes when estimating genomic relationships between two ancestries.
NOTE. -thread n command can parallel the computation with n CPUs  
Please see example4 to be downloaded in the web. 
Please download the latest version from the web if you have one before Aug/22.
https://sites.google.com/view/s-hong-lee-homepage/mtg2

# Reference 
Momin, M. M., Shin, J., Lee, S., Truong, B., Benyamin, B., & Lee, S. H. (2021). A novel method for an unbiased estimate of cross-ancestry genetic correlation using individual-level data. bioRxiv.

# Contact information
Please contact Hong Lee (hong.lee@unisa.edu.au) or Md Moksedul Momin (cvasu.momin@gmail.com) if you have any queries.
