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
Genomic relationship matrix (GRM) could be estimated using -rtmx2 and/or -rtmx3 function using mtg2 software. Please see mtg2 manuals (https://www.dropbox.com/s/4fw5vu7seol4xbh/MTG2%20manual.pdf?dl=1) for details.

# Reference 
Momin, M. M., Shin, J., Lee, S., Truong, B., Benyamin, B., & Lee, S. H. (2021). A novel method for an unbiased estimate of cross-ancestry genetic correlation using individual-level data. bioRxiv.

# Contact information
Please contact Hong Lee (hong.lee@unisa.edu.au) or Md Moksedul Momin (cvasu.momin@gmail.com) if you have any queries.
