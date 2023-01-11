nr.qtl=1000 #no. of SNPs
var_g1<-1
var_g2<-1
nr.trait<-2
rg<-1 #genetic correlation between ancestry
alpha1=0.5 #scale factor for ancestry 1
alpha2=0.5 ##scale factor for ancestry 2


cov_g<-rg*sqrt(var_g1*var_g2)
mu_g<-array(0,nr.trait)
covar<-c(var_g1,cov_g,cov_g,var_g2)
gmat<-matrix(covar,nr.trait,nr.trait)
gmat


plinkbim=read.table("combine2.bim") ##combine2 bim file has common snp between population
v1=seq(1:nrow(plinkbim))
v2=sort(sample(v1,nr.qtl,replace = F))
write.table(plinkbim[v2,2],file='snp_seltemp',col.names=FALSE,row.names=FALSE,append=FALSE,quote=FALSE,sep='\t')


###extract the true eQTL

system(paste('./plink1.9 --bfile combine --extract snp_seltemp --make-bed --out subplinktemp > subplink.log',sep='')) ##master file with common SNP only
system(paste('./plink1.9 --bfile test1 --extract snp_seltemp --make-bed --out test1a'))
system(paste('./plink1.9 --bfile test2 --extract snp_seltemp --make-bed --out test2a'))


###switch the alleles or keep the allele order as in the master (combine plink file)
system(paste('./plink1.9 --bfile test1a --keep-allele-order --a1-allele subplinktemp.bim 5 2 --make-bed --out test1b'))
system(paste('./plink1.9 --bfile test2a --keep-allele-order --a1-allele subplinktemp.bim 5 2 --make-bed --out test2b'))


#get frequency of subset
#system('./mtg2 -plink subplinktemp -frq 1')
system('./mtg2 -plink test1b -frq 1')
system('./mtg2 -plink test2b -frq 1')

##Assign effect of true QTL
#bim<-read.table('subplinktemp.freq',header=FALSE)
bim_1b=read.table('test1b.freq',header=FALSE)
bim_2b=read.table('test2b.freq',header=FALSE)

library(MASS)
gv<-mvrnorm(nr.qtl,t(mu_g),gmat)
v12<-gv*sqrt(1/(bim_1b$V3*(1-bim_1b$V3/2))^(alpha1*2)) 
v13<-gv*sqrt(1/(bim_2b$V3*(1-bim_2b$V3/2))^(alpha2*2))
v3<-cbind(1:length(v2),as.character(bim_1b$V2),v12[,1],bim_1b$V3)
v4<-cbind(1:length(v2),as.character(bim_2b$V2),v13[,2],bim_2b$V3)
write.table(v3, "v3.lst", row.names = F, col.names = F, quote = F)
write.table(v4, "v4.lst", row.names = F, col.names = F, quote = F)

### to get breeding value for each individual
system(paste('./mtg2 -plink test1b -simreal v3.lst'))
system(paste('mv test1b.bv phen1.bv'))

system(paste('./mtg2 -plink test2b -simreal v4.lst'))
system(paste('mv test2b.bv phen2.bv'))




### simulate residual
bv1<-read.table("phen1.bv")
bv2<-read.table("phen2.bv")

##value needs to change according to simulation
vg1<-1
vg2<-1
re<-0
h2_1=0.5
h2_2=0.5
ve1<-vg1*(1-h2_1)/h2_1
ve2<-vg2*(1-h2_2)/h2_2
cov_e<-re*sqrt(ve1*ve2)

library(MASS)
cov.was<-c(ve1,cov_e,cov_e,ve2)
emat<-matrix(cov.was,nr.trait,nr.trait)
mu<-array(0,nr.trait)

ev<-mvrnorm(dim(bv1)[1],t(mu),emat)
out_phen1<-bv1$V1+ev[,1]
out_phen2<-bv2$V1+ev[,2]
fam1=read.table("test1b.fam")
fam2=read.table("test2b.fam")
phen1=cbind(fam1$V1,out_phen1)
phen2=cbind(fam2$V1,out_phen2)

phen=merge(phen1,phen2, by=c("V1"), all=T) #merging phenotype
data=cbind(as.character(phen$V1),as.character(phen$V1),phen$out_phen1, phen$out_phen2)


write.table(data, "phenotype", row.names=F,col.names=F,quote=F)


