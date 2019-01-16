#simulate genotype-phenotype data

# 10th Sept 2018 - Paul T

#http://www.stat.purdue.edu/bigtap/online/docs/simple-association-test.html
# Looks like multiple single regressions are fitted, then the p-value for parameters are extracted.
#Covariate adjustment is possible, but we do not use this in theis analysis

# simulate data for analysis with N SNPs with some linkage disequibrium, and 3 correlated phenotyes
# simulate 3 situations
# 1. no relationship between geno and pheno
# 2. strong relationship (r = .4) for small N SNPs with pheno
# 3. weak relationship (r = .1) for large N SNPs with pheno

#Aim: start by using GSCA approach and see if it captures the associations

require(doBy)
library(tidyr)
require(tidyverse)
require(MASS)
require(stats)
library(doSNOW)
library(foreach)
library(ASGSCA)
library(matrixcalc)
#set.seed(1981)
set.seed(Sys.time(),intern=TRUE)
options(scipen = 999) #turn off scientific notation

gsca_sim_multigene_plink<-function(nsnp,nsnpeff,ngenes=1,nsub=120,ncases=100000,gpcorr=0.4,n2sim=3,plotter=TRUE,cb=cb,cluster=cluster,correct=c('bonferroni','holm','BH','BY'))
{
  gene_break<-floor(nsnp/ngenes)
  
  #set the correlation 
  gpcov<-gpcorr
  
  snpnames<-paste0('snp', 1:nsnp)
  maf<-runif(nsnp,.25,.5) #minor allele freq set to value from .25 to .5
  
  #Setup a dataframe to hold to population simulated data.
  mydata<-data.frame(matrix(nrow=ncases,ncol=(3+sum(nsnp))))
  
  
  myfilenames <- c('DatasetRand_N',
                   paste0('DatasetUncorr_',nsnpeff,'SNP_',10*gpcorr),
                   'Dataset4Block_N',
                   paste0('Dataset',ngenes,'Block_',nsnpeff,'SNP_',10*gpcorr))
  thisfile<-4 #User specified: Indicates which type of correlation structure in datafile (could use loop here)
  mydatafile<- myfilenames[thisfile]
  
  #SNP values are 0, 1 or 2, but we start with random normal numbers
  
  #-----------------------------------------------------------------
  #simulate the SNPS as correlated zscores, correlation is mycorr
  #-----------------------------------------------------------------
  #Where there is correlation, assume correlation depends on how close SNPs are, 
  #achieved by making size of correlation proportional to distance
  # (where sequence order is proxy for distance)
  
  h <-nsnp
  
  mymean<-rep(0,h) #vector of means to use in mvrnorm function
  mycorr<-0 #default is uncorrelated
  mycov2<-matrix(mycorr,nrow=h,ncol=h) #cov matrix for mvrnorm
  
  diag(mycov2)<-rep(1,h) #ones on diagonal of cov matrix
  if(thisfile>2){ #only add correlation for conditions where thisfile is 3 or 4
    for (i in 1:h){
      irange<-1+as.integer((i-1)/gene_break) #irange specifies haplotype block for i
      for (j in 1:h){
        jrange<- 1+as.integer((j-1)/gene_break) #jrange specifies haplotype block for j
        if(irange==jrange){
          k<- abs(i-j)
          thisr<-ifelse(nsnp<=20,1-gene_break*k/100,(1-gene_break*k/((nsnp^2)/(10*ngenes)))) #tweaked so magnitude of correlation declines with distance between SNPs (PT: further tweek to allow different nsnp without problems.)
          if(thisr<0){thisr<-0}
          mycov2[i,j]<-thisr
          mycov2[j,i]<-mycov2[i,j] #correlation determined by distance between SNPs!
          if(k==0){mycov2[i,j]<-1}
        }
      }
    }
  }
  
  #################################################
  
  mycov_PT<-matrix(mycorr,nrow=nsnp+n2sim,ncol=nsnp+n2sim)
  mycov_PT[1:nsnp,1:nsnp]<-mycov2
  m<-diag(n2sim)
  pheno_mat<-m[lower.tri(m)|upper.tri(m)]<-0.75
  mycov_PT[(nsnp+1):(nsnp+n2sim),(nsnp+1):(nsnp+n2sim)]<-pheno_mat
  
  #snpeff_pos<-sample(1:nsnp,nsnpeff,replace=F)
  
  if(!nsnpeff==0){
    if(cluster==0){snpeff_pos<-sample(1:nsnp,nsnpeff,replace=F)} else {snpeff_pos<-c((nsnp-(nsnpeff-1)):nsnp)}
    
    mycov_PT[snpeff_pos,(nsnp+1):(nsnp+n2sim)] <- runif(nsnpeff*n2sim,gpcov-0.05,gpcov+0.05)
    mycov_PT[(nsnp+1):(nsnp+n2sim),snpeff_pos] <- t(mycov_PT[snpeff_pos,(nsnp+1):(nsnp+n2sim)])
  }
  
  mycov<-mycov_PT
  diag(mycov)<-rep(1,dim(mycov)[1])
  #################################################
  
  #then set diagonal values to 1 for all
  #diag(mycov)<-rep(1,n2sim)
  if(matrixcalc::is.positive.definite(mycov)==FALSE)
  {mycov<-Matrix::nearPD(mycov,keepDiag=TRUE)$mat}
  
  
  
  mymean<-rep(0,dim(mycov)[1])
  mydata=mvrnorm(n = ncases, mymean, mycov)
  
  mydata<-as.data.frame(mydata)
  
  colnames(mydata)[1:nsnp]<-snpnames
  colnames(mydata)[(nsnp+1):(nsnp+3)]<-c('NwdRepPheno','LangPheno','NeurodevPheno')
  
  
  
  #-------------------------------------------------------------------------
  #Convert gene scores to integers: 0-2 for autosomal
  
  firstcol<-1
  lastcol<-nsnp
  p<-c(0,0,0) #initialise a vector to hold p values for different N alleles
  for (i in 1:nsnp){
    p[1]<-(1-maf[i])^2
    p[2]<-2*(1-maf[i])*maf[i]
    p[3]<-maf[i]^2
    
    #now use p-values to assign z-score cutoffs that convert to 0,1,2 or 3 minor alleles
    temp<-mydata[,i]
    w<-which(temp<qnorm(p[1]))
    mydata[w,i]<-0
    w<-which(temp>qnorm(p[1]))
    mydata[w,i]<-1
    w<-which(temp>qnorm(p[2]+p[1]))
    mydata[w,i]<-2
    
  }
  
  myr<-cor(mydata) #to check you have desired correlation structure, View(myr)
  
  
  #----------------------------------------------------------------------
  if(plotter==TRUE){
    library(tidyr)
    
    cor.data<-as.matrix(abs(myr))
    cor.data[lower.tri(cor.data)] <- NA
    
    cor_tri_N <- as.data.frame(cor.data) %>% 
      mutate(Var1 = factor(row.names(.), levels=row.names(.))) %>% 
      gather(key = Var2, value = value, -Var1, na.rm = TRUE, factor_key = TRUE) 
    
    
    ggplot(data = cor_tri_N, aes(Var2, Var1, fill = value)) + 
      geom_tile()+scale_fill_gradient(low = "yellow", high = "red")+theme_bw()+ theme(axis.text.x = element_text(angle = 90, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.title = element_blank())
    
    #ggsave(paste0("h:/DVMB/Genetic_analysis_SCT_twin/simulations_results/plots/plot_PLINK_",cb,".pdf"))
    ggsave(paste0("c:/Users/pthompson/Desktop/plot_PLINK_",cb,".pdf"))
    
  }
  
  
  
  #----------------------------------------------------------------------
  # For each of nrun runs take a sample and analyse
  #----------------------------------------------------------------------
  nrun <- 100 #N runs to simulate
  
  mybigdata<-mydata
  
  
  
  ###########################################################################
  # loop over iterations for power
  ###########################################################################
  
  
  #  pb <- txtProgressBar(max=100, style=3)
  #  progress <- function(n) setTxtProgressBar(pb, n)
  #  opts <- list(progress=progress)
  
  SNP_p<-array(0,c(nsnp,nrun,3)) #array to hold the 3 matrices (Nsnps [rows] x Niterations [cols]) 
  for(myn in 1:nrun)
  {
    myrows<-sample(ncases,nsub) #sample of row numbers for sampling equal to N=nsub
    mysample<-mybigdata[myrows,] #take sample of size Nsub from large simulated data.
    
    #PLINK - type regressions
    
    x<-mysample[,1:nsnp] #extract the data for the SNPs as predictors variables.
    y<-mysample[,(nsnp+1):(nsnp+3)] #extract the data for the phenotypes as dependent variable.
  
    #loop for each phenotype
    for(i in 1:3)
    {
      SNP_p[,myn,i]<-t(apply(x, 2, function(x.col) summary(lm(y[,i]~x.col))$coefficients[2,4])) #extract p-value for each SNP parameter.
    }
    
  }
  #close(pb)
  gene_power<-data.frame(y1=rep(0,nsnp),y2=rep(0,nsnp),y3=rep(0,nsnp),y4=rep(0,nsnp)) #dataframe to hold the extracted number of significant runs from the 100, to calc power estimate.
  #alpha_B<-0.05/(nsnp) #Bonferrroni corrected alpha.
  alpha<-0.05
  #loop over phenotypes
  for(i in 1:4)
  {
    gene_power[,i]<-apply(SNP_p[,,1],1,function(x) length(which(p.adjust(x,method=correct[i],n=length(x))<alpha))) #extract the N sig. runs.
  }
  
  snpeff_power<-gene_power[snpeff_pos,]
  
  colnames(snpeff_power)<-colnames(gene_power)<-c('Bonferroni','Holm','BH','BY')
  
  return(list(gene_power=gene_power,SNP_p=SNP_p,snpeff_pos=snpeff_pos,snpeff_power=snpeff_power))
}


#----------------------------------------------------------------------------
# User specifies how many SNPs have effect on phenotype and how big an effect here
# These values are ignored if thisfile is 1 or 3 (no effect)
#nsnpeff<-5 #N snps exerting effect on phenotype
#gpcorr<-gpcov<-.4 # effect size correlation
#n2sim<-3 #N phenotypes
#ngenes<-1 # N genes
#nsub<-120 #N subjects - set to resemble our study
#----------------------------------------------------------------------------


#combos_plink<-c(nsnpeff=2,nsnp=50,ngenes=1,nsub=500,cluster=0,gpcorr=0.4,cb=3) #Arguments supplied to the function.


combos_plink<-expand.grid(nsnpeff=1,nsnp=c(20,40,100),ngenes=1,nsub=c(100,500),cluster=0,gpcorr=c(0.1,0.4))

test_combos<-data.frame(combos_plink,Bonferroni=rep(NA,length(combos_plink[,1])),Holm=rep(NA,length(combos_plink[,1])),BH=rep(NA,length(combos_plink[,1])),BY=rep(NA,length(combos_plink[,1])))

for(cb in 1:length(combos_plink[,1]))
{
  print(cb)
test_combos[cb,7:10]<-gsca_sim_multigene_plink(nsnp=combos_plink[cb,2],nsnpeff=combos_plink[cb,1],ngenes=combos_plink[cb,3],nsub=combos_plink[cb,4],ncases=50000,gpcorr=combos_plink[cb,6],n2sim=3,plotter = TRUE,cb=cb,cluster=combos_plink[cb,5],correct=c('bonferroni','holm','BH','BY'))$snpeff_power
}

save.image()
savehistory()

write.csv(test_combos,"c:/Users/pthompson/Desktop/test_combos_regression_plink_MC.csv",row.names=FALSE)
#------------------------------------------------------------------------
# If we manipulate the nsub, nsnp, or gpcorr, then the power changes dramatically.
#------------------------------------------------------------------------


