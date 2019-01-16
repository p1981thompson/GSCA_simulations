#simulate genotype-phenotype data

# 21st Sep 2018
# no linkage disequilibrium ie no correl between adjacent SNPS
#- NB not sure if previously this was only specified in clustered condition?
# 
# But we will ignore it and assume researchers select SNPs that are not in LD

# cycle by *proportion* of SNPs with effect, rather than abs N (though retain old nsnp variable for total snps in analysis)
# Altered how correlation with phenotype is computed
#  Correlation in mvrnorm is with single latent phenotype.
# Then create indicators of the phenotype, which are correlated with it.
# Doing it this way makes it less likely to get matrix which is not positive definite - esp if there is LD

require(doBy)
library(tidyr)
require(tidyverse)
require(MASS)
require(stats)
library(doSNOW)
library(foreach)
library(ASGSCA)
library(matrixcalc)
library(tictoc)

options(scipen = 999) #turn off scientific notation
nrun<-100
startrow<-1 #row of big summary to write to
mynperm=100 #for GSCA - probably need more than this but to test prog use 100 for speed

npop<-10000 #size of population to simulate
n2sim<-3 #n phenos indexing latent pheno factor : nb SNPS correlate withlatent factor

LDbase <- 0 #This is the highest correlation between adjacent SNPs - correl will decline with distance
#original file had max LD of .9. Later I tried with reduced value of .75

phencorr<-.75 #correlation of 3 phenotypes with latent variable
#with this setting, intercorrelation between phenotypes is around .56

effsizelist<-c(.1,.2)
psnpefflist<-c(.2,.4,.8)
ngenelist<-c(1,2)
nsublist<-c(100,200,400)
nsnp_per_gene<-c(10,20)

bigsumname<-'bigsum_nsnp1020_eff12_peff24_ngen12_nsub_12_perm100_run100_combined' #name for this big summary file
bigsumNrow<-nrun*length(effsizelist)*length(psnpefflist)*length(nsublist)*sum(ngenelist)*length(nsnp_per_gene)
#NB we have to do 'sum' for ngenelist, as there is a row for each gene
bigsummary<-data.frame(matrix(NA,ncol=13,nrow=bigsumNrow))
colnames(bigsummary)<-c('run','cluster','nsnp','ngenes','psnpeff','nsnpeff','effsize','nsub','condition','gene','path','p','sig.p')

cluster<-0  #previously in a loop - now redundant: only looking at cases where equal effect in all genes for now
 mycondition<-0 #initialise
 tic() #start timer
 for (nsnpg in nsnp_per_gene){
    for (effsize in effsizelist){
      for(psnpeff in psnpefflist){ #cycle by *proportion* of SNPs with effect rather than N SNPs
        for (ngenes in ngenelist){
          nsnp<-nsnpg*ngenes #compatibility with previous version where nsnp was all in analysis summed over genes
          nsnpeff<-psnpeff*nsnp #N snps with an effect
          print(paste('cluster',cluster,'nsnp',nsnp,'ngenes',ngenes,'effsize',effsize,'psnpeff',psnpeff,'nsnpeff',nsnpeff))
          myskip<-ngenes*cluster #no difference between cluster/no cluster if just one gene, so skip cluster version
          if(myskip!=1){
            gene_break<-floor(nsnp/ngenes)
            gpcov<-effsize#effect size of SNPs with effect, set as correlation 
            snpnames<-paste0('snp', 1:nsnp)
            maf<-runif(nsnp,.25,.5) #minor allele freq set to random value from .25 to .5
            
            #Setup a dataframe to hold to hold simulated data.
            mydata<-data.frame(matrix(nrow=npop,ncol=(4+sum(nsnp)))) #set to hold npop cases - we will sample from these for each run
            #4 additional columns with hold (a)latent phenotype and (b) 3 component variables
            #SNP values are 0, 1 or 2, but we start with random normal numbers
            
            #-----------------------------------------------------------------
            #simulate the SNPS as correlated zscores, correlation is mycorr
            #-----------------------------------------------------------------
            #Where there is correlation, assume correlation depends on how close SNPs are, 
            #achieved by making size of correlation proportional to distance
            # (where sequence order is proxy for distance)
            # NB we are going to ignore linkage disequilibrium for now, so just have mycorr as zero
            
            h <-nsnp+1
            mymean<-rep(0,h) #vector of means to use in mvrnorm function
            mycorr<-0 #default is uncorrelated
            mycov2<-matrix(mycorr,nrow=h,ncol=h) #cov matrix for mvrnorm
            
            diag(mycov2)<-rep(1,h) #ones on diagonal of cov matrix
     #Now specify correlation between SNPs, i.e. linkage disequilibrium       
            if(LDbase>0){ #need to compute correlation for conditions where SNP effects clustered in genes
              for (i in 1:h){
                irange<-1+as.integer((i-1)/gene_break) #irange specifies haplotype block for i
                for (j in 1:h){
                  jrange<- 1+as.integer((j-1)/gene_break) #jrange specifies haplotype block for j
                  if(irange==jrange){
                    k<- abs(i-j)
                    thisr<-ifelse(nsnp<=20,LDbase-gene_break*k/100,(LDbase-gene_break*k/((nsnp^2)/(10*ngenes)))) #tweaked so magnitude of correlation declines with distance between SNPs (PT: further tweek to allow different nsnp without problems.)
                    if(thisr<0){thisr<-0}
                    mycov2[i,j]<-thisr
                    mycov2[j,i]<-mycov2[i,j] #correlation determined by distance between SNPs!
                    if(k==0){mycov2[i,j]<-1}
                  }
                }
              }
            }
            
            #Identify which SNPs are correlated with latent phenotype
            #Here assume effects divided equally betewen genes
            
              nsnp_g<-round(nsnp/ngenes,0)
              nsnpeff_g<-round(nsnpeff/ngenes,0)
              snpeff_pos<-vector() #initialise vector for holding indices of SNPs with effect
                for (g in 1:ngenes){
                  myrange<-(1+(g-1)*nsnp_g):(g*nsnp_g) #range of indices for SNPs in this gene
                  snpeff_pos<-c(snpeff_pos,sample(myrange,nsnpeff_g,replace=F)) #randomly selected SNPs in this range
                }
            mycov2[snpeff_pos,h]<-effsize #for SNPs with effect, add effect size in cov matrix
            mycov2[h,snpeff_pos]<-effsize
            
                        if(matrixcalc::is.positive.definite(mycov2)==FALSE)
{mycov2<-Matrix::nearPD(mycov2,keepDiag=TRUE)$mat}
            
            #create simulated data for SNPs +latent phenotype only
            mymean<-rep(0,h)
            mydata=data.frame(mvrnorm(n = npop, mymean, mycov2))
            colnames(mydata)[h]<-'latent'
            #################################################
      
            #values for phenotype created as correlated with latent: r = phencorr
            mydata$nonword<-phencorr*mydata$latent+sqrt(1-phencorr^2)*rnorm(nrow(mydata))
            mydata$elang<-phencorr*mydata$latent+sqrt(1-phencorr^2)*rnorm(nrow(mydata))
            mydata$rlang<-phencorr*mydata$latent+sqrt(1-phencorr^2)*rnorm(nrow(mydata))
            
            #################################################
            
            colnames(mydata)[1:nsnp]<-snpnames
            
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
            plotter=FALSE
            if(plotter==TRUE){
              library(tidyr)
              
              cor.data<-as.matrix(abs(myr))
              cor.data[lower.tri(cor.data)] <- NA
              
              cor_tri_N <- as.data.frame(cor.data) %>% 
                mutate(Var1 = factor(row.names(.), levels=row.names(.))) %>% 
                gather(key = Var2, value = value, -Var1, na.rm = TRUE, factor_key = TRUE) 
              
              
              ggplot(data = cor_tri_N, aes(Var2, Var1, fill = value)) + 
                geom_tile()+scale_fill_gradient(low = "yellow", high = "red")+theme_bw()+ theme(axis.text.x = element_text(angle = 90, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.title = element_blank())
              
              #  ggsave(paste0("h:/DVMB/Genetic_analysis_SCT_twin/simulations_results/plots/plot_",cb,".pdf"))
              
            }
            
            
            #----------------------------------------------------------------------
            # For each of nrun runs take a sample and analyse
            #----------------------------------------------------------------------
            for (nsub in nsublist){
              mycondition<-mycondition+1
            for (myrun in 1:nrun){
            
                
                myrows<-sample(npop,nsub) 
                mysample<-mydata[myrows,-h] #skip latent phenotype
                
                #use ASCSCA to analyse
                snprange<-1:nsnp
                phenrange<-((nsnp+1):(nsnp+3))
                ObservedVar=colnames(mysample)[c(snprange,phenrange)] #skip latent phenotype - not used in analysis
          
                LatentVar=c(paste0("gene",1:ngenes),"Neurodev")
                
                #W0 is I x L matrix where rows are genotypes and traits, columns are genes and latent phenotype
                W0=matrix(rep(0,length(LatentVar)*(n2sim+nsnp)),nrow=n2sim+nsnp,ncol=length(LatentVar), dimnames=list(ObservedVar,LatentVar))
                
                if(ngenes==1){W0[1:nsnp,1] = 1}
                
                ind<-1:ngenes*gene_break
                #placer=ifelse(length(ind)<2,NA,c(1:ind[1],seq((ind[i-1]+1):ind[2],(ind[2]+1):ind[3],(ind[3]+1):ind[4]))
                
                for(w in 1:ngenes)
                {
                  if(w==1)
                  {W0[1:ind[1],w] = 1}
                  if(w>1){W0[(ind[w-1]+1):ind[w],w] = 1}
                }
                W0[(nsnp+1):(nsnp+3),(ngenes+1)]=1
                
                
                #B0 is L x L matrix with 0s and 1s, where 1 indicates arrow from latent geno in row to latent phenotype in column
                B0=matrix(rep(0,(ngenes+1)*(ngenes+1)),nrow=(ngenes+1),ncol=(ngenes+1), dimnames=list(LatentVar,LatentVar))
                B0[1:ngenes,(ngenes+1)]=1
                
                # GSCA(mysample,W0, B0,latent.names=LatentVar, estim=TRUE,path.test=FALSE,path=NULL,nperm=1000)
                # for quick scrutiny of weights use this -but for pvalues need slow version using path.test
                
                myfit<-GSCA(mysample,W0, B0,latent.names=LatentVar, estim=TRUE,path.test=TRUE,path=NULL,nperm=mynperm)
                
                runsummary<-c(myfit$Path[1:ngenes,ngenes+1],myfit$pvalues[1:ngenes,ngenes+1])
                runsummary<-data.frame(runsummary)
                endrow<-startrow+ngenes-1
                
                bigsummary[startrow,1:9]<-c(myrun,cluster,nsnp,ngenes,psnpeff,nsnpeff,gpcov,nsub,mycondition)
                bigsummary$gene[startrow:endrow]<-1:ngenes
                bigsummary$path[startrow:endrow]<-runsummary[1:ngenes,]
                bigsummary$p[startrow:endrow]<-runsummary[(ngenes+1):(2*ngenes),]
                startrow<-endrow+1 #update startrow for next run
              }
            }
            rm(mydata)
          }
        }
      }
    }
 }
 t2<-toc()
bigsummary$sig.p<-0
w<-which(bigsummary$p<.05)
bigsummary$sig.p[w]<-1
bigsummary<-bigsummary[1:endrow,]
#fill in blank rows
for (r in 1:endrow){
if(is.na(bigsummary$run[r]))
  {bigsummary[r,1:9]<-bigsummary[(r-1),1:9]}
}

saveRDS(bigsummary,bigsumname)
bigsummary$sig.p.cor<-0
w<-which(bigsummary$p<.025 & bigsummary$ngene==2)
bigsummary$sig.p.cor[w]<-1
w2<-which(bigsummary$p<.05 & bigsummary$ngene==1)
bigsummary$sig.p.cor[w2]<-1

by_N_neff_nsub <- bigsummary %>% group_by(nsub,ngenes,gene,nsnp,psnpeff,effsize,condition)


groupedsummary_new<-by_N_neff_nsub %>% summarise(power = mean(sig.p.cor))

#groupedsummary_new$cond2<-c(1,2,1,2,3,4,3,4,5,6,5,6,7,8,7,8,9,10,9,10,
#                        11,12,11,12,13,14,13,14,15,16,15,16,17,18,17,18,
#                        19,20,19,20,21,22,21,22,23,24,23,24)

groupedsummary_new$psnpeffplot<-groupedsummary_new$psnpeff+groupedsummary_new$ngenes/100-.02

groupedsummary_new$nsnp<-as.factor(groupedsummary_new$nsnp)
groupedsummary_new$ngenes<-as.factor(groupedsummary_new$ngenes)
groupedsummary_new$nsub<-as.factor(groupedsummary_new$nsub)
groupedsummary_new$gene<-as.factor(groupedsummary_new$gene)

library(ggplot2)

ggplot(groupedsummary_new,aes(x=psnpeffplot,y=power,colour=ngenes,shape=nsub))+geom_point()+geom_line(aes(linetype=gene))+facet_grid(nsnp~effsize)+theme_bw()+theme(legend.position = "bottom")
       

 # saveRDS(groupedsummary,paste0('grp_',bigsumname)) #didnt work?
  write_csv(groupedsummary_new,paste0('grp_',bigsumname,'.csv'))
  write_csv(bigsummary,paste0(bigsumname,'.csv'))
 # 
