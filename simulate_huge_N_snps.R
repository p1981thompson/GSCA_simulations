#####################################################################################
#
#   Adapted script based on Dorothy's simulation code.
#
#####################################################################################

#12-02-2018

#PT fixed error with positive definite matrix with large number of SNPs
#PT fixed error with categorisation of counts with small p in qnorm.

require(doBy)
require(tidyverse)
require(MASS)
require(stats)
options(scipen = 999) #turn off scientific notation
ngene<-4
genenames<-c('CNTNAP2','ATP2C2','FOXP2','UBR1')
nsnp<-c(4286,613,401,143) #n SNPs per gene (currently all same but option to vary this)

ncases<-135
mydata<-data.frame(matrix(nrow=ncases,ncol=(5+4+sum(nsnp))))
#We will simulate data for 3 phenotypes and 38 SNPs for each of 4 genes, plus 4 PCAs
#For the NLGN genes, values can be 0, 1, 2, 3 for trisomies (twins not yet included in simulation)
#For the other two autosomal genes, values are 0, 1 or 2.
colnames(mydata)[1:5]<-c('Karyotype','Ethnicity','NwdRepPheno','LangPheno','NeurodevPheno')
 
 #start by simulating Karyotype: assume 40 of each (irrelevant to analysis at present)  
mydata$Karyotype<-as.factor(c(rep(1,45),rep(2,45),rep(3,45)))
mydata$Ethnicity<-as.factor(c(rep(1,105),rep(2,15),rep(3,15)))
levels(mydata$Karyotype)<-c('XXX','XXY','XYY')
levels(mydata$Ethnicity)<-c('White','Asian','Black')#assume 100,20,20
#-----------------------------------------------------------------

#Simulate the 3 phenotypes as correlated zscores, correlation is mycorr
mymean<-rep(0,3)
mycorr<-.75
mycov<-matrix(mycorr,nrow=3,ncol=3)
diag(mycov)<-rep(1,3)
mydata[,3:5]=mvrnorm(n = ncases, mymean, mycov)
#-----------------------------------------------------------------
#3rd phenotype (column 5) will be ordinal with values 0, 1, 2.
#Assume 25% in categories 0 and 2; corresponds to zscore of +/- .67
temp<-mydata[,5]
mydata[,5]<-2 #default
w<-which(temp< .67)
mydata[w,5]<-1
w<-which(temp< -.67)
mydata[w,5]<-0
#-----------------------------------------------------------------
#Now for each gene, simulate the 38 SNPS as correlated zscores, correlation is mycorr
#We will assume correlation depends on how close SNPs are, ranging from .9 for adjacent and stepping 
#down by .1 for each additional distance (where sequence order is proxy for distance)

startcol<-6
for (g in 1 :ngene){
  if (g>1){ startcol<-startcol+nsnp[g-1]}
  endcol<-startcol+nsnp[g]-1
  h <-nsnp[g]
    mymean<-rep(0,h)
   mycorr<-0 #default is uncorrelated
   mycov<-matrix(mycorr,nrow=h,ncol=h)
   diag(mycov)<-rep(1,h)
   for (i in 1:h){
     for (j in 1:h){
       k<- abs(i-j)
       mycov[i,j]<-(1-2*k/100) #correlation determined by distance between SNPs!
       
  
       if(k==0){mycov[i,j]<-1}
     }
   }
   mycov <- (mycov * lower.tri(mycov)) + t(mycov * lower.tri(mycov)) 
diag(mycov) <- 1 
mycov<-Matrix::nearPD(mycov)$mat
   
      mydata[,startcol:endcol]=mvrnorm(n = ncases, mymean, mycov)
   colnames(mydata)[startcol:endcol]<-paste0(genenames[g],'_',seq(1:h))
}


#-------------------------------------------------------------------------
#Convert gene scores to integers: 0-3 for NLGNs and to 0-2 for autosomal

firstcol<-6
lastcol<-length(mydata)
p<-c(0,0,0,0) #initialise a vector to hold p values for different N alleles
autosomestart<-which(colnames(mydata)=='CNTNAP2_1')
for (g in firstcol:lastcol){
  maf<-runif(1,min=.2,max=.5) #Let MAF for each SNP range from .2 to .5
  
  #Now find probability of given N minor alleles
  #NB because R doesn't allow zero indices, p[1] is probability of zero, p[2] is probability of 1 and so on

  if (g<autosomestart){ #triploid genes
  p[1]<-(1-maf)^3
  p[4]<-maf^3
  p[2]<-3*(1-maf)^2*maf
  p[3]<-3*(maf^2*(1-maf))
  }
  
  if (g>(autosomestart-1)){ #diploid genes
    p[1]<-(1-maf)^2
    p[4]<-0.00000000001
    p[2]<-2*(1-maf)*maf
    p[3]<-maf^2
  }
  #now use p-values to assign z-score cutoffs that convert to 0,1,2 or 3 minor alleles
  # temp<-mydata[,g]
  # w<-which(temp<qnorm(p[1]))
  # mydata[w,g]<-0
  # w<-which(temp>qnorm(p[1]))
  # mydata[w,g]<-1
  # w<-which(temp>qnorm(p[2]+p[1]))
  # mydata[w,g]<-2
  # w<-which(temp>qnorm(p[1]+p[2]+p[3]))
  # mydata[w,g]<-3
  
  mydata[,g]<-car::recode(qcat(pnorm(mydata[,g]),prob=p),"1=0;2=1;3=2;4=3")
}
#----------------------------------------------------------------------------------------
#Explore correlation structure of genes with PCA
#Additional columns at end of mydata hold PCA scores based on all SNPs in that gene
for (mygene in 1:ngene){
  thiscol<-5+sum(nsnp)+mygene
  colnames(mydata)[thiscol]<-paste0('PCA_',genenames[mygene])
  startcol<-6+(mygene-1)*nsnp[mygene]
  endcol<-6+(mygene*nsnp[mygene])-1
  mypca<-prcomp(mydata[,startcol:endcol])
  plot(mypca)
  myloadings <- as.matrix(mypca$rotation[,1]) #loadings for 1st component
  for (i in 1:ncases){
    thisrow<-as.matrix(t(mydata[i,startcol:endcol]))#transpose for matrix crossproduct
    mydata[i,thiscol]<-crossprod(myloadings,thisrow)
  }

#If looking at single SNP can modify this bit to print table of means and SDs by N alleles
# mytable<-summaryBy(LangPheno ~ N.alleles, data = mydata,
#                     FUN = function(x) { c(N=length(x),m = mean(x), s = sd(x)) } )
# mytable
myvar<-colnames(mydata)[thiscol]
mymodel<-lm(LangPheno~get(myvar,mydata),data=mydata)
#check this is right variable : mymodel1<-lm(LangPheno~PCA_NLGN3,data=mydata)
summary(mymodel)
mymodel$coefficients
}
#write.table(mycov, "trygenevoc.txt", sep=",",row.names=FALSE) 
#write.table(mydata, "dummydata.txt", sep=",",row.names=FALSE) 