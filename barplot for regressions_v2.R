# library
library(tidyverse)
library(viridis)
library(ggplot2)
library(grid)
library(gridExtra)



# 1 gene, cor=0.4
plink_res<-read.csv("c:/Users/pthompson/Desktop/test_combos_regression_plink_MC.csv")

plink_long<-gather(test_combos,key="correction",value="power",Bonferroni:BY)

names(plink_long)[6]<-"eff_size"
names(plink_long)[4]<-"Npats"


data<-plink_long

data<-data.frame(Sample_size=data$Npats,Effect_size=data$eff_size,Nsnps=data$nsnp, power=data$power,correction=data$correction)


data$Sample_size<-as.factor(data$Sample_size)
data$Effect_size<-as.factor(data$Effect_size)
data$Nsnps<-as.factor(data$Nsnps)
data$correction<-factor(data$correction,levels = c("Bonferroni","Holm","BH","BY"))


levels(data$Effect_size)<-c("Effect size = 0.1","Effect size = 0.4")



# Make the plot
r1 = ggplot(data,aes(x=Sample_size, y=power, fill=Nsnps,label=power)) +      
  
  # Add the stacked bar
  geom_bar(stat="identity", alpha=0.5) +
  scale_fill_viridis(discrete=TRUE) + facet_grid(Effect_size~correction)+
  geom_text(size = 3, position = position_stack(vjust = 0.5))+
  theme_bw()+theme(axis.title.y = element_blank(),axis.text.y = element_blank(),legend.position = "top",text = element_text(size=14))+xlab("Sample size")

png(filename="c:/Users/pthompson/Dropbox/project SCT analysis/GSCA validation/draft article//plink_reg_MC.png",height=6,width=8,units="in",res=100)  
r1
dev.off()
