
combos<-expand.grid(nsnpeff=c(10,5,0),nsnp=c(20,40,100),ngenes=c(1,2,4),nsub=c(100,500),cluster=c(0,1),gpcorr=c(0.1,0.4))

for(cb in 1:length(combos[,1]))
{
myfitsummaryPT<-read.csv(paste0("h:/DVMB/Genetic_analysis_SCT_twin/simulations_results/simulation_",cb,"_results.csv")) 


if(combos[cb,3]==1){gene_power<-length(which(myfitsummaryPT[,2]<.05))} else{gene_power<-apply(myfitsummaryPT[,(combos[cb,3]+1):(combos[cb,3]*2)],2,function(x2) length(which(x2<(.05/combos[cb,3]))))}

gene_power2<-rep(NA,4)
gene_power2[1:length(gene_power)]<-gene_power

model_out_summary2<-c(combos[cb,2],combos[cb,1],combos[cb,3],combos[cb,6],gene_power2,100,combos[cb,5],combos[cb,4])


write.table(t(model_out_summary2), file = paste0("h:/DVMB/Genetic_analysis_SCT_twin/simulations_results2/PT_simulation_results_MC_adj.csv"), sep = ",", append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE)

}