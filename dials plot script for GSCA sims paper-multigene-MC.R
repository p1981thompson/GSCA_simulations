# library
library(tidyverse)
library(viridis)
library(ggplot2)
library(grid)
library(gridExtra)



# 1 gene, cor=0.4
plot.data<-read.csv("H:/DVMB/Genetic_analysis_SCT_twin/simulations_results2/PT_simulation_results_MC_adj_formatted.csv")

data<-plot.data[plot.data$ngenes==1 &plot.data$eff_size==0.4&plot.data$cluster==0,]

data<-data.frame(individual.lab=data$Npats,group=data$nseff,observation=data$nsnp, value=data$nsig1)

data$individual<-paste0(data$individual.lab,":",data$group,":",data$observation)
# Transform data in a tidy format (long format)
#data = data %>% gather(key = "observation", value="value", -c(1,2))

data$group<-as.factor(data$group)
data$observation<-as.factor(data$observation)
data$individual<-as.factor(data$individual)
data$group<-as.factor(data$group)

# Set a number of 'empty bar' to add at the end of each group
empty_bar=2
nObsType=nlevels(as.factor(data$observation))
to_add = data.frame( matrix(NA, empty_bar*nlevels(data$group)*nObsType, ncol(data)) )
colnames(to_add) = colnames(data)
to_add$group=rep(levels(data$group), each=empty_bar*nObsType )
data=rbind(data, to_add)
data=data %>% arrange(group, individual)
data$id=rep( seq(1, nrow(data)/nObsType) , each=nObsType)

# Get the name and the y position of each label
label_data= data %>% group_by(id, individual) %>% summarize(tot=sum(value))
number_of_bar=nrow(label_data)
angle= 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
#label_data$hjust<-ifelse( angle < -90, 1, 0)
label_data$hjust <-1
label_data$angle<-ifelse(angle < -90, angle+180, angle)

label_data$individual <- substring(label_data$individual,1,3)


# prepare a data frame for base lines
base_data=data %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data = base_data
grid_data$end = grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start = grid_data$start - 1
grid_data=grid_data[-1,]

# Make the plot
p1 = ggplot(data) +      
  
  # Add the stacked bar
  geom_bar(aes(x=as.factor(id), y=value, fill=observation), stat="identity", alpha=0.5) +
  geom_text(aes(x=as.factor(id), y=value, fill=observation,label=value),size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_viridis(discrete=TRUE) +
  
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 100, xend = start, yend = 100), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 200, xend = start, yend = 200), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 300, xend = start, yend = 300), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  
  # Add text showing the value of each 100/75/50/25 lines
  #annotate("text", x = rep(max(data$id),4), y = c(0, 100, 200, 300), label = c("0", "100", "200", "300") , color="grey", size=6 , angle=0, fontface="bold", hjust=1) +
  
  ylim(-160,330) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() +
  
  # Add labels on top of each bar
  geom_text(data=label_data, aes(x=id, y=330, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=6, angle= label_data$angle, inherit.aes = FALSE ) +
  
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -18, label=group), hjust=c(1,1,0), colour = "black", alpha=0.8, size=6, fontface="bold", inherit.aes = FALSE) +
  annotate("text", x = 0, y = -160, label = 'atop(bold("Nsnps_effect"))', parse = TRUE,size=5)+
  annotate("text", x = 0, y = 320, label = 'atop(bold("NParticipants"))', parse = TRUE,size=6)+
  guides(fill=guide_legend(title="Nsnps"))+theme(text = element_text(size=25))
#p1


#################################################################################

# 2 gene, cor=0.4
plot.data<-read.csv("H:/DVMB/Genetic_analysis_SCT_twin/simulations_results2/PT_simulation_results_MC_adj_formatted.csv")

data<-plot.data[plot.data$ngenes==2 &plot.data$eff_size==0.4&plot.data$cluster==0,]

data<-data[,-c(7,8)]

data2<-tidyr::gather(data,"gene_ind","nsig",nsig1:nsig2)

data2<-data.frame(individual.lab=data2$Npats,group=data2$nseff,observation=data2$nsnp, value=data2$nsig,gene_ind=data2$gene_ind)

data2$individual<-paste0(data2$individual.lab,":",data2$group,":",data2$gene_ind,":",data2$observation)
# Transform data in a tidy format (long format)
#data = data %>% gather(key = "observation", value="value", -c(1,2))

data2$group<-as.factor(data2$group)
data2$observation<-as.factor(data2$observation)
data2$individual<-as.factor(data2$individual)
data2$group<-as.factor(data2$group)
data2$gene_ind<-as.factor(data2$gene_ind)



# Set a number of 'empty bar' to add at the end of each group
empty_bar=2
nObsType=nlevels(as.factor(data2$observation))
to_add = data.frame( matrix(NA, empty_bar*nlevels(data2$group)*nObsType, ncol(data2)) )
colnames(to_add) = colnames(data2)
to_add$group=rep(levels(data2$group), each=empty_bar*nObsType )
data2=rbind(data2, to_add)
data2=data2 %>% arrange(group, individual)
data2$id=rep( seq(1, nrow(data2)/nObsType) , each=nObsType)

# Get the name and the y position of each label
label_data= data2 %>% group_by(id, individual) %>% summarize(tot=sum(value))
number_of_bar=nrow(label_data)
angle= 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
#label_data$hjust<-ifelse( angle < -90, 1, 0)
label_data$hjust<-1
label_data$angle<-ifelse(angle < -90, angle+180, angle)

label_data$individual <- substring(label_data$individual,1,3)

is.even <- function(x) x %% 2 == 0

for(k in 1:length(label_data$individual))
{ 
  if(!is.na(label_data$individual[k]))
  {label_data$individual[k]<-ifelse(is.even(label_data$id[k])==FALSE,paste0(label_data$individual[k],"(1)"),paste0(label_data$individual[k],"(2)"))}
  
}

# prepare a data frame for base lines
base_data=data2 %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data = base_data
grid_data$end = grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start = grid_data$start - 1
grid_data=grid_data[-1,]

# Make the plot
p2 = ggplot(data2) +      
  
  # Add the stacked bar
  geom_bar(aes(x=as.factor(id), y=value, fill=observation), stat="identity", alpha=0.5) +
  geom_text(aes(x=as.factor(id), y=value, fill=observation,label=value),size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_viridis(discrete=TRUE) +
  
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 100, xend = start, yend = 100), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 200, xend = start, yend = 200), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 300, xend = start, yend = 300), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  
  # Add text showing the value of each 100/75/50/25 lines
  #annotate("text", x = rep(max(data$id),4), y = c(0, 100, 200, 300), label = c("0", "100", "200", "300") , color="grey", size=6 , angle=0, fontface="bold", hjust=1) +
  
  ylim(-160,330) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() +
  
  # Add labels on top of each bar
  geom_text(data=label_data, aes(x=id, y=330, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=6, angle= label_data$angle, inherit.aes = FALSE ) +
  
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -18, label=group), hjust=c(1,1,0), colour = "black", alpha=0.8, size=6, fontface="bold", inherit.aes = FALSE)+
  annotate("text", x = 0, y = -160, label = 'atop(bold("Nsnps_effect"))', parse = TRUE,size=5)+
  annotate("text", x = 0, y = 320, label = 'atop(bold("NParticipants"))', parse = TRUE,size=6)+
  guides(fill=guide_legend(title="Nsnps"))+theme(text = element_text(size=25))
#p2

#################################################################################

# 4 gene, cor=0.4
plot.data<-read.csv("H:/DVMB/Genetic_analysis_SCT_twin/simulations_results2/PT_simulation_results_MC_adj_formatted.csv")

data<-plot.data[plot.data$ngenes==4 &plot.data$eff_size==0.4&plot.data$cluster==0,]

data2<-tidyr::gather(data,"gene_ind","nsig",nsig1:nsig4)

data2<-data.frame(individual.lab=data2$Npats,group=data2$nseff,observation=data2$nsnp, value=data2$nsig,gene_ind=data2$gene_ind)

data2$individual<-paste0(data2$individual.lab,":",data2$group,":",data2$gene_ind,":",data2$observation)
# Transform data in a tidy format (long format)
#data = data %>% gather(key = "observation", value="value", -c(1,2))

data2$group<-as.factor(data2$group)
data2$observation<-as.factor(data2$observation)
data2$individual<-as.factor(data2$individual)
data2$group<-as.factor(data2$group)
data2$gene_ind<-as.factor(data2$gene_ind)



# Set a number of 'empty bar' to add at the end of each group
empty_bar=2
nObsType=nlevels(as.factor(data2$observation))
to_add = data.frame( matrix(NA, empty_bar*nlevels(data2$group)*nObsType, ncol(data2)) )
colnames(to_add) = colnames(data2)
to_add$group=rep(levels(data2$group), each=empty_bar*nObsType )
data2=rbind(data2, to_add)
data2=data2 %>% arrange(group, individual)
data2$id=rep( seq(1, nrow(data2)/nObsType) , each=nObsType)

# Get the name and the y position of each label
label_data= data2 %>% group_by(id, individual) %>% summarize(tot=sum(value))
number_of_bar=nrow(label_data)
angle= 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
#label_data$hjust<-ifelse( angle < -90, 1, 0)
label_data$hjust<-1
label_data$angle<-ifelse(angle < -90, angle+180, angle)

label_data$individual <- substring(label_data$individual,1,3)

# is.even <- function(x) x %% 2 == 0
# 
# for(k in 1:length(label_data$individual))
# { 
#   if(!is.na(label_data$individual[k]))
#   {label_data$individual[k]<-ifelse(is.even(label_data$id[k])==FALSE,paste0(label_data$individual[k],"(1)"),paste0(label_data$individual[k],"(2)"))}
#   
# }

# prepare a data frame for base lines
base_data=data2 %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data = base_data
grid_data$end = grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start = grid_data$start - 1
grid_data=grid_data[-1,]

# Make the plot
p3 = ggplot(data2) +      
  
  # Add the stacked bar
  geom_bar(aes(x=as.factor(id), y=value, fill=observation), stat="identity", alpha=0.5) +
  geom_text(aes(x=as.factor(id), y=value, fill=observation,label=value),size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_viridis(discrete=TRUE) +
  
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 100, xend = start, yend = 100), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 200, xend = start, yend = 200), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 300, xend = start, yend = 300), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  
  # Add text showing the value of each 100/75/50/25 lines
  #annotate("text", x = rep(max(data$id),4), y = c(0, 100, 200, 300), label = c("0", "100", "200", "300") , color="grey", size=6 , angle=0, fontface="bold", hjust=1) +
  
  ylim(-160,330) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() +
  
  # Add labels on top of each bar
  geom_text(data=label_data, aes(x=id, y=330, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=6, angle= label_data$angle, inherit.aes = FALSE ) +
  
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -18, label=group), hjust=c(1,1,0), colour = "black", alpha=0.8, size=6, fontface="bold", inherit.aes = FALSE)+
  annotate("text", x = 0, y = -160, label = 'atop(bold("Nsnps_effect"))', parse = TRUE,size=5)+
  annotate("text", x = 0, y = 320, label = 'atop(bold("NParticipants"))', parse = TRUE,size=6)+
  guides(fill=guide_legend(title="Nsnps"))+theme(text = element_text(size=25))
#p3

#################################################################################

# 1 gene, cor=0.1
plot.data<-read.csv("H:/DVMB/Genetic_analysis_SCT_twin/simulations_results2/PT_simulation_results_MC_adj_formatted.csv")

data<-plot.data[plot.data$ngenes==1 &plot.data$eff_size==0.1&plot.data$cluster==0,]

data<-data.frame(individual.lab=data$Npats,group=data$nseff,observation=data$nsnp, value=data$nsig1)

data$individual<-paste0(data$individual.lab,":",data$group,":",data$observation)
# Transform data in a tidy format (long format)
#data = data %>% gather(key = "observation", value="value", -c(1,2))

data$group<-as.factor(data$group)
data$observation<-as.factor(data$observation)
data$individual<-as.factor(data$individual)
data$group<-as.factor(data$group)

# Set a number of 'empty bar' to add at the end of each group
empty_bar=2
nObsType=nlevels(as.factor(data$observation))
to_add = data.frame( matrix(NA, empty_bar*nlevels(data$group)*nObsType, ncol(data)) )
colnames(to_add) = colnames(data)
to_add$group=rep(levels(data$group), each=empty_bar*nObsType )
data=rbind(data, to_add)
data=data %>% arrange(group, individual)
data$id=rep( seq(1, nrow(data)/nObsType) , each=nObsType)

# Get the name and the y position of each label
label_data= data %>% group_by(id, individual) %>% summarize(tot=sum(value))
number_of_bar=nrow(label_data)
angle= 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
#label_data$hjust<-ifelse( angle < -90, 1, 0)
label_data$hjust<-1
label_data$angle<-ifelse(angle < -90, angle+180, angle)

label_data$individual <- substring(label_data$individual,1,3)

# prepare a data frame for base lines
base_data=data %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data = base_data
grid_data$end = grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start = grid_data$start - 1
grid_data=grid_data[-1,]

# Make the plot
p4 = ggplot(data) +      
  
  # Add the stacked bar
  geom_bar(aes(x=as.factor(id), y=value, fill=observation), stat="identity", alpha=0.5) +
  geom_text(aes(x=as.factor(id), y=value, fill=observation,label=value),size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_viridis(discrete=TRUE) +
  
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 100, xend = start, yend = 100), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 200, xend = start, yend = 200), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 300, xend = start, yend = 300), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  
  # Add text showing the value of each 100/75/50/25 lines
  #annotate("text", x = rep(max(data$id),4), y = c(0, 100, 200, 300), label = c("0", "100", "200", "300") , color="grey", size=6 , angle=0, fontface="bold", hjust=1) +
  
  ylim(-160,330) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() +
  
  # Add labels on top of each bar
  geom_text(data=label_data, aes(x=id, y=330, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=6, angle= label_data$angle, inherit.aes = FALSE ) +
  
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -18, label=group), hjust=c(1,1,0), colour = "black", alpha=0.8, size=6, fontface="bold", inherit.aes = FALSE)+
  annotate("text", x = 0, y = -160, label = 'atop(bold("Nsnps_effect"))', parse = TRUE,size=5)+
  annotate("text", x = 0, y = 320, label = 'atop(bold("NParticipants"))', parse = TRUE,size=6)+
  guides(fill=guide_legend(title="Nsnps"))+theme(text = element_text(size=25))

#p4


#################################################################################

# 2 gene, cor=0.1
plot.data<-read.csv("H:/DVMB/Genetic_analysis_SCT_twin/simulations_results2/PT_simulation_results_MC_adj_formatted.csv")

data<-plot.data[plot.data$ngenes==2 &plot.data$eff_size==0.1&plot.data$cluster==0,]

data<-data[,-c(7,8)]

data2<-tidyr::gather(data,"gene_ind","nsig",nsig1:nsig2)

data2<-data.frame(individual.lab=data2$Npats,group=data2$nseff,observation=data2$nsnp, value=data2$nsig,gene_ind=data2$gene_ind)

data2$individual<-paste0(data2$individual.lab,":",data2$group,":",data2$gene_ind,":",data2$observation)
# Transform data in a tidy format (long format)
#data = data %>% gather(key = "observation", value="value", -c(1,2))

data2$group<-as.factor(data2$group)
data2$observation<-as.factor(data2$observation)
data2$individual<-as.factor(data2$individual)
data2$group<-as.factor(data2$group)
data2$gene_ind<-as.factor(data2$gene_ind)



# Set a number of 'empty bar' to add at the end of each group
empty_bar=2
nObsType=nlevels(as.factor(data2$observation))
to_add = data.frame( matrix(NA, empty_bar*nlevels(data2$group)*nObsType, ncol(data2)) )
colnames(to_add) = colnames(data2)
to_add$group=rep(levels(data2$group), each=empty_bar*nObsType )
data2=rbind(data2, to_add)
data2=data2 %>% arrange(group, individual)
data2$id=rep( seq(1, nrow(data2)/nObsType) , each=nObsType)

# Get the name and the y position of each label
label_data= data2 %>% group_by(id, individual) %>% summarize(tot=sum(value))
number_of_bar=nrow(label_data)
angle= 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
#label_data$hjust<-ifelse( angle < -90, 1, 0)
label_data$hjust<-1
label_data$angle<-ifelse(angle < -90, angle+180, angle)

label_data$individual <- substring(label_data$individual,1,3)

is.even <- function(x) x %% 2 == 0

for(k in 1:length(label_data$individual))
{ 
  if(!is.na(label_data$individual[k]))
  {label_data$individual[k]<-ifelse(is.even(label_data$id[k])==FALSE,paste0(label_data$individual[k],"(1)"),paste0(label_data$individual[k],"(2)"))}
  
}

# prepare a data frame for base lines
base_data=data2 %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data = base_data
grid_data$end = grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start = grid_data$start - 1
grid_data=grid_data[-1,]

# Make the plot
p5 = ggplot(data2) +      
  
  # Add the stacked bar
  geom_bar(aes(x=as.factor(id), y=value, fill=observation), stat="identity", alpha=0.5) +
  geom_text(aes(x=as.factor(id), y=value, fill=observation,label=value),size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_viridis(discrete=TRUE) +
  
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 100, xend = start, yend = 100), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 200, xend = start, yend = 200), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 300, xend = start, yend = 300), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  
  # Add text showing the value of each 100/75/50/25 lines
  #annotate("text", x = rep(max(data$id),4), y = c(0, 100, 200, 300), label = c("0", "100", "200", "300") , color="grey", size=6 , angle=0, fontface="bold", hjust=1) +
  
  ylim(-160,330) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() +
  
  # Add labels on top of each bar
  geom_text(data=label_data, aes(x=id, y=330, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=6, angle= label_data$angle, inherit.aes = FALSE ) +
  
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -18, label=group), hjust=c(1,1,0), colour = "black", alpha=0.8, size=6, fontface="bold", inherit.aes = FALSE)+
  annotate("text", x = 0, y = -160, label = 'atop(bold("Nsnps_effect"))', parse = TRUE,size=5)+
  annotate("text", x = 0, y = 320, label = 'atop(bold("NParticipants"))', parse = TRUE,size=6)+
  guides(fill=guide_legend(title="Nsnps"))+theme(text = element_text(size=25))
#p5
#################################################################################

# 4 gene, cor=0.1
plot.data<-read.csv("H:/DVMB/Genetic_analysis_SCT_twin/simulations_results2/PT_simulation_results_MC_adj_formatted.csv")

data<-plot.data[plot.data$ngenes==4 &plot.data$eff_size==0.1&plot.data$cluster==0,]

data2<-tidyr::gather(data,"gene_ind","nsig",nsig1:nsig4)

data2<-data.frame(individual.lab=data2$Npats,group=data2$nseff,observation=data2$nsnp, value=data2$nsig,gene_ind=data2$gene_ind)

data2$individual<-paste0(data2$individual.lab,":",data2$group,":",data2$gene_ind,":",data2$observation)
# Transform data in a tidy format (long format)
#data = data %>% gather(key = "observation", value="value", -c(1,2))

data2$group<-as.factor(data2$group)
data2$observation<-as.factor(data2$observation)
data2$individual<-as.factor(data2$individual)
data2$group<-as.factor(data2$group)
data2$gene_ind<-as.factor(data2$gene_ind)



# Set a number of 'empty bar' to add at the end of each group
empty_bar=2
nObsType=nlevels(as.factor(data2$observation))
to_add = data.frame( matrix(NA, empty_bar*nlevels(data2$group)*nObsType, ncol(data2)) )
colnames(to_add) = colnames(data2)
to_add$group=rep(levels(data2$group), each=empty_bar*nObsType )
data2=rbind(data2, to_add)
data2=data2 %>% arrange(group, individual)
data2$id=rep( seq(1, nrow(data2)/nObsType) , each=nObsType)

# Get the name and the y position of each label
label_data= data2 %>% group_by(id, individual) %>% summarize(tot=sum(value))
number_of_bar=nrow(label_data)
angle= 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
#label_data$hjust<-ifelse( angle < -90, 1, 0)
label_data$hjust<-1
label_data$angle<-ifelse(angle < -90, angle+180, angle)

label_data$individual <- substring(label_data$individual,1,3)

# is.even <- function(x) x %% 2 == 0
# 
# for(k in 1:length(label_data$individual))
# { 
#   if(!is.na(label_data$individual[k]))
#   {label_data$individual[k]<-ifelse(is.even(label_data$id[k])==FALSE,paste0(label_data$individual[k],"(1)"),paste0(label_data$individual[k],"(2)"))}
#   
# }

# prepare a data frame for base lines
base_data=data2 %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data = base_data
grid_data$end = grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start = grid_data$start - 1
grid_data=grid_data[-1,]

# Make the plot
p6 = ggplot(data2) +      
  
  # Add the stacked bar
  geom_bar(aes(x=as.factor(id), y=value, fill=observation), stat="identity", alpha=0.5) +
  geom_text(aes(x=as.factor(id), y=value, fill=observation,label=value),size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_viridis(discrete=TRUE) +
  
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 100, xend = start, yend = 100), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 200, xend = start, yend = 200), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 300, xend = start, yend = 300), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  
  # Add text showing the value of each 100/75/50/25 lines
  #annotate("text", x = rep(max(data$id),4), y = c(0, 100, 200, 300), label = c("0", "100", "200", "300") , color="grey", size=6 , angle=0, fontface="bold", hjust=1) +
  
  ylim(-160,330) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() +
  
  # Add labels on top of each bar
  geom_text(data=label_data, aes(x=id, y=330, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=6, angle= label_data$angle, inherit.aes = FALSE ) +
  
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -18, label=group), hjust=c(1,1,0), colour = "black", alpha=0.8, size=6, fontface="bold", inherit.aes = FALSE)+
  annotate("text", x = 0, y = -160, label = 'atop(bold("Nsnps_effect"))', parse = TRUE,size=5)+
  annotate("text", x = 0, y = 320, label = 'atop(bold("NParticipants"))', parse = TRUE,size=6)+
  guides(fill=guide_legend(title="Nsnps"))+theme(text = element_text(size=25))
#p6

#################################################################################################
#################################################################################################
#################################################################################################


png(filename="c:/Users/pthompson/Dropbox/project SCT analysis/GSCA validation/draft article//dials_nocluster_MC.png",height=20,width=20,units="in",res=100)
#windows()
big_p<-ggpubr::ggarrange(p1,p2,p3,p4,p5,p6,ncol=3,nrow=2, labels = c("A (Ngene=1; r=0.4)", "B (Ngene=2; r=0.4)","C (Ngene=4; r=0.4)","D (Ngene=1; r=0.1)", "E (Ngene=2; r=0.1)","F (Ngene=4; r=0.1)"),
                         common.legend = TRUE, legend = "bottom",hjust=-0.6,font.label = list(size = 25, face = "bold"))


library(grid)
print(big_p, vp=viewport(angle=0))
dev.off()

file.show("c:/Users/pthompson/Dropbox/project SCT analysis/GSCA validation/draft article//dials_nocluster_MC.png")

