#----Supplementary code accompanying Part1 Chapter2 of PVR thesis----#

# Script for plotting colocalisation and cluster size of model 1 and model 2 from the synaptosome model generated in testSTORM. 

# The input is a csv file with the results from SynaptoAnalysis, containing the WOC, Manders M1 and M2 coefficients, and the Ripley's K cluster size results. 
# A .csv file with the results from RMSD analysis is also required to plot the RMSD values

# Pedro Vallejo Ramirez
# Laser Analytics Group, University of Cambridge
# Created: 07/03/2019


## This script uses a standard theme for the plots, with a grey-scale background, 10 pt font sizes, and margins for all text

## Formatting is standard for biology-style plots, with grayscale scatter plots and 
## a minimalist axis.

library(tidyverse)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
set.seed(2019) # for reproducible plots
prism_multi <- function() {
  
  # Generate the colors for the chart procedurally with RColorBrewer
  palette <- brewer.pal("Greys", n=9)
  color.background = "white"
  color.grid.major = "white"
  color.axis.text = "black"
  color.axis.title = "black"
  color.title = palette[9]
  font_size <- 10
  title_size <- 12
  
  # Begin construction of chart
  theme_bw(base_size=11) +
    
    # Format the legend, uncomment position if need be
    theme(legend.text = element_text(size=font_size, color=color.axis.title)) +
    theme(legend.title =element_text(size=font_size, face = "bold", color=color.axis.title)) +
    theme(legend.position="bottom")+
    
    # Set title and axis labels, and format these and tick marks
    theme(plot.title=element_text(color=color.title, size=font_size, vjust=1.25, hjust = 0.4, face = "bold")) +
    theme(axis.text.x=element_text(size=font_size,color=color.axis.text, margin = unit(c(0.1, 0.5, 0.3, 0.5), "cm"))) +
    theme(axis.text.y=element_text(size=font_size,color=color.axis.text, margin = unit(c(0.5, 0.1, 0.5, 0.5), "cm"))) +
    theme(axis.title.x=element_text(size=font_size,color=color.axis.title, 
                                    vjust=1.25, face = "bold", margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))) +
    theme(axis.title.y=element_text(size=font_size,color=color.axis.title, 
                                    vjust=-4, hjust = 0.5, face = "bold", margin = unit(c(0.2, 0.2, 0.2, 
                                                                                          0.2), "cm"))) 
  #theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 
  #0.25)) 
}

# Import data with colocalisation and Ripley's K results
# df <- read.csv("D:\\Experiments\\synaptosomes\\testSTORM20190218\\summary_results20190307.csv")
data_model1 <- read.csv("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/figure synapto simulation colocalisation/summary_SynaptoAnalysis_results_model1.csv")
data_model3 <- read.csv("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/figure synapto simulation colocalisation/summary_SynaptoAnalysis_results_model3.csv")

## Colocalisation analysis 
{
# organize colocalisation results in long form (tidy format)
data_model1 <- subset(data_model1,condition == "phys")

df_model1 <- select(data_model1,WeightedOverlapWithGreen, WeightedOverlapWithBlue, MandersRG, MandersGR,MandersRB,MandersBR)
df_model3 <- select(data_model3,WeightedOverlapWithGreen, WeightedOverlapWithBlue, MandersRG, MandersGR,MandersRB,MandersBR)

df_model1$model <- c(rep("model1",nrow(df_model1)))
df_model3$model <- c(rep("model3",nrow(df_model3)))

df_model_master <- full_join(df_model1,df_model3)
# rename variables for plotting
df_model_master2 <- rename(df_model_master,"WOC mCLING-asyn" = WeightedOverlapWithGreen,"WOC mCLING-VAMP2" = WeightedOverlapWithBlue,"M1 mCLING-asyn" = MandersRG,"M2 asyn-mCLING" = MandersGR,"M1 mCLING-VAMP2" = MandersRB,"M2 VAMP2-mCLING" = MandersBR)
df_model_melt <- melt(df_model_master2) # works so well! I love this function


# Using ggplot2, generate dot plot
ggplot(df_model_melt,aes(x = model, y = value, color = variable,shape = variable)) + 
  geom_quasirandom( size=1.0, alpha = 0.8, dodge.width = 0.8)+
  #scale_color_manual(values = c("blue","green4"))+
  xlab("Model")+
  ylab("Spatial co-occurrence")+
 # ylim(0,1) + 
  prism_multi()+
  scale_fill_brewer(palette = "Set1")+
  scale_color_brewer(palette = "Set1")
  #theme(legend.position = "None")

# save to drive
ggsave("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/figure synapto simulation colocalisation/colocalisationSimulation_2Models.pdf", width = 140, height = 80, units = "mm")
}

## Cluster size analysis
{
# import data from RMSD analysis

# office desktop location : data_table <- read.csv("D:\\Experiments\\synaptosomes\\testSTORM20190218\\model 1\\rmsd\\Repeat_A\\results_rmsd.csv")
data_rmsd <- read.csv("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/figure synapto simulation cluster size/results_rmsd_model1.csv")
data_rmsd <- subset(data_rmsd,condition=="phys") # calcium conditions here are all the same, select first

# tidy up and combine into master data frame
df_rmsd <- select(data_rmsd,rmsd_RC,rmsd_GC, rmsd_BC)
df_rmsd$metric = c(rep("RMSD",nrow(df_rmsd)))
df_rmsd <- rename(df_rmsd, "mCLING"=rmsd_RC,"a-syn" =rmsd_GC, "VAMP-2"=rmsd_BC)

df_ripley <-  select(data_model1,clustersizeRC, clustersizeGC,clustersizeBC)
df_ripley$metric =  c(rep("Ripley's L(r)-r",nrow(df_model3))) # add column with label for cluster size metric
df_ripley <- rename(df_ripley, "mCLING"=clustersizeRC,"a-syn" =clustersizeGC, "VAMP-2"=clustersizeBC)

df_clustersize_master = full_join(df_rmsd,df_ripley)
df_clustersize_melt = melt(df_clustersize_master)

df$channel <- factor(df$channel, levels = c("mCLING","asyn","VAMP2"))


ggplot(df_clustersize_melt,aes(x = variable, y = value, color = metric,shape = metric)) + 
  geom_quasirandom(alpha = 0.8, size = 2.0,dodge.width=0.8)+
  geom_boxplot(alpha = 0.2,outlier.shape = NA)+
  ylab("Cluster size (nm)")+
  prism_multi()+
  scale_fill_brewer(palette = "Set1")+
  scale_color_brewer(palette = "Set1")+
  theme(legend.position = "None")
  
  
ggsave("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/figure synapto simulation cluster size/clustersize_Ripley_RMSD.pdf",width = 140, height = 80, units = "mm")

# mean ripley values for mcling, asyn, and VAMP-2

mean_mCLING <- mean(df_ripley[,1])
mean_asyn <- mean(df_ripley[,2])
mean_VAMP2 <- mean(df_ripley[,3])
}

