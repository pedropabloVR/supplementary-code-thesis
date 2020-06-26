# Script for plotting colocalisation and cluster analaysis results from the compare_conditions.m script from SynaptoAnalysis. 

# This script can be used to replicate the figures in Part 1 Chapter 2 of my thesis 

# The input is a csv file with the results from the overlap calculations from synapto-analysis.  

# Pedro Vallejo Ramirez

# Created: 10/03/2019
# Updated: 28/05/2020


## This script uses a standard theme for the plots, with a grey-scale background, 10 pt font sizes, and margins for all text

## Formatting is standard for biology-style plots, with grayscale scatter plots and 
## a minimalist axis.

library(tidyverse)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ggbeeswarm)
library(dplyr)
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
  tick_size <- 8
  
  # Begin construction of chart
  theme_bw(base_size=11) +
    
    # Format the legend, uncomment position if need be
    theme(legend.text = element_text(size=font_size, color=color.axis.title)) +
    theme(legend.title =element_text(size=font_size, face = "bold", color=color.axis.title)) +
    theme(legend.position="bottom")+
    
    # Set title and axis labels, and format these and tick marks
    theme(plot.title=element_text(color=color.title, size=font_size, vjust=1.25, hjust = 0.4, face = "bold")) +
    theme(axis.text.x=element_text(size=tick_size,color=color.axis.text, margin = unit(c(0.1, 0.5, 0.3, 0.5), "cm"))) +
    theme(axis.text.y=element_text(size=tick_size,color=color.axis.text, margin = unit(c(0.5, 0.1, 0.5, 0.5), "cm"))) +
    theme(axis.title.x=element_text(size=tick_size,color=color.axis.title, 
                                    vjust=1.25, margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))) +
    theme(axis.title.y=element_text(size=tick_size,color=color.axis.title, 
                                    vjust=-4, hjust = 0.5, margin = unit(c(0.2, 0.2, 0.2, 
                                                                                          0.2), "cm"))) 
  #theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 
  #0.25)) 
} # plotting function


# Import data tables post-overlap threshold
data_table_37C_old <- read.csv("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/results datasets for synaptosomes for colocalisation and ripley analysis/results_combined_after_overlap_threshold_37C_20190917_trimmed.csv")
data_table_4C_old <- read.csv("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/results datasets for synaptosomes for colocalisation and ripley analysis/results_combined_after_overlap_threshold_4C_20190915_trimmed.csv")

data_table_37C <- read.csv("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/results datasets for synaptosomes for colocalisation and ripley analysis/results_combined_after_overlap_threshold_37C_20200526_trimmed.csv")
data_table_4C <- read.csv("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/results datasets for synaptosomes for colocalisation and ripley analysis/results_combined_after_overlap_threshold_4C_20200526_trimmed.csv")

# Plotting overall overlap coefficient (WOC) for a-syn and VAMP2 at different temperatures
{
## To plot weighted averages at 4C ####

# select data columns
df <- select(data_table_4C,Repeat,condition,WeightedOverlapWithGreen,WeightedOverlapWithBlue)
# Organize into tidy form
df_melt = melt(df)
# Change labels for plotting 
df2_4C <- data.frame(Repeat = df_melt[,1],condition = df_melt[,2],variable = c(rep("mCLING on asyn",nrow(data_table_4C)),rep("mCLING on VAMP2",nrow(data_table_4C))),value = df_melt[,4])

# Using ggplot2, generate dot plot
ggplot(df2_4C,aes(x = condition, y = value, color = variable,fill = variable)) + 
  geom_jitter( alpha = 0.7)+
  #geom_violin(alpha = 0.3)+
  geom_boxplot(alpha = 0.2)+
  #scale_color_manual(values = c("blue","green4"))+
  xlab("Compared data sets")+
  ylab("Spatial co-ocurrence")+
  ylim(0,1) + 
  ggtitle("Incubation at 4C")+
  prism_multi()

## To plot weighted averages at 37C ####

# select data columns
df <- select(data_table_37C,Repeat,condition,WeightedOverlapWithGreen,WeightedOverlapWithBlue)
# Organize into tidy form
df_melt = melt(df)
# Change labels for plotting 
df2_37C <- data.frame(Repeat = df_melt[,1],condition = df_melt[,2],variable = c(rep("mCLING on asyn",nrow(data_table_37C)),rep("mCLING on VAMP2",nrow(data_table_37C))),value = df_melt[,4])

# Using ggplot2, generate dot plot
ggplot(df2_37C,aes(x = condition, y = value, color = variable,fill = variable)) + 
  geom_jitter( alpha = 0.7)+
  #geom_violin(alpha = 0.3)+
  geom_boxplot(alpha = 0.2)+
  #scale_color_manual(values = c("blue","green4"))+
  xlab("Compared data sets")+
  ylab("Spatial co-ocurrence")+
  ylim(0,1) + 
  ggtitle("Incubation at 37C")+
  prism_multi()

## To pool results at 4C and 37C and plot together ####

# Select data values from 4C and 37C data

## Can't find a nice way to copy values from the variable columns of df2 4c and df2 37c
#df3 <- data.frame(Temperature = c(rep("4C",nrow(data_table_4C)*2),rep("37C",nrow(data_table_37C)*2)),variable = c(df2_4C[,3],df2_37C[,3]),value =c(df2_4C[,4],df2_37C[,4]) )
df3 <- data.frame(Temperature = c(rep("4C",nrow(data_table_4C)*2),rep("37C",nrow(data_table_37C)*2)),
                  variable = c(rep("mCLING on asyn",nrow(data_table_4C)),rep("mCLING on VAMP2",nrow(data_table_4C)),
                               rep("mCLING on asyn",nrow(data_table_37C)),rep("mCLING on VAMP2",nrow(data_table_37C))),
                  value =c(df2_4C[,4],df2_37C[,4]))


# Using ggplot2, generate dot plot (plotting parameters for paper)
ggplot(df3,aes(x = Temperature, y = value, color = variable,fill = variable,shape = variable)) + 
  geom_quasirandom(alpha = 1.0, size = 1,dodge.width = 0.8)+
  geom_boxplot(alpha = 0.6,outlier.shape = NA)+
  xlab("Temperature conditions")+
  ylab("Overall co-occurrence coefficient")+
  ylim(0,1) + 
  prism_multi()+
  scale_color_manual(values = c("#619CFF","#F8766D"))+
  scale_fill_manual(values = c("#619CFF","#F8766D"))+
  theme(legend.position = "none")
  #stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
  #             width = 0.75, size = 1,linetype = "solid", position = position_dodge2(width = 0.8))+
  #scale_fill_brewer(palette = "Set2")+
  #scale_color_brewer(palette = "Set2")+

  
ggsave("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/SUPERPLOTS for colocalisation and cluster size/WOC_diffColor.pdf",width = 75, height = 63, units = "mm")


## T-test to test the probability of the data with the null hypothesis that the difference between means is zero. 

# compare asyn at 4C to asyn at 37c
asyn   = subset(df3,variable=="mCLING on asyn")

x_asyn = subset(asyn, Temperature == "4C")
y_asyn = subset(asyn, Temperature == "37C")

x = x_asyn %>% select(value) # all values at 4C for mcling on asyn
y = y_asyn %>% select(value) # all values at 37C for mcling on asyn

t.test(x,y)

# compare vamp2 at 4C to vamp2 at 37c
vamp2  = subset(df3,variable=="mCLING on VAMP2")
x_vamp2 = subset(vamp2, Temperature == "4C")
x = vamp2 %>% select(value) # all values at 4C for mcling on asyn

y_vamp2 = subset(vamp2, Temperature == "37C")
y = y_vamp2 %>% select(value) # all values at 37C for mcling on asyn

t.test(x,y)
}

# Plot the individual mander's coefficients for mCLING-asyn, asyn-mCLING, mCLING-VAMP2, VAMP2-mCLING

{
# set up data frames
df_manders_37C <- select(data_table_37C,Repeat,condition,MandersRG,MandersGR,MandersRB,MandersBR,pairwiseRG,pairwiseRB,pairwiseGB)
df_manders_4C  <- select( data_table_4C,Repeat,condition,MandersRG,MandersGR,MandersRB,MandersBR,pairwiseRG,pairwiseRB,pairwiseGB)
# add column for temperature
df_manders_37C$Temperature <- c(rep("37C",nrow(df_manders_37C)))
df_manders_4C$Temperature  <- c(rep("4C", nrow(df_manders_4C)))
# create master data frame
df_manders_master <- full_join(df_manders_37C,df_manders_4C)

# get individual colours from the Rcolorbrewer Set2 palette
my_palette = c(brewer.pal(5, "Set2")[c(1,2,3,4,5)])
#grid::grid.raster(my_palette, int=F)

scale_colour_discrete = function(...) scale_colour_manual(..., values = palette())

palette(my_palette)
p # custom colors


# mCLING-asyn as a function of temperature, single colour 
  ggplot(df_manders_master,aes(x = Temperature, y = MandersGR, colour = Temperature, shape = Temperature)) +
    geom_boxplot(alpha = 0.6,outlier.shape = NA,fatten=NULL)+
    geom_quasirandom(alpha = 1, size = 2.0,dodge.width = 0.8)+
    scale_fill_manual(values = c(my_palette[1],my_palette[1]))+
    scale_color_manual(values = c(my_palette[1],my_palette[1]))+
    xlab("Temperature")+
    ylab("Fractional overlap asyn-mCLING")+
    prism_multi()+
    stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.75, size = 1,linetype = "solid", position = position_dodge2(width = 0.8))+
    theme(legend.position = "none")
  
  # as an estimation graphic
  unpaired_mean_diff_M1_asyn <- dabest(df_manders_master,Temperature,MandersRG,
                                        idx = c("37C","4C"), # add the egta K for a multi group plot
                                        paired = FALSE)
  
  mk_size = 1;
  plot(unpaired_mean_diff_M1_asyn,
       rawplot.ylabel = "Fractional overlap mCLING-asyn",
       rawplot.markersize = mk_size,
       palette = "Set1",
       tick.fontsize  = 8, axes.title.fontsize = 8, effsize.markersize = mk_size)
       #swarmplot.params = list(shape = c(1) ) )
  
  ggsave("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/SUPERPLOTS for colocalisation and cluster size/M1_mCLING_asyn_ggsave.pdf",width = 75, height = 60, units = "mm")
  

# mCLING-VAMP2 as a function of temperature, colored by calcium
  
  ggplot(df_manders_master,aes(x = Temperature, y = MandersBR,color = Temperature,shape = Temperature)) +
    geom_quasirandom(alpha = 1, size = 2.0,dodge.width = 0.8)+
    geom_boxplot(alpha = 0.6,outlier.shape = NA,fatten=NULL)+
    scale_fill_manual(values = c(my_palette[2],my_palette[2]))+
    scale_color_manual(values = c(my_palette[2],my_palette[2]))+
    xlab("Temperature")+
    ylab("Fractional overlap VAMP2-mCLING")+
    prism_multi()+
    stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.75, size = 1,linetype = "solid", position = position_dodge2(width = 0.8))+
    theme(legend.position = "none")

  unpaired_mean_diff_M1_VAMP2 <- dabest(df_manders_master,Temperature,MandersRB,
                                       idx = c("37C","4C"), # add the egta K for a multi group plot
                                       paired = FALSE)
  

  plot(unpaired_mean_diff_M1_VAMP2, 
       rawplot.ylabel = "Fractional overlap mCLING-VAMP2",
       rawplot.markersize = mk_size,
       palette = "Accent",
       tick.fontsize  = 8, axes.title.fontsize = 8, effsize.markersize = mk_size)
  #swarmplot.params = list(shape = c(1) ) ) # trying to plot these with different shapes but it doesn't work!!
  
  ggsave("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/SUPERPLOTS for colocalisation and cluster size/M1_mCLING_VAMP2_ggsave.pdf",width = 75, height = 60, units = "mm")
} 
  
# Plot the intercluster distances between mCLING-asyn, mCLING-VAMP2 (measured as the pairwise distance between cluster centroids)

{
# Pairwise distances for mCLING-asyn
  ggplot(df_manders_master,aes(x = Temperature, y = pairwiseRG,color = Repeat)) +
    geom_quasirandom(alpha = 1, size = 2.0,dodge.width = 0.8)+
    geom_boxplot(alpha = 0.6,outlier.shape = NA,fatten=NULL)+
    xlab("Temperature")+
    ylab("inter-cluster distance mCLING-asyn")+
    prism_multi()+
    stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.75, size = 1,linetype = "solid", position = position_dodge2(width = 0.8))+
    scale_fill_brewer(palette = "Set2")+
    scale_color_brewer(palette = "Set2")+
    ylim(0,400)
  
  unpaired_mean_diff_intercluster_asyn <- dabest(df_manders_master,Temperature,pairwiseRG,
                                       idx = c("37C","4C"), # add the egta K for a multi group plot
                                       paired = FALSE)
  
  mk_size = 1;
  plot(unpaired_mean_diff_intercluster_asyn,
       rawplot.ylabel = "Intercluster distance mCLING-asyn",
       rawplot.markersize = mk_size,
       palette = "Accent",
       tick.fontsize  = 8, axes.title.fontsize = 8, effsize.markersize = mk_size)

  
  # Pairwise distances for mCLING-VAMP2
  ggplot(df_manders_master,aes(x = Temperature, y = pairwiseGB,color = Repeat)) +
    geom_quasirandom(alpha = 1, size = 2.0,dodge.width = 0.8)+
    geom_boxplot(alpha = 0.6,outlier.shape = NA,fatten=NULL)+
    xlab("Temperature")+
    ylab("inter-cluster distance mCLING-VAMP2")+
    prism_multi()+
    stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.75, size = 1,linetype = "solid", position = position_dodge2(width = 0.8))+
    scale_fill_brewer(palette = "Set2")+
    scale_color_brewer(palette = "Set2")+
    ylim(0,400)
  
  unpaired_mean_diff_intercluster_VAMP2 <- dabest(df_manders_master,Temperature,pairwiseRB,
                                                 idx = c("37C","4C"), # add the egta K for a multi group plot
                                                 paired = FALSE)
  
  mk_size = 1;
  plot(unpaired_mean_diff_intercluster_VAMP2,
       rawplot.ylabel = "Intercluster distance mCLING-VAMP2",
       rawplot.markersize = mk_size,
       palette = "Set1",
       tick.fontsize  = 8, axes.title.fontsize = 8, effsize.markersize = mk_size)
  #swarmplot.params = list(shape = c(1) ) )
  
  
  # the measurement of intercluster distance using the pairwise distances from centroid are not informative 
  }
  
  
#### Plotting the ripley's K cluster size results from the master data set 
{
# Make master data frame with info from ripley's K function
df_ripley_37C <- select(data_table_37C,Repeat,condition,clustersizeGC,clustersizeBC,interclusterdistRG,interclusterdistRB,interclusterdistGB)
df_ripley_4C  <- select( data_table_4C,Repeat,condition,clustersizeGC,clustersizeBC,interclusterdistRG,interclusterdistRB,interclusterdistGB)

df_ripley_37C$Temperature <- c(rep("37C",nrow(df_ripley_37C)))
df_ripley_4C$Temperature  <- c(rep("4C", nrow(df_ripley_4C)))

df_ripley_master <- full_join(df_ripley_37C,df_ripley_4C)

# Order the different conditions to make it easier to display control vs treatments
df_ripley_master$condition <- factor(df_ripley_master$condition, levels = c("phys","egta","egtak"))


# A-syn cluster size as a function of calcium, colored by temperature 
{

  # Superplot 1: Group by condition, dodge and colour by temperature
  ggplot(df_ripley_master,aes(x = condition, y = clustersizeGC,color = factor(Temperature))) +
    geom_quasirandom(alpha = 0.4, size = 2.0,dodge.width = 0.8)+
    #geom_boxplot(alpha = 0.6,outlier.shape = NA,fatten=NULL)+
    xlab("Conditions")+
    ylab("Cluster size")+
    prism_multi()+
    stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = 0.75, size = 1,linetype = "solid", position = position_dodge2(width = 0.8))

  # Superplot 2: Group by condition, dodge scatterplots by repeat, and colour them by temperature
  ggplot(df_ripley_master,aes(x = condition, y = clustersizeGC, group=factor(Repeat),color = Temperature ))  + 
    geom_quasirandom(alpha = 0.4, size = 2.0,dodge.width=0.8)+
    #scale_shape_manual(values = c(2,0,21))+
    #geom_boxplot(alpha = 0.2,outlier.shape = NA,fatten = NULL)+
    xlab("Calcium treatment")+
    ylab("RMS distance (nm) from centroid") + 
    prism_multi()+
    stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.75, size = 1,linetype = "solid", position = position_dodge2(width = 0.8))
  

# Estimation graphics to visualize effect size 
library(dabestr)      # for plotting the differences between distributions as an estimation graphic


unpaired_mean_diff_asynRipley <- dabest(df_ripley_master,condition,clustersizeGC,
                                  idx = c("phys","egta","egtak"), # add the egta K for a multi group plot
                                  paired = FALSE)

plot(unpaired_mean_diff_asynRipley, rawplot.ylabel = "cluster size value",effsize.ylim = c(-100, 100), palette = "Set2")

p1_ripley <- plot(unpaired_mean_diff_asynRipley,
     rawplot.ylabel = "Cluster radius (nm)",
     rawplot.markersize = mk_size,
     palette = "Set2",
     tick.fontsize  = 8, axes.title.fontsize = 8,effsize.markersize = 2)
p1_ripley
#swarmplot.params = list(shape = c(1) ) )

ggsave("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/results cluster size/clusterSize_asyn_RipleyAll.pdf",width = 80, height = 80, units = "mm")



# t test
x_cluster <- subset(df_ripley_master,condition == 'phys')
y_cluster <- subset(df_ripley_master,condition == 'egta')
x = x_cluster %>% select(clustersizeGC)
y = y_cluster %>% select(clustersizeGC)
t.test(x,y) # p-value = 0.0331

#  mean of x mean of y 
#  451.0448  512.9412

}


# VAMP2 cluster size as a function of calcium, colored by temperature 
{
  # VAMP2
  # Superplot 1: Group by condition, dodge and colour by temperature
  ggplot(df_ripley_master,aes(x = condition, y = clustersizeBC,color = factor(Temperature))) +
    geom_quasirandom(alpha = 0.4, size = 2.0,dodge.width = 0.8)+
    #geom_boxplot(alpha = 0.6,outlier.shape = NA,fatten=NULL)+
    xlab("Conditions")+
    ylab("Cluster size")+
    prism_multi()+
    stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.75, size = 1,linetype = "solid", position = position_dodge2(width = 0.8))
  
  # Superplot 2: Group by condition, dodge scatterplots by repeat, and colour them by temperature
  ggplot(df_ripley_master,aes(x = condition, y = clustersizeBC, group=factor(Repeat),color = Temperature ))  + 
    geom_quasirandom(alpha = 0.4, size = 2.0,dodge.width=0.8)+
    #scale_shape_manual(values = c(2,0,21))+
    #geom_boxplot(alpha = 0.2,outlier.shape = NA,fatten = NULL)+
    xlab("Calcium treatment")+
    ylab("RMS distance (nm) from centroid") + 
    prism_multi()+
    stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.75, size = 1,linetype = "solid", position = position_dodge2(width = 0.8))
  
  # Estimation graphics to visualize effect size 
  unpaired_mean_diff_VAMP2ripley <- dabest(df_ripley_master,condition,clustersizeBC,
                                    idx = c("phys","egta","egtak"), # add the egta K for a multi group plot
                                    paired = FALSE)
  
  plot(unpaired_mean_diff_VAMP2ripley, rawplot.ylabel = "cluster size value",palette= "Set2")
  
  p2_ripley <- plot(unpaired_mean_diff_VAMP2ripley,
       rawplot.ylabel = "Cluster radius (nm)",
       rawplot.markersize = mk_size,
       palette = "Accent",
       tick.fontsize  = 8, axes.title.fontsize = 8,
       effsize.markersize = 2)
  #swarmplot.params = list(shape = c(1) ) )
  p2_ripley
  
  ggsave("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/results cluster size/clusterSize_VAMP2_RipleyAll.pdf",width = 80, height = 80, units = "mm")
  
  
  # t test
  x_cluster <- subset(df_ripley_master,condition == 'phys')
  y_cluster <- subset(df_ripley_master,condition == 'egta')
  x = x_cluster %>% select(clustersizeBC)
  y = y_cluster %>% select(clustersizeBC)
  t.test(x,y) # p-value = 0.01048
  
  # mean of x mean of y 
  # 495.6716  568.2353
  
}

# using cow plot to arrange Gardner-Altman plots in a grid

plot_grid(
  p1_ripley, p2_ripley,
  labels = "AUTO", ncol = 1,vjust = 0
)
ggsave("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/results cluster size/clusterSize_Ripley_Grid.pdf",width = 150, height = 150, units = "mm")



}

# Peak from the bivariate Ripley's K analysis (not used currently)
{
ggplot(df_ripley_master,aes(x = Temperature, y = interclusterdistRB,color = factor(Temperature))) +
  geom_quasirandom(alpha = 0.4, size = 2.0,dodge.width = 0.8)+
  #geom_boxplot(alpha = 0.6,outlier.shape = NA,fatten=NULL)+
  xlab("Conditions")+
  ylab("Interaction length mCLING-asyn")+
  prism_multi()+
  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = 0.75, size = 1,linetype = "solid", position = position_dodge2(width = 0.8))


df_37C_interclusterdistRG = subset(df_ripley_master, Temperature == '37C')

mean( subset(df_ripley_master, Temperature == '37C')[,'interclusterdistRG'] )
mean( subset(df_ripley_master, Temperature == '4C')[,'interclusterdistRG'] )

mean( subset(df_ripley_master, Temperature == '37C')[,'interclusterdistRB'] )
mean( subset(df_ripley_master, Temperature == '4C')[,'interclusterdistRB'] )

}

# Deprecated code
{
#data_table_4C <- read.csv("E:\\Experiments\\synaptosomes\\analysis_20190305\\Results_combined_4C\\Repeats_combined\\results_combined_after_overlap_threshold_trimmed4C.csv")
#data_table_37C <- read.csv("E:\\Experiments\\synaptosomes\\analysis_20190305\\Results_combined_37C\\Repeats_combined\\results_combined_after_overlap_threshold_trimmed37C.csv")


# Using ggplot2, generate dot plot (plotting parameters for poster)
ggplot(df3,aes(x = Temperature, y = value, color = variable,fill = variable,shape = variable)) + 
  geom_quasirandom(alpha = 0.4, size = 1.5,dodge.width = 0.8)+
  geom_boxplot(alpha = 0.8,lwd = 2)+
  xlab("Temperature conditions")+
  ylab("Spatial co-ocurrence")+
  ylim(0,1) + 
  prism_multi()
}

