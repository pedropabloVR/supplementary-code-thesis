# Script for plotting results from RMSD analysis pooled from all synaptosome data colleted from analysis on 2019_03_05
# The data set is the same as that used to test for co-occurrence of a-syn and VAMP2
# Pedro Vallejo Ramirez
# Created: 07/08/2019
# Updated: 25/05/2020


## This script uses a standard theme for the plots, with a grey-scale background, 10 pt font sizes, and margins for all text

## Formatting is standard for biology-style plots, with grayscale scatter plots and 
## a minimalist axis.

library(tidyverse)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ggbeeswarm)
library(dplyr)
library(car) # this one contains the Levene test to verify the assumptions in the ANOVA test. 

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
                                    vjust=1.25, margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))) +
    theme(axis.title.y=element_text(size=font_size,color=color.axis.title, 
                                    vjust=-4, hjust = 0.5, margin = unit(c(0.2, 0.2, 0.2, 
                                                                           0.2), "cm"))) 
  #theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 
  #0.25)) 
} # plotting function


## Import data tables

# 2019 03 data set
#data_table_rmsd <- read.csv("E:\\OneDrive\\OneDrive - University Of Cambridge\\lag\\microscopy work\\synaptosomes\\rmsd analysis 20190807 summary\\37C_4C_rmsd_combined.csv")
#data_table_rmsd <- read.csv("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/rmsd analysis 20190807 summary/37C_4C_rmsd_combined.csv")

# 2019 09 15 data set
#data_table_rmsd <- read.csv("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/rmsd analysis 20190915 data/pooled_results_4C_37C_20190915.csv")
data_table_rmsd <- read.csv("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/results rmsd analysis 20190915 data/pooled_results_4C_37C_20190915_trimmed.csv")


## SUPERPLOTS: explore how the different repeats (or temperature) affect the data distribution for a-syn and VAMP-2

# Transform into master data frame 
df_all <- select(data_table_rmsd,Temperature,Repeat,condition,mCLING,asyn, VAMP2) # include repeat column

# tidy up master data frame
df_all_melt = melt(df_all)

# Order the different conditions to make it easier to display control vs treatments
df_all_melt$condition <- factor(df_all_melt$condition, levels = c("phys","egta","egtak"))

# A-syn
{
  
  df_asyn <- subset(df_all_melt,variable == 'asyn')
  
  # Superplot 1: a-syn RMSD as a function of calcium
  ggplot(df_asyn,aes(x = condition, y = value, color = factor(Repeat))) + 
    geom_quasirandom(alpha = 0.4, size = 2.0,dodge.width=0.8)+
    scale_shape_manual(values = c(2,0,21))+
    geom_boxplot(alpha = 0.2,outlier.shape = NA,fatten = NULL)+
    xlab("Calcium treatment")+
    ylab("RMS distance (nm) from centroid") + 
    prism_multi()+
    stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.75, size = 1,linetype = "solid", position = position_dodge2(width = 0.8))
  ggsave("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/rmsd SUPERPLOTS/RMSD_asyn_repeat_temperature.pdf",width = 110, height = 80, units = "mm" )
  
  # Superplot 2: Group by condition, dodge scatterplots by repeat, and colour them by temperature
  ggplot(df_asyn,aes(x = condition, y = value, group=factor(Repeat),color = Temperature ))  + 
    geom_quasirandom(alpha = 0.4, size = 2.0,dodge.width=0.8)+
    #scale_shape_manual(values = c(2,0,21))+
    #geom_boxplot(alpha = 0.2,outlier.shape = NA,fatten = NULL)+
    xlab("Calcium treatment")+
    ylab("RMS distance (nm) from centroid") + 
    prism_multi()+
    stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.75, size = 1,linetype = "solid", position = position_dodge2(width = 0.8))
  ggsave("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/rmsd SUPERPLOTS/RMSD_asyn_repeat_temperature2.pdf",width = 110, height = 80, units = "mm" )
  
  # Superplot 3:  Group by condition, dodge and colour by Temperature
  ggplot(df_asyn,aes(x = condition, y = value,color = factor(Temperature))) +
    geom_quasirandom(alpha = 0.4, size = 2.0,dodge.width = 0.8)+
    #geom_boxplot(alpha = 0.6,outlier.shape = NA,fatten=NULL)+
    xlab("Conditions")+
    ylab("Cluster size")+
    prism_multi()+
    stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.75, size = 1,linetype = "solid", position = position_dodge2(width = 0.8))
  
  # t test between the phys-egta conditions 
  x_rmsd <- subset(df_asyn,condition == 'phys')
  y_rmsd<- subset(df_asyn,condition == 'egta')
  x = x_rmsd %>% select(value)
  y = y_rmsd %>% select(value)
  t.test(x,y) # 0.05486

}

# VAMP-2
{
  df_VAMP2 <- subset(df_all_melt,variable == 'VAMP2')
  
  # Superplot 1: VAMP2 RMSD as a function of calcium
  ggplot(df_VAMP2,aes(x = condition, y = value, color = factor(Repeat))) + 
    geom_quasirandom(alpha = 0.4, size = 2.0,dodge.width=0.8)+
    scale_shape_manual(values = c(2,0,21))+
    geom_boxplot(alpha = 0.2,outlier.shape = NA,fatten = NULL)+
    xlab("Calcium treatment")+
    ylab("RMS distance (nm) from centroid") + 
    prism_multi()+
    stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.75, size = 1,linetype = "solid", position = position_dodge2(width = 0.8))
  ggsave("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/rmsd SUPERPLOTS/RMSD_VAMP2_repeat_temperature.pdf",width = 110, height = 80, units = "mm" )
  
  # Superplot 2: Group by condition, dodge scatterplots by repeat, and colour them by temperature
  ggplot(df_VAMP2,aes(x = condition, y = value, group=factor(Repeat),color = Temperature )) + 
    geom_quasirandom(alpha = 0.4, size = 2.0,dodge.width=0.8)+
    #geom_boxplot(alpha = 0.2,outlier.shape = NA,fatten = NULL)+
    xlab("Calcium treatment")+
    ylab("RMS distance (nm) from centroid") + 
    prism_multi()+
    stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.75, size = 1,linetype = "solid", position = position_dodge2(width = 0.8))
  ggsave("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/rmsd SUPERPLOTS/RMSD_VAMP2_repeat_temperature2.pdf",width = 110, height = 80, units = "mm" )
  
  # Superplot 3:  Group by condition, dodge and colour by Temperature
  ggplot(df_VAMP2,aes(x = condition, y = value,color = factor(Temperature))) +
    geom_quasirandom(alpha = 0.4, size = 2.0,dodge.width = 0.8)+
    #geom_boxplot(alpha = 0.6,outlier.shape = NA,fatten=NULL)+
    xlab("Conditions")+
    ylab("Cluster size")+
    prism_multi()+
    stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.75, size = 1,linetype = "solid", position = position_dodge2(width = 0.8))
  
  # t test between phys-egta
  x_rmsd <- subset(df_VAMP2,condition == 'phys')
  y_rmsd<- subset(df_VAMP2,condition == 'egta')
  x = x_rmsd %>% select(value)
  y = y_rmsd %>% select(value)
  t.test(x,y) # 0.01407


}

  
# ESTIMATION GRAPHICS to visualize the change in RMSD between treatments 
# source: https://cran.r-project.org/web/packages/dabestr/vignettes/using-dabestr.html

{
  # visualize effect size instead of using p values
  
  library(dabestr)      # for plotting the differences between distributions as an estimation graphic
  
  # ASYN
  unpaired_mean_diff_asynRMSD <- dabest(df_asyn,condition,value,
                                   idx = c("phys","egta","egtak"), # add the egta K for a multi group plot
                                   paired = FALSE)
  
  plot(unpaired_mean_diff_asynRMSD, rawplot.ylabel = "RMSD value",rawplot.ylim = c(0, 400),palette = 'Set2')
  
  # for figure
  plot(unpaired_mean_diff_asynRMSD,
       rawplot.ylabel = "RMSD (nm)",
       rawplot.markersize = mk_size,
       palette = "Set2",
       tick.fontsize  = 8, axes.title.fontsize = 8,effsize.markersize = 2)
  #swarmplot.params = list(shape = c(1) ) )
  
  ggsave("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/results cluster size/clusterSize_asyn_RMSDAll.pdf",width = 80, height = 80, units = "mm")
  
  
  
  # VAMP2
  unpaired_mean_diff_vampRMSD <- dabest(df_VAMP2,condition,value,
                                    idx = c("phys","egta","egtak"), # add the egta K for a multi group plot
                                    paired = FALSE)
  
  plot(unpaired_mean_diff_vampRMSD, rawplot.ylabel = "RMSD value",rawplot.ylim = c(0, 400),palette = 'Set2')
  
  # for figure
  plot(unpaired_mean_diff_vampRMSD,
       rawplot.ylabel = "RMSD (nm)",
       rawplot.markersize = mk_size,
       palette = "Accent",
       tick.fontsize  = 8, axes.title.fontsize = 8,effsize.markersize = 2)
  #swarmplot.params = list(shape = c(1) ) )
  
  ggsave("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/results cluster size/clusterSize_VAMP2_RMSDAll.pdf",width = 80, height = 80, units = "mm")
  
  
  
  
  # plot them together in a multi-group plot
  multi.group <- 
    df_all_melt%>%
    dabest(condition, value, 
           idx = list(c("phys", "egta"), 
                      c("phys", "egta")),
           paired = FALSE
    )
  
  multi.group 
  plot(multi.group,rawplot.ylabel = "RMSD value")
  # Calling the object automatically prints out a summary.
}

## Deprecated code

# select data columns
df     <- select(data_table_rmsd,Temperature,condition,mCLING,asyn, VAMP2)

# Organize into tidy form
df_melt = melt(df)

# remove outliers (RMSD values of zero)
df_melt <- filter(df_melt, value > 0)

# Plotting all conditions together for thesis and paper 
{
  # Using ggplot2, generate dot plot based on calcium condition for all conditions
  ggplot(df_melt,aes(x = condition, y = value, color = variable,fill = variable, shape = variable)) + 
    geom_quasirandom(alpha = 0.4, size = 2.0,dodge.width = 0.8)+
    scale_shape_manual(values = c(2,0,21))+
    #geom_violin(alpha = 0.3)+
    geom_boxplot(alpha = 0.2,outlier.shape = NA)+
    #scale_color_manual(values = c("blue","green4"))+
    xlab("Calcium treatment")+
    ylab("RMS distance (nm) from centroid") + 
    prism_multi()
  ggsave("E:\\OneDrive\\OneDrive - University Of Cambridge\\lag\\microscopy work\\synaptosomes\\manuscript figures\\RMSD_summary.pdf",width = 110, height = 80, units = "mm")
  ggsave("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/rmsd analysis 20190807 summary/RMSD_summary_new.pdf",width = 150, height = 120, units = "mm" )
  
  
  # Using ggplot2, generate dot plot based on temperature
  ggplot(df_melt,aes(x = Temperature, y = value, color = variable,fill = variable)) + 
    geom_quasirandom(alpha = 0.4, size = 0.6,dodge.width = 0.8)+
    #geom_violin(alpha = 0.3)+
    geom_boxplot(alpha = 0.2)+
    #scale_color_manual(values = c("blue","green4"))+
    xlab("Compared data sets")+
    ylab("RMSD") + 
    prism_multi()
  ggsave("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/rmsd analysis 20190807 summary/RMSD_summary_byTemp.pdf",width = 150, height = 120, units = "mm" )
}

##----------Statistics----------##
{
##------One-way ANOVA------##


# Use a one-way ANOVA to test the differences among the means of the RMSD values for each channel in the three different calcium conditions
df_rc = subset(df_melt, variable == 'mCLING')
df_gc = subset(df_melt, variable == 'asyn')
df_bc = subset(df_melt, variable == 'VAMP2')

res.aov <- aov(value ~ condition, data = df_gc)
summary(res.aov) # obtain summary from one-way ANOVA
TukeyHSD(res.aov)# interpret F value using the Tukey metric

#1. Check homogeneity of variances 
plot(res.aov,1)
# Use Levene's test to check the homogeneity of the variances in the data
leveneTest(value ~ condition, data = df_gc)

#2. Check normality condition 
plot(res.aov,2)

##------Kruskal-Wallis------##

# if the data do not satisfy the conditions for a one-way ANOVA, then we can try the non-parametric alternative (Kruskal-Wallis test)
# along with a wilcoxon multiple pairwise comparison test to get specific p-values among conditions. 

kruskal.test(value ~ condition, data = df_gc) # if the p-value resulting from this test is <0.05, then we can conclude there are significant differences between groups

# We can calculate p-values between group levels with corrections for multiple testing 
pairwise.wilcox.test(df_gc$value, df_gc$condition,
                     p.adjust.method = "BH")

kruskal.test(value ~ condition, data = df_bc) # if the p-value resulting from this test is <0.05, then we can conclude there are significant differences between groups

# We can calculate p-values between group levels with corrections for multiple testing 
pairwise.wilcox.test(df_bc$value, df_bc$condition,
                     p.adjust.method = "BH")



# get the means of each condition
mean_rc_phys = subset(df_rc, condition == 'phys')
mean(mean_rc_phys[,"value"])
mean_rc_egta = subset(df_rc, condition == 'egta')
mean(mean_rc_egta[,"value"])
mean_rc_egtak = subset(df_rc, condition == 'egtak')
mean(mean_rc_egtak[,"value"])

mean_gc_phys = subset(df_gc, condition == 'phys')
mean(mean_gc_phys[,"value"])
mean_gc_egta = subset(df_gc, condition == 'egta')
mean(mean_gc_egta[,"value"])
mean_gc_egtak = subset(df_gc, condition == 'egtak')
mean(mean_gc_egtak[,"value"])

mean_bc_phys = subset(df_bc, condition == 'phys')
mean(mean_bc_phys[,"value"])
mean_bc_egta = subset(df_bc, condition == 'egta')
mean(mean_bc_egta[,"value"])
mean_bc_egtak = subset(df_bc, condition == 'egtak')
mean(mean_bc_egtak[,"value"])

}



# from the superplots paper tutorial (Lord, Velle, et al. 2019)
{
ReplicateAverages <- df_asyn %>% group_by(condition, Repeat) %>%
  summarise_all(list(mean))

ggplot(df_asyn,aes(x = condition, y = value, color = factor(Repeat))) + 
  geom_beeswarm(cex=3) + 
  scale_colour_brewer(palette = "Set1") + 
  geom_beeswarm(data=ReplicateAverages, size=8) + 
  stat_compare_means(data=df_asyn, 
                     comparisons = list(c("A", "B")), 
                     method="t.test", paired=TRUE) 
}

