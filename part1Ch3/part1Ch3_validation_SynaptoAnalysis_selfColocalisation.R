# Supplementary code accompanying Chapter 3 of PVR thesis to generate the plot in Fig. 3.14

# Script for plotting results from the self-colocalisation tests on the synaptosome analysis scripts. 

# A localisation data table was split in even and odd halves, and then were used as different channels in the synaptosome
# analysis to check how the data colocalises with itself,

# The input is a csv file with the results from the overlap calculations from Synapto-analysis.  

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

# Import data
data_table_selfColoc <- read.csv("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/figure self colocalisation/results_summary_updated.csv")

# organize into data frame with separate columns for the WOC and Manders M1 coefficient on each data set
df_selfColoc <- data.frame(channels = c(rep('even-even',nrow(data_table_selfColoc)), 
                              rep('even-odd',nrow(data_table_selfColoc)),
                              rep('even-Rotated',nrow(data_table_selfColoc))
                              ), 
                 WOC = c(data_table_selfColoc[,5],
                           data_table_selfColoc[,6],
                           data_table_selfColoc[,7]),
                 M1 = c(data_table_selfColoc[,8],
                        data_table_selfColoc[,9],
                        data_table_selfColoc[,11]))

df_selfColoc_melt <- melt(df_selfColoc)
  
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
    #theme(legend.position="bottom")+
    
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

# Using ggplot2, generate dot plot
ggplot(df_selfColoc_melt,aes(x = channels, y = value,color = variable,shape = variable)) + 
  geom_quasirandom(size=0.7, alpha = 1,dodge.width = 0.8)+
  #scale_color_manual(values = c("blue","green4"))+
  xlab("Compared data sets")+
  ylab("Spatial co-ocurrence")+
  #ylim(0,1) + 
  prism_multi()+
  scale_color_brewer(palette = "Set2")+
  scale_fill_brewer(palette = "Set2")+
  theme(legend.position = "none")

# save to drive
ggsave("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/figure self colocalisation/self_colocalisationPlot.pdf",width = 105, height = 60, units = "mm")
