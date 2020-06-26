
#----Supplementary code accompanying Part2 Chapter2 of PVR thesis----#

# Script to plot the FRC value as a function of projection number for an OPT acquisition

# The input is a csv file with obtained after splitting an OPT data set, and calculating the FRC value from two data sets using NanoJ Squirrell.


# Pedro Vallejo Ramirez
# Laser Analytics Group, University of Cambridge
# Created: 13/08/2019


library(tidyverse)
library(ggplot2)
library(ggbeeswarm)
# For a full page a4 figure, w = 160 mm, h = ?
# For a half-page a4 figure, w = 80 mm, h = 70 mm (no legend)
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
    theme(legend.title =element_text(size=font_size, color=color.axis.title))+
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

# import data
datatab <- read.csv("E:\\OneDrive\\OneDrive - University Of Cambridge\\lag\\opt\\frc study\\results_overall_frcMeanvsframe.csv")
# tidy up into data frame
df <- data.frame(projections = datatab[,1],
                 value = datatab[,2]/1000)

ggplot(df, aes(x = projections,y = value))+
  geom_quasirandom(alpha = 0.8,size = 0.3)+
  geom_boxplot(aes(group = projections),alpha = 0.3,outlier.shape = NA)+
  labs(x='Number of projections',y='FRC value (um)')+
  scale_x_continuous(breaks=c(0,32,64,128,256))+
  prism_multi()
ggsave("E:\\OneDrive\\OneDrive - University Of Cambridge\\lag\\opt\\frc study\\frc_v_projections.pdf",width = 110, height = 80, units = "mm")
