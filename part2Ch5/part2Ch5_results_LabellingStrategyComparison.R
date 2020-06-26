
#----Supplementary code accompanying Part2 Chapter2 of PVR thesis----#

# Script to plot a line profile through different regions of slice 402 of the right medial lung lobe reconstruction (OptiJ)

# The inputs are csv files from line profiles drawn in ImageJ/Fiji

# Pedro Vallejo Ramirez
# Laser Analytics Group, University of Cambridge
# Created: 13/08/2019


library(tidyverse)
library(ggplot2)

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


#data_lineProfile <- read.csv("E:\\OneDrive\\OneDrive - University Of Cambridge\\lag\\opt\\Lung data\\quantification of autofluorescence 20190813\\LineProfile_antiSurfactant_contrastRatio.csv")
data_lineProfile <- read.csv("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/opt/Lung data/quantification of autofluorescence 20190813/LineProfile_antiSurfactant_contrastRatio.csv")
data_lineProfile <- read.csv("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/opt/Lung data/quantification of autofluorescence 20190813/LineProfile_rightmedial_contrastRatio.csv")

df_lineProfile <- data.frame(x = data_lineProfile[,1]*12.90,
                             value = data_lineProfile[,2]-450)

# subtract 300 from antisurfactant
#df_lineprof <- subset(df_lineProfile, value >= 5)

ggplot(df_lineProfile,aes(x = x, y = value))+
  geom_line()+
  labs(x='Distance (um)',y='Fluorescence intensity (a.u.)')+
  prism_multi()
ggsave("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/opt/Lung data/quantification of autofluorescence 20190813/LineProfilePlot_rmb.pdf",width = 86, height = 62.5, units = "mm")
