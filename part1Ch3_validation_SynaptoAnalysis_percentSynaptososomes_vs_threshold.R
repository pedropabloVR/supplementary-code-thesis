# Supplementary code accompanying Part1 Chapter2 of PVR thesis.

# plotting % of synaptosomes detected vs overlap thresholds from SynaptoAnalysis results

# Pedro Vallejo Ramirez
# Laser Analytics Group, University of Cambridge
# Created: 2020-02-20


# Imports a .csv output by a matlab 

library(tidyverse)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
theme_update(plot.title = element_text(hjust = 0.5))

prism_multi <- function() {
  
  # Generate the colors for the chart procedurally with RColorBrewer
  palette <- brewer.pal("Greys", n=9)
  color.background = "white"
  color.grid.major = "white"
  color.axis.text = "black"
  color.axis.title = "black"
  color.title = palette[9]
  font_size <- 10
  title_size <- 10
  
  # Begin construction of chart
  theme_bw(base_size=11) +
    
    # Format the legend, uncomment position if need be
    theme(legend.text  = element_text(size=font_size, color=color.axis.title))+
    theme(legend.title = element_text(size=font_size, color=color.axis.title))+
    
    # Set title and axis labels, and format these and tick marks
    theme(plot.title=element_text(color=color.title, size=font_size, vjust=1.25, hjust = 0.4, face = "bold")) +
    theme(axis.text.x=element_text(size=font_size,color=color.axis.text, margin = unit(c(0.1, 0.5, 0.3, 0.5), "cm"))) +
    theme(axis.text.y=element_text(size=font_size,color=color.axis.text, margin = unit(c(0.5, 0.1, 0.5, 0.5), "cm"))) +
    theme(axis.title.x=element_text(size=font_size,color=color.axis.title, 
                                    vjust=1.25, margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))) +
    theme(axis.title.y=element_text(size=font_size,color=color.axis.title, 
                                    vjust=-4, hjust = 0.5, margin = unit(c(0.2, 0.2, 0.2, 
                                                                           0.2), "cm"))) 
  #theme(axis.text.x = element_text(angle = 45, hjust = 0.1, vjust = 0.10)) 
} # plotting function


dataset <-'37Cb'
#data <- read.csv(paste("F:\\synaptosomes\\analysis\\",dataset,"\\Results_20181205\\Results_combined\\detection_vs_threshold_results.csv",sep=""))

data <- read.csv(paste("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/figure detectedSynaptosomes_vs_threshold plots/detection_vs_threshold_results_",dataset,".csv",sep=""))

df <- data.frame(threshold=data[,1],
                 Synaptosomes_detected=c(data[,2], data[,3]),
                 Type=c(rep('R-G overlap threshold',nrow(data)),rep('R-B overlap threshold',nrow(data)))
)

g<-ggplot(df,aes(x=threshold,y=Synaptosomes_detected,color=Type))+geom_line(size = 0.5)+xlab("Threshold value")+ylab("% synaptosomes detected")+labs(title=dataset)+prism_multi()+theme(legend.position = "bottom")
#g+ scale_x_continuous(breaks = seq(0, 100, by = 10),expand = c(0, 0))#scale_y_continuous(limits = c(0, 100),breaks = seq(0, 100, by = 10),expand = c(0, 0))


#ggsave(paste("F:\\synaptosomes\\analysis\\",dataset,"\\Results_20181205\\Results_combined\\detectedSynaptosomes_vs_threshold_",dataset,".pdf",sep=""), width = 75, height = 55, units = "mm")
ggsave(paste("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/detectedSynaptosomes_vs_threshold plots/detection_vs_threshold_results_new",dataset,".pdf",sep=""),width = 75, height = 55, units = "mm" )
