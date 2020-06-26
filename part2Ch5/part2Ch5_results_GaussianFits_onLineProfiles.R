
#----Supplementary code accompanying Part2 Chapter2 of PVR thesis----#

# Script for plotting a Gaussian fit on a line profile made in Fiji on slice 499 of the reconstructed right medial lobe
# A gaussian is also fitted to a line profile through the image of a bead reconstructed using OptiJ

# Pedro Vallejo Ramirez
# Laser Analytics Group, University of Cambridge
# Created: 29/03/2019


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
} # plotting function


# Import data tables post-overlap threshold calculation and post gaussian fit in Fiji
#data_table_1    <- read.csv("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/opt/frc study/crossSection_ROI_Slice499_longer_ROI.csv")
#data_table_2    <- read.csv("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/opt/frc study/crossSection_ROI_2_Slice499.csv")
data_table_1 <- read.csv("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/opt/frc study/gaussian profile from slice 499/crossSection_ROI_Slice499_longer_ROI.csv")
data_table_2 <- read.csv("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/opt/frc study/gaussian profile from slice 499/crossSection_ROI_2_Slice499.csv")

data_table_1_fit <- read.csv("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/opt/frc study/gaussian profile from slice 499/crossSection_ROI_Slice499_GaussianFit.csv")
data_table_2_fit <- read.csv("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/opt/frc study/gaussian profile from slice 499/crossSection_ROI_2_Slice499_GaussianFit.csv")

foox <- data_table_1_fit$X1[!is.na(data_table_1_fit$X1)]
fooy <- data_table_1_fit$Y1[!is.na(data_table_1_fit$Y1)]

foox2 <- data_table_2_fit$X1[!is.na(data_table_2_fit$X1)]
fooy2 <- data_table_2_fit$Y1[!is.na(data_table_2_fit$Y1)]

pdf("rplot.pdf") 
plot(foox,fooy, xlim = c(0,120),ylim = c(150,900),col = "#00AFBB",pch = 1, xlab="Distance (um)",ylab = "Fluorescence intensity (a.u.)")
par(new=T)
plot(foox2,fooy2, xlim = c(0,120),ylim = c(150,900), col = "#FC4E07",pch = 2, xlab="Distance (um)",ylab = "Fluorescence intensity (a.u.)")

lines(data_table_1_fit$X0,data_table_1_fit$Y0, lwd = 3.0,col = "#00AFBB", xlim = c(0,120),ylim = c(150,900),xlab="",ylab = "")
lines(data_table_2_fit$X0,data_table_2_fit$Y0, lwd = 3.0,col = "#FC4E07",xlim = c(0,120),ylim = c(150,900),xlab="",ylab = "")

# 80x60 mm = 3.15 x 2.36 in
# 120x100 mm = 4.72x3.94 in 
# if the landscape mark is ticked, then invert the dimensions (small number x large number)
dev.off()

plot(foox2,fooy2, xlim = c(0,120),ylim = c(150,800),xlab="Distance (um)",ylab = "Fluorescence intensity (a.u.)")
par(new=T)
lines(data_table_2_fit$X0,data_table_2_fit$Y0, xlim = c(0,120),ylim = c(150,900),xlab="",ylab = "",axes = F)

# Save them via the Export tab in the gui window. 

# repeating procedure with the bead data
data_table_bead_fit <- read.csv("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/opt/frc study/bead_data/GaussianFit_Values_ROI2.csv")
foox3 <- data_table_bead_fit$X1[!is.na(data_table_bead_fit$X1)]
fooy3 <- data_table_bead_fit$Y1[!is.na(data_table_bead_fit$Y1)]

pdf("rplot.pdf") 
plot(foox3,fooy3, xlim = c(0,180),ylim = c(200,1700),xlab="Distance (um)",ylab = "Fluorescence intensity (a.u.)")
par(new=T)
lines(data_table_bead_fit$X0,data_table_bead_fit$Y0, xlim = c(0,180),ylim = c(0,1700),xlab="",ylab = "")


