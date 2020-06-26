#----Supplementary code accompanying Part1 Chapter2 of PVR thesis----#

# The input is a csv file with the values from a line profile taken in ImageJ/Fiji.
# The script fits a Gaussian to the data and plots it

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
data_table    <- read.csv("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/figure synapto simulation cluster size/crossSection_blob_GC_equal_MeasuredAgain.csv")
data_table_RC <- read.csv("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/figure synapto simulation cluster size/crossSection_blob_RC_equal.csv")
data_table_BC <- read.csv("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/figure synapto simulation cluster size/crossSection_blob_BC_equal.csv")


# Fit Gaussian to the data
fitG =
  function(x,y,mu,sig,scale){
    
    f = function(p){
      d = p[3]*dnorm(x,mean=p[1],sd=p[2])
      sum((d-y)^2)
    }
    
    optim(c(mu,sig,scale),f)
  }

# Use initial parameters from the fit in Fiji (mu = c, sigma = d in the Fiji gaussian fit ) 
fitGC = fitG(data_table$X,data_table$Y,0.62,.13,1)
p = fitGC$par
fwhm_GC = 2.3548*p[2]

# for RC , the Gaussian fit in ImageJ yields c = 0.64 and d = 0.12
# for BC , the Gaussian fit in ImageJ yields c = 0.6 and d = 0.1
fitRC   = fitG(data_table_RC$X,data_table_RC$Y,0.6,.1,1)
p_RC    = fitRC$par
fwhm_RC = 2.3548*p_RC[2]


fitBC   = fitG(data_table_BC$X,data_table_BC$Y,0.6,.1,1)
p_BC    = fitBC$par
fwhm_BC = 2.3548*p_BC[2]


# Tidy up data
df <- data.frame(x = rep(data_table$X,3),
                 y = c(data_table$Y,data_table_RC$Y,data_table_BC$Y))

# data frame with y values for each channel
df_Gauss <- data.frame(x = rep(data_table$X,3),
                 GC = data_table$Y,
                 RC = data_table_RC$Y,
                 BC = data_table_BC$Y)
df_Gauss_melt = melt(df_Gauss)

point_sz = 0.5
lwd = 1.7

ggplot(df_Gauss,aes(x= x,y = GC))+
  geom_point(size = point_sz, aes(y=GC),colour = "#00BA38") +
  geom_point(size = point_sz, aes(y=RC),colour = "#F6766D") +
  geom_point(size = point_sz, aes(y=BC),colour = "#619CFF") +
  
  geom_line(aes(x = x,y = p[3]*dnorm(x,p[1],p[2]) ), # mu is p[1] and sigma is p[2]
            colour = "#00BA38", alpha = 0.5, size = lwd)+
  geom_line(aes(x = x,y = p_RC[3]*dnorm(x,p_RC[1],p_RC[2]) ), # mu is p[1] and sigma is p[2]
            colour = "#F6766D", alpha = 0.5, size = lwd)+
  geom_line(aes(x = x,y = p_BC[3]*dnorm(x,p_BC[1],p_BC[2]) ), # mu is p[1] and sigma is p[2]
            colour = "#619CFF", alpha = 0.5, size = lwd)+
  xlab("Distance (um)")+
  ylab("Localisation density")+
  prism_multi()


ggsave("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/figure synapto simulation cluster size/clustersize_GaussianFWHM.pdf",width = 100, height = 60, units = "mm")

