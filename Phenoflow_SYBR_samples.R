#Phenoflow_package tutorial
#For SYBR stained samples. 
#Created: October 23, 2018
#Last run: October 28, 2018

#### Install packages ####
#install.packages("devtools")
#install_github("lievenclement/flowFDA", build_vignettes=TRUE)
#devtools::install_github("rprops/Phenoflow_package")
#install.packages("hexbin")

library(devtools)
library("Phenoflow") # for fingerprinting
library("flowViz") # for plotting
library("ggplot2") # for plotting
library("flowAI") # for denoising

set.seed(777)

# Clear all existing data
rm(list=ls())

# Set working directory 
setwd("/Users/kathrynschmidt/Downloads/FCS.files.Chlorella.Microcystis.Experiment/Ruben_visit_10.26.18/")

#### Load data ####
#data(flowData)
#head(flowData[[1]]@exprs)
path = "SYBR_stained"
algaeData <- read.flowSet(path = path, emptyValue = FALSE, transformation = FALSE, pattern=".fcs")

#### Denoise data ####
# Next: you need to select a biological gate
# transform data (data that is exported is not transformed)
# Select phenotypic features of interest and transform parameters, arcsin hyperbolic transformation
algaeData_transformed <- transform(algaeData,`BL1-A`=asinh(`BL1-A`), 
                                   `SSC-A`=asinh(`SSC-A`), 
                                   `BL3-A`=asinh(`BL3-A`), 
                                   `VL1-A`=asinh(`VL1-A`), 
                                   `VL4-A`=asinh(`VL4-A`), 
                                   `FSC-A`=asinh(`FSC-A`))
param=c("BL1-A", "FL3-H","SSC-A","BL3-A","VL1-A","VL4-A","FSC-A") #only used these parameters for analysis
remove(algaeData) #to make sure you are looking at the transformed data

### Create a PolygonGate for denoising the dataset
### Define coordinates for gate in sqrcut1 in format: c(x,x,x,x,y,y,y,y)
sqrcut1 <- matrix(c(8.75,8.75,14,14,
                    3,7.5,14,3),ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL1-H","FL3-H")
polyGate1 <- polygonGate(.gate=sqrcut1, filterId = "Total Cells")

###  Gating quality check
sampleNames(flowData_transformed)
xyplot(`VL4-A` ~ `BL1-A`, data=algaeData_transformed[sample(10:11)],
       scales=list(y=list(limits=c(0,14)),
                   x=list(limits=c(6,16))),
       axis = axis.default, nbin=125, 
       par.strip.text=list(col="white", font=2, cex=2), smooth=FALSE)
#filter=polyGate1,
#can input "data=flowData_transformed[sample(1:41)[1:6]]" which will randomly select 6 to visualize out of the 41 experiments

### Isolate only the cellular information based on the polyGate1
flowData_transformed <- Subset(flowData_transformed, polyGate1)
