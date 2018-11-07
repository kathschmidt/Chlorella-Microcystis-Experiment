#Phenoflow_package tutorial
#Ruben Props; edited by Kathryn Schmidt
#Experiment: Chlorella-Microcystis Invasion
#Date: October 23, 2018
#LAST RUN: November 6, 2018

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

#### Load data ####
data(flowData)
head(flowData[[1]]@exprs)
#path = "test_data"
#flowData <- read.flowSet(path = path, transformation = FALSE, pattern=".fcs")

#### Denoise data ####
# Next: you need to select a biological gate
# transform data (data that is exported is not transformed)
# Select phenotypic features of interest and transform parameters, arcsin hyperbolic transformation
flowData_transformed <- transform(flowData,`FL1-H`=asinh(`FL1-H`), 
                                  `SSC-H`=asinh(`SSC-H`), 
                                  `FL3-H`=asinh(`FL3-H`), 
                                  `FSC-H`=asinh(`FSC-H`))
param=c("FL1-H", "FL3-H","SSC-H","FSC-H") #only used these four parameters for analysis, for us use area "A"
remove(flowData) #to make sure you are looking at the transformed data

### Create a PolygonGate for denoising the dataset
### Define coordinates for gate in sqrcut1 in format: c(x,x,x,x,y,y,y,y)
sqrcut1 <- matrix(c(8.75,8.75,14,14,
                    3,7.5,14,3),ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL1-H","FL3-H")
polyGate1 <- polygonGate(.gate=sqrcut1, filterId = "Total Cells")

###  Gating quality check
sampleNames(flowData_transformed)
xyplot(`FL3-H` ~ `FL1-H`, data=flowData_transformed[1:6], filter=polyGate1,
       scales=list(y=list(limits=c(0,14)),
                   x=list(limits=c(6,16))),
       axis = axis.default, nbin=125, 
       par.strip.text=list(col="white", font=2, cex=2), smooth=FALSE)
      #can input "data=flowData_transformed[sample(1:41)[1:6]]" which will randomly select 6 to visualize out of the 41 experiments

### Isolate only the cellular information based on the polyGate1
flowData_transformed <- Subset(flowData_transformed, polyGate1)

### Extract metadata from sample names (metadata will be different per. project)
metadata <- data.frame(do.call(rbind, lapply(strsplit(flowCore::sampleNames(flowData),"_"), rbind)))
colnames(metadata) <- c("Cycle_nr", "Location", "day", "timepoint", "Staining", "Reactor_phase", "replicate")

summary <- fsApply(x = flowData_transformed, FUN = function(x) apply(x, 2, max), use.exprs = TRUE)
head(summary) #finding the maximum
maxval <- max(summary[,"FL1-H"]) #replace with the column representing the green fluorescence channel (e.g. "FITC-H"), take maximum value for parameter of interest
mytrans <- function(x) x/maxval #normalize from 0-1 with the maximum data
flowData_transformed <- transform(flowData_transformed,`FL1-H`=mytrans(`FL1-H`),
                                  `FL3-H`=mytrans(`FL3-H`), 
                                  `SSC-H`=mytrans(`SSC-H`),
                                  `FSC-H`=mytrans(`FSC-H`))
###### Fingerprinting ######
###density estimation of microbial community, wil do for all bivaried parameters
### Randomly resample to the lowest sample size, use if samples have below 10,000 cells
#flowData_transformed <- FCS_resample(flowData_transformed, replace=TRUE)
### Calculate fingerprint with bw = 0.01, nbin 128 is the ideal trade-off between speed and accuracy
### Check wikipedia page for explaination 
fbasis <- flowBasis(flowData_transformed, param, nbin=128, 
                    bw=0.01,normalize=function(x) x)

### Calculate Diversity from normalized fingerprint (d can be increased, but not lowered)
Diversity.fbasis <- Diversity(fbasis,d=3,plot=TRUE, R=999)

# Diversity assessment with cleaning --> use this function if possible for final result, takes in account sample size
Diversity.clean <- Diversity_rf(flowData_transformed, param = param, R = 3, R.b = 3,
                                cleanFCS = FALSE)
p2 <- ggplot(data = Diversity.clean, aes(x = as.numeric(as.character(metadata$day)), y = D2, color = metadata$Reactor_phase))+
  geom_point(size = 8, alpha = 0.7)+
  geom_line()+
  scale_color_manual(values = c("#a65628", "red", 
                                "#ffae19", "#4daf4a", "#1919ff", "darkorchid3", "magenta"))+
  theme_bw()+
  labs(color = "Reactor phase", y = "Phenotypic diversity (D2)", x = "Days")+
  geom_errorbar(aes(ymin=D2-sd.D2, ymax=D2+sd.D2), width=0.05, color="black")
print(p2)

### Export diversity estimates to .csv file in the chosen directory
#write.csv2(file="results.metrics.csv", Diversity.clean)

###### BETA DIVERSITY #####
# Beta-diversity assessment of fingerprint
beta.div <- beta_div_fcm(fbasis, ord.type="PCoA")

# Plot ordination
plot_beta_fcm(beta.div, color = metadata$Reactor_phase, labels="Reactor phase") + 
  theme_bw() +
  geom_point(size = 8, alpha = 0.5)

#### Extract cell counts ####
### Creating a rectangle gate for counting HNA and LNA cells
rGate_HNA <- rectangleGate("FL1-H"=c(asinh(20000), 20)/maxval,"FL3-H"=c(0,20)/maxval, 
                           filterId = "HNA bacteria")
### Normalize total cell gate
sqrcut1 <- matrix(c(8.75,8.75,14,14,3,7.5,14,3)/maxval,ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL1-H","FL3-H")
polyGate1 <- polygonGate(.gate=sqrcut1, filterId = "Total Cells")

### Check if rectangle gate is correct, if not, adjust rGate_HNA
xyplot(`FL3-H` ~ `FL1-H`, data=flowData_transformed[1], filter=rGate_HNA,
       scales=list(y=list(limits=c(0,1)),
                   x=list(limits=c(0.4,1))),
       axis = axis.default, nbin=125, par.strip.text=list(col="white", font=2, 
                                                          cex=2), smooth=FALSE)
### Extract the cell counts
a <- flowCore::filter(flowData_transformed, rGate_HNA) 
HNACount <- summary(a);HNACount <- toTable(HNACount)
s <- flowCore::filter(flowData_transformed, polyGate1)
TotalCount <- summary(s);TotalCount <- toTable(TotalCount)

### Extract the volume
vol <- as.numeric(flowCore::fsApply(flowData_transformed, FUN = function(x) x@description$`$VOL`))/1000

### Store the data
results_counts <- data.frame(Samples=flowCore::sampleNames(flowData_transformed), 
                             Total.cells = TotalCount$true/vol, HNA.cells = HNACount$true/vol)

### Exporting cell counts to .csv file to working directory
write.csv2(file="results.counts.csv", results_counts)

### Plot cell density
ggplot(data = results_counts, aes(x = as.numeric(as.character(metadata$day)), y = Total.cells, color = metadata$Reactor_phase))+
  geom_point(size = 8, alpha = 0.9)+
  scale_color_manual(values = c("#a65628", "red", 
                                "#ffae19", "#4daf4a", "#1919ff", "darkorchid3", "magenta"))+
  theme_bw()+
  labs(color = "Reactor phase", y = "Total cell density (cells/µL)", x = "Days")  

?RandomF_FCS #will use machine learning to discriminate between two populations
