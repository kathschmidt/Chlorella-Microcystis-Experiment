## Merging metadata from FCM and experiment
## Kathryn Schmidt
## Date: November 12, 2018

## Input files: Well_IDs, Microcystis_Chlorella_FCM_Experiment
## Output file: MC_metadata

# Clear all existing data
rm(list=ls())

# Close graphics devices
graphics.off()

# set working directory
setwd("/Volumes/GoogleDrive/My Drive/Denef Lab/")

# read in data
Well_IDS<-read.csv("Fall 2017/Experiment (1.2.18)/Well_IDs.11.12.18.csv")
exp_metadata<-read.csv("Fall 2018/FCM with experimental samples/Microcystis_Chlorella_FCM_Experiment.11.12.18.csv")

# merge files
MC_metadata<- cbind(Well_IDS,exp_metadata)

# export csv
write.csv(MC_metadata, file="MC_metadata.csv")
