################################################################################
## This script is to generate a null distribution to test confidence of unique segments
## Started: 12/02/25
## Updated: 20/10/25
##
## I used all single-cell data processed in this work for generating the null distribution
## For each dataset/cohort/study, I first run the "sc_confidence_scores.R" script
## This script outputs a file per sample with segTables of each cell in a sample 
## including the confidence score respect to the other sample cells and also the minimal 
## changepoint respect to other cells
##
## Then, I join confidence scores of all segments according to their segment size
## to create a total of 20 null distributions (one per segment size component)
## 
################################################################################

### LIBRARIES
rm(list=ls(all=TRUE))

library(this.path)
library(data.table)
library(QDNAseq)
library(GenomicRanges)
library(foreach)
library(doMC)
library(YAPSA)
library(ComplexHeatmap)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(lemon)
library(RColorBrewer)
library(igraph)
library(qgraph)
library(lsa)
library(dplyr)
# Poisson
library(flexmix)
# Gaussian
library(mclust)
library(cowplot)


### PATHS
BASE=dirname(dirname(this.path()))
SCELL_DIR=dirname(this.path())
INPUT_DIR=file.path(SCELL_DIR,"intermediate")
OUTPUT_DIR=file.path(INPUT_DIR,"confidence")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

### PLOTTING THEME
theme_set(theme_tufte(base_size = 6, base_family = "ArialMT"))
theme_update(text = element_text(size = 6),
             axis.text = element_text(size = 6),
             axis.line = element_line(size = 0.5), 
             axis.ticks = element_line(size = 0.5),
             axis.ticks.length = unit(.1, "cm"))


### LOAD INPUT DATA
# We need as input the output from "sc_confidence_scores.R" script for each dataset
# and also the ExC matrix to know the segsize of each segment

# Minussi
confidence.minussi=dir(paste0(INPUT_DIR,"/confidence/minussi"),full=T,pattern="Confidence")
ExC.minussi=readRDS(paste0(INPUT_DIR,"/1.1.1_Minussi_ExC.rds"))

# BRCA1/2 KOs 
confidence.brca_ko=dir(paste0(INPUT_DIR,"/confidence/brca_ko"),full=T,pattern="Confidence")
ExC.brca_ko=readRDS(paste0(INPUT_DIR,"/3.1.1_BRCA_KO_ExC.rds"))

# PDAC PDOs
confidence.pdac_pdos=dir(paste0(INPUT_DIR,"/confidence/pdac_pdos"),full=T,pattern="Confidence")
ExC.pdac_pdos=readRDS(paste0(INPUT_DIR,"/4.1.1_PDAC-PDOs_ExC.rds"))

# RPE1 treated
confidence.rpe.drug=dir(paste0(INPUT_DIR,"/confidence/drug_rpe1"),full=T,pattern="Confidence")
ExC.rpe.drug=readRDS(paste0(INPUT_DIR,"/5.1.1_Drug_ExC.rds"))

# # Minussi 
# confidence.minussi_sub=dir(paste0(INPUT_DIR,"/confidence/minussi_sub"),full=T,pattern="Confidence")
# ExC.minussi=readRDS(paste0(INPUT_DIR,"/6.1.1_subMinussi_ExC.rds"))


################################################################################

## 1. Get all confidence files and join
confidence.files=c(confidence.minussi,confidence.brca_ko,confidence.pdac_pdos,confidence.rpe.drug)
lConfidence=lapply(confidence.files,function(thisFile){
    segs=readRDS(thisFile)
    return(segs)
})
names(lConfidence)=gsub("_Segments_Confidence.rds","",basename(confidence.files))

Confidence=rbindlist(lConfidence)
Confidence=Confidence[Confidence$segVal!=2,] #81464 events


## 1. Get all ExC files and join => same order than conidence files
ExC=rbind(ExC.minussi,rbind(ExC.brca_ko,rbind(ExC.pdac_pdos,ExC.rpe.drug)))
ExC.segsize=ExC[,1:20]
max.comp=colnames(ExC.segsize)[apply(ExC.segsize,1,which.max)]

lDistributions=lapply(colnames(ExC.segsize),function(thisComp){
    idxs=which(max.comp==thisComp)
    confs=Confidence$confidence[idxs]
    confs=confs[confs<=0.5]
    return(confs)
})
names(lDistributions)=colnames(ExC.segsize)
saveRDS(lDistributions,paste0(OUTPUT_DIR, "/NullDistribution_Confidence_Segsize.rds"))


