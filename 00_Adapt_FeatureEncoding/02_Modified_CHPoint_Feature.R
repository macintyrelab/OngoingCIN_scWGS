################################################################################
## This script documents the modifications in the feature encoding for being
## able to map signatures to individual copy number events
## Started: 08/03/25
## Updated: 08/03/25
##
## To improve signature mapping to individual events, we also refined the 
## original changepoint feature to avoid linking an event next to a focal 
## amplification to a replication stress signature if it has low copy number value
##
## Here, we compare changepoint values for events surrounding focal amplifications
## using the original and the modified encoding. We select surrounding events that 
## have < 4 copies and are located after a focal amplification.
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
FEATS_DIR=dirname(this.path())
INPUT_DIR=file.path(BASE,"data")
OUTPUT_DIR=file.path(FEATS_DIR,"outputs")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
PLOTS_DIR=file.path(FEATS_DIR,"plots")
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)

### PLOTTING THEME
theme_set(theme_tufte(base_size = 6, base_family = "ArialMT"))
theme_update(text = element_text(size = 6),
             axis.text = element_text(size = 6),
             axis.line = element_line(size = 0.5), 
             axis.ticks = element_line(size = 0.5),
             axis.ticks.length = unit(.1, "cm"))

### FUNCTIONS
source(file.path(BASE, "scripts/newCNSigs_functions.R"))
source(file.path(BASE, "scripts/main_functions.R"))


### PARAMETERS
REMOVEQUIET=TRUE # Just to check that all samples have at least 20 CNAs. If not remove.
RMNORM=TRUE
CORES=1
WIGGLE=0.1
IGNOREDELS=FALSE
SIZE=5000000
UNINFPRIOR=TRUE
AMP_THRESH=8 #to define focal amplifications
LOWCN_THRESH=4 #to define surrounding events to be explored

### LOAD INPUT DATA
## TCGA data downloaded from https://github.com/markowetzlab/Drews2022_CIN_Compendium
# TCGA copy number profiles
TCGA.TABLE=readRDS(paste0(INPUT_DIR, "/tcga/0_TCGA_Segments_dCIN.rds")) #6335 samples with CIN (>20 CNAs)
TCGA.TABLE$sample=factor(TCGA.TABLE$sample)
lTCGA.TABLE=split(TCGA.TABLE, TCGA.TABLE$sample)
# TCGA signature activities
TCGA.EXPS.GOLD=readRDS(paste0(INPUT_DIR,"/tcga/Signature_Compendium_v5_Cosine-0.74_Activities_THRESH95_NAMESAPRIL21.rds"))
# Signature definitions
DEFS.GOLD=readRDS(paste0(INPUT_DIR, "/tcga/Signature_Compendium_v5_Cosine-0.74_Signatures_NAMESAPRIL21.rds"))


################################################################################
############# EXPLORE CHANGEPOINT FEATURE IN TCGA ##############################

### 1) Extract values of the Changepoint feature for all copy number events in the TCGA

## 1.1. Original encoding
CHPOINT.org=getChangepointCN2(abs_profiles = lTCGA.TABLE, rmNorm = RMNORM)

## 1.2. Modified encoding
CHPOINT.new=getChangepointCN.new(abs_profiles = lTCGA.TABLE, rmNorm = RMNORM)


### 2) Extract idx of after focal amplifications
TCGA.TABLE$ID=paste0(TCGA.TABLE$sample, "|", TCGA.TABLE$chromosome, ":", TCGA.TABLE$start, "-", TCGA.TABLE$end)

# Get table with events that are not normal
TCGA.EVENTS=TCGA.TABLE[TCGA.TABLE$segVal!=2,]
# Avoiding artefacts
TCGA.EVENTS$segVal[TCGA.EVENTS$segVal<0]=0
TCGA.EVENTS$length=TCGA.EVENTS$end-TCGA.EVENTS$start
TCGA.EVENTS=TCGA.EVENTS[TCGA.EVENTS$length!=0,]

idx = which(TCGA.EVENTS$segVal>=AMP_THRESH)
NEXT = TCGA.EVENTS[idx+1,]
NEXT = NEXT[NEXT$segVal<LOWCN_THRESH,]

idx = which(TCGA.EVENTS$ID%in%NEXT$ID)

## 2.1. Original encoding
next_segs.org = CHPOINT.org[idx,]
amps_segs.org = CHPOINT.org[idx-1,]
prev_segs.org$type = "prev"
next_segs.org$type = "next"
amps_segs.org$type = "amp"
chp.org = rbind(next_segs.org,amps_segs.org)
chp.org$version="old"

## 2.2. Modified encoding
next_segs.new = CHPOINT.new[idx,]
amps_segs.new = CHPOINT.new[idx-1,]
prev_segs.org$type = "prev"
next_segs.new$type = "next"
amps_segs.new$type = "amp"
chp.new = rbind(next_segs.new,amps_segs.new)
chp.new$version="new"

### 3) Plotting distributions of changepoint values --> Supplementary Note Figure 2
df=rbind(chp.org,chp.new)
df$value=as.numeric(df$value)
df$type=factor(df$type,levels = c("amp","next"))
df$version=factor(df$version,levels = c("old","new"))

pdf(file=paste0(PLOTS_DIR, "/02_Distribution_CHP_original_vs_new.pdf"), width = 90/25.4, height = 45/25.4)
p = ggplot(df,aes(x=type,y=value))+
    ylim(0,100) +
    geom_boxplot(outlier.shape = NA, color="grey30", fill=NA)+
    geom_point(aes(fill=version), col="grey50", pch = 21, position = position_jitter(0.2), alpha = 0.2, size = 0.8, stroke = 0.1) +
    scale_fill_manual(values = c("skyblue1","tomato1")) +
    facet_wrap(~version, scales = "free") +
    labs(x = "", y = "Changepoint feature values") +
    theme(axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          legend.position = "none")
print(p)
dev.off()


############# END OF THE SCRIPT ! ##############################################
