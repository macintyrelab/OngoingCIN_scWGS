################################################################################
## This script documents the modifications in the feature encoding for being
## able to map signatures to individual copy number events
## Started: 08/03/25
## Updated: 08/03/25
##
## The original feature encoding of our signature technology does not allow us to 
## evaluate genomic scars at copy number event level. This is because it contains 
## 3 copy number features (oscCIN, bpMB and bpChr) that were encoded per genomic 
## region rather than per individual copy number event
##
## Here, we introduce a new feature (breakpoints per window side) that replaces 
## these 3 features. We evaluate if the new feature capture the signal from the
## original features and it's able to capture changes in copy number distributions
## across the genome across samples
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


### LOAD INPUT DATA
## TCGA data downloaded from https://github.com/markowetzlab/Drews2022_CIN_Compendium
# TCGA copy number profiles
TCGA.TABLE=readRDS(paste0(INPUT_DIR, "/tcga/0_TCGA_Segments_dCIN.rds")) #6335 samples with CIN (>20 CNAs)
TCGA.TABLE$sample=as.character(TCGA.TABLE$sample)
lTCGA.TABLE=split(TCGA.TABLE, as.character(TCGA.TABLE$sample))
# TCGA signature activities
TCGA.EXPS.GOLD=readRDS(paste0(INPUT_DIR,"/tcga/Signature_Compendium_v5_Cosine-0.74_Activities_THRESH95_NAMESAPRIL21.rds"))
# Signature definitions
DEFS.GOLD=readRDS(paste0(INPUT_DIR, "/tcga/Signature_Compendium_v5_Cosine-0.74_Signatures_NAMESAPRIL21.rds"))


################################################################################
############# EXPLORE NEW FEATURE IN TCGA ######################################

### 1) Extract values of the BPWindow feature for all copy number events in the TCGA

BPWINDOW=getBPwindow(abs_profiles = lTCGA.TABLE, rmNorm = RMNORM, size = SIZE)


### 2) Explore distributions --> Supplementary Note Figure 2
# We selected TCGA samples with known CIN types (determined via CIN signatures)
# We selected samples with the highest activity for signatures where the replaced features had strong contributions

## 2.1. CX1-highest TCGA sample
samp=names(which.max(TCGA.EXPS.GOLD[,1]))
dat=BPWINDOW[BPWINDOW$ID==samp,]
pdf(paste0(PLOTS_DIR,"/01_Histogram_BPW_TCGA-CX1.pdf"), width = 45/25.4, height = 45/25.4)
p = ggplot(dat, aes(dat$value)) +
    geom_histogram(aes(y = after_stat(density)), boundary = 0, bins=100, fill = "grey70") + 
    # geom_density(color = "grey30", adjust = 2) +
    # ylim(0,3) +
    scale_x_continuous(limits=c(0,15), breaks=seq(0, 15, 2)) +
    labs(title = "CX1-enriched sample", subtitle=samp,
         x = "Breakpoints per 5Mb side", y = "") +
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 6),
          axis.text.y = element_text(size = 6),
          legend.position = "none")
print(p)
dev.off()


## 2.2. CX13-highest TCGA sample
samp=names(which.max(TCGA.EXPS.GOLD[,13]))
dat=BPWINDOW[BPWINDOW$ID==samp,]
pdf(paste0(PLOTS_DIR,"/01_Histogram_BPW_TCGA-CX13.pdf"), width = 45/25.4, height = 45/25.4)
p = ggplot(dat, aes(dat$value)) +
    geom_histogram(aes(y = after_stat(density)), boundary = 0, bins=100, fill = "grey70") + 
    # geom_density(color = "grey30", adjust = 2) +
    # ylim(0,3) +
    scale_x_continuous(limits=c(0,15), breaks=seq(0, 15, 2)) +
    labs(title = "CX13-enriched sample", subtitle=samp,
         x = "Breakpoints per 5Mb side", y = "") +
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 6),
          axis.text.y = element_text(size = 6),
          legend.position = "none")
print(p)
dev.off()


## 2.3. CX2-highest TCGA sample
samp=names(which.max(TCGA.EXPS.GOLD[,2]))
dat=BPWINDOW[BPWINDOW$ID==samp,]
pdf(paste0(PLOTS_DIR,"/01_Histogram_BPW_TCGA-CX2.pdf"), width = 45/25.4, height = 45/25.4)
p = ggplot(dat, aes(dat$value)) +
    geom_histogram(aes(y = after_stat(density)), boundary = 0, bins=100, fill = "grey70") + 
    # geom_density(color = "grey30", adjust = 2) +
    # ylim(0,3) +
    scale_x_continuous(limits=c(0,15), breaks=seq(0, 15, 2)) +
    labs(title = "CX2-enriched sample", subtitle=samp,
         x = "Breakpoints per 5Mb side", y = "") +
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 6),
          axis.text.y = element_text(size = 6),
          legend.position = "none")
print(p)
dev.off()


## 2.4. CX3-highest TCGA sample
samp=names(which.max(TCGA.EXPS.GOLD[,3]))
dat=BPWINDOW[BPWINDOW$ID==samp,]
pdf(paste0(PLOTS_DIR,"/01_Histogram_BPW_TCGA-CX3.pdf"), width = 45/25.4, height = 45/25.4)
p = ggplot(dat, aes(dat$value)) +
    geom_histogram(aes(y = after_stat(density)), boundary = 0, bins=100, fill = "grey70") + 
    # geom_density(color = "grey30", adjust = 2) +
    # ylim(0,3) +
    scale_x_continuous(limits=c(0,15), breaks=seq(0, 15, 2)) +
    labs(title = "CX3-enriched sample", subtitle=samp,
         x = "Breakpoints per 5Mb side", y = "") +
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 6),
          axis.text.y = element_text(size = 6),
          legend.position = "none")
print(p)
dev.off()


############# END OF THE SCRIPT ! ##############################################
