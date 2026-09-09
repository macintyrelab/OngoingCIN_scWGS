################################################################################
## This script is to derive feature encoding for single-event signatures for
## exploring ongoingCIN pan-cancer
## Started: 21/09/24
## Updated: 24/01/25
##
## I derived a feature encoding able to map signatures to individual copy number events.
## This feature encoding was used in the forecasting manuscript and I demonstrated its
## robustness for identifying and quantifying the activity of signatures in a sample.
##
## However, for moving this technology to the single-cell scenario, I need to define 
## signatures in the single-cell space. To do so, I need to adapt the feature component 
## distributions to the resolution of the scWGS data and to integer copy number values. 
## 
## Here, I derived feature encoding and signatures using copy number profiles from 
## TCGA, PCAWG and HMF cohorts. I extracted the 3 features for each copy number event 
## of the profiles. I then use these feature distributions to derive components 
## by mixture modeling, which will be used to derive single-event clonal signatures.
## For the segsize feature, I match segment sizes with single-cell resolution by 
## defining the minimal value as the bin resolution used to derive scWGS profiles (100kb)
## For the changepoint feature, I added an extra function to not allow changepoints < 1 
## (this is not possible at single-cell level) and rounded values to the closest integer 
## (mimicking the single-cell scenario where there is not subclonal changes)
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
SC_DIR=dirname(this.path())
INPUT_DIR=file.path(BASE,"data")
INTER_DIR=file.path(SC_DIR,"intermediate")
dir.create(file.path(INTER_DIR,"mixture"), showWarnings = FALSE, recursive = TRUE)
OUTPUT_DIR=file.path(SC_DIR,"outputs")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
PLOTS_DIR=file.path(SC_DIR,"plots")
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(SC_DIR,"data"), showWarnings = FALSE, recursive = TRUE)


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
WINDOW=100e3 #resolution of single-cell data x2 for calling a segment
RMNORM=TRUE
CORES=1
WIGGLE=0.1
IGNOREDELS=FALSE
SIZE=5000000


### LOAD INPUT DATA
# All data needed for running this script can be found downloaded from https://github.com/markowetzlab/Drews2022_CIN_Compendium
# Please, download all files and include them in a /data folder within the BASE path 

# TCGA profiles => we use all segments to derive feature space to be also applied at early evolutionary stages
TCGA.TABLE=rbind(readRDS(paste0(INPUT_DIR, "/tcga/0_TCGA_Segments_dCIN.rds")),
                 readRDS(paste0(INPUT_DIR, "/tcga/0_TCGA_Segments_noCIN.rds")))
TCGA.META=readRDS(paste0(INPUT_DIR, "/tcga/Metadata_TCGA_ASCAT_penalty70.rds")) #33 ttypes & 10793 samples


################################################################################
#### 1. EXTRACT COPY NUMBER FEATURES ###########################################

### 1.1. TCGA
# I extract features for all segments => segsize and changepoints are adapted to mimic a single-cell scenario
lProfiles=split(TCGA.TABLE, as.factor(TCGA.TABLE$sample))
ECNF=extractThreeCopynumberFeatures.sc(lProfiles, cores=CORES, prePath=INPUT_DIR, rmNorm=RMNORM, size=SIZE, res=WINDOW)
saveRDS(ECNF, paste0(INTER_DIR, "/1.1.1_TCGA_Segments_100kb_ECNF.rds")) #837428 CNAs & 10584 samples


#### 2. EXPLORE DISTRIBUTION OF COPY NUMBER FEATURES ###########################
# Here, I explore the distributions of the three features across the datasets
# These distributions will be then used for derive feature components

### 2.1. Explore distribution of features in a single-cell line scenario
# Run the "01a_Explore_scECNF_Segments.Rmd"


#### 3. FIT AND CHOOSE OPTIMAL MIXTURE COMPONENTS FOR FEATURES #################
ECNF.TCGA=readRDS(paste0(INTER_DIR, "/1.1.1_TCGA_Segments_100kb_ECNF.rds"))

### 3.1. SEGMENT SIZE
# It may follow a Gaussian distribution
# Apart from not having segments shorter than 100kb, the distribution of this feature 
# is similar to the one originally observed => 22 components in Drews et al. 2022
# Therefore, I just exclude the 2 components capturing shorter segments


### 3.2. CHANGEPOINT
# In this case, it may follow a Poisson distribution because it's an integer change
# The distribution of this feature is similar to the one originally observed
# However, only components representing integer copy number changes may be included.
# 10 components in Drews et al. 2022 (3 subclonal changes) => try a max components of 10

FEATURE="changepoint"
MINCOMP=1
MAXCOMP=10
FITDIST="pois"
PRIOR=0.001
NITER=1000
DECISION="BIC"
SEED=77777

## 3.2.1. TCGA
OUTNAME=paste0(INTER_DIR,"/mixture/fit_changepoint_tcga")

# Fit model => running with all data we get 2 components => we need representation of changes from 1 to 5
dat=as.numeric(ECNF.TCGA[[FEATURE]][,2])

# Fixed categories for dat<=5
Lambda=c(1,2,3,4,5)

# For data>5 (2%) => fit model for the tail of the distribution
dt=dat[dat>5] #remove zeros

fittedDist = fitComponent2(dt, dist=FITDIST, min_prior=PRIOR, min_comp=MINCOMP, 
                           max_comp=(MAXCOMP-length(Lambda)), niter=NITER, model_selection=DECISION, seed=SEED)

# Convert cluster components to a table for easier validation
if(FITDIST == "pois") {
    # Poisson distribution -> flexmix object
    listComp = unlist(fittedDist@components)
    allComps = sapply(listComp, function(thisComp) {unlist(thisComp@parameters)})
    dtComp = data.table(Lambda = sort(allComps))
} else {
    # Normal distribution -> Mclust object
    dtComp = data.table(Mean = fittedDist$parameters$mean,
                        SD = sqrt(fittedDist$parameters$variance$sigmasq))
    dtComp = dtComp[order(Mean),]
}
dtComp=rbind(as.data.frame(Lambda),dtComp)
saveRDS(dtComp, paste0(OUTNAME, ".rds"))


### 3.3. BREAKPOINTS PER SIDE
# It may follow a Poisson distribution
# The distribution of this feature is similar to the one originally observed in the
# forecasting work => 9 components in Hernando et al. 2025
# I just use the same mixture components


### 3.4. COMBINATION OF FEATURE COMPONENTS 
# Finally, I combine the solutions of the mixture modeling for each feature in one R object
MIXMODELS=readRDS(file.path(INPUT_DIR,"NEW_MIXMODELS.rds"))
scMIXMODELS=list()
scMIXMODELS[["segsize"]]=MIXMODELS[["segsize"]][3:22,] #remove segments shorter than 100kb
scMIXMODELS[["changepoint"]]=readRDS(paste0(INTER_DIR, "/mixture/fit_changepoint_tcga.rds"))
scMIXMODELS[["bpwindow"]]=MIXMODELS[["bpwindow"]]
saveRDS(scMIXMODELS, file=paste0(INPUT_DIR,"/SC_MIXMODELS.rds"))

############ END SCRIPT ! ######################################################




