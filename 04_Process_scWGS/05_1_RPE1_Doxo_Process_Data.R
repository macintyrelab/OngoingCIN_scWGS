################################################################################
## This script is to process single-cell sequencing data treated with doxorubicin
## and hydroxyurea for 10 days
## Started: 09/06/25
## Updated: 26/06/25
## 
## We treated RPE1-TP53 cells (107 cycles) with doxorubicin, which is known to induce micronuclei.
## We also treated with hydroxyrea, which is known to induce replication stress, but not 
## amplifications are really seen (not correct concentrations?)
## Cells were cultured for some days to allow the acquisition of copy number alterations
## caused by these treatments. 
## Given that these cells were under drug treatment and looked not really healthy, 
## we have a very low efficiency of the DLP runs. Therefore, I just manually selected 
## good cells (no overfragmented and low noise)
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
INPUT_DIR=file.path(BASE,"data")
INTER_DIR=file.path(SCELL_DIR,"intermediate")
dir.create(INTER_DIR, showWarnings = FALSE, recursive = TRUE)
# Create directories for outputs
OUTPUT_DIR=file.path(SCELL_DIR,"outputs")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
PLOTS_DIR=file.path(SCELL_DIR,"plots/drug_rpe1")
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


# For evaluation of observed noise
sdDiffTrim <- function(x, ..., trim=0.001, scale=TRUE) {
    if (scale)
        x <- x / mean(x, na.rm=TRUE)
    matrixStats::sdDiff(x, ..., trim=trim)
}
sdFUN <- match.fun(sdDiffTrim)


### PARAMETERS
SAVEINTERMEDIATE=TRUE
REMOVEQUIET=TRUE # Just to check that all samples have at least 20 CNAs. If not remove.
WINDOW=100e3 #resolution of single-cell data x2 for calling a segment
NOISE.PPLEX=0.26 #obtained in the "03_1a_BRCA_KO_Explore_Noise.R"
NOISE.DLP=0.21 #obtained in the "03_1a_BRCA_KO_Explore_Noise.R"
RMNORM=TRUE
CORES=1
WIGGLE=0.1
IGNOREDELS=FALSE
SIZE=5000000
UNINFPRIOR=TRUE
RESOLUTION="100kb"
BIN.SIZE=as.numeric(substr(RESOLUTION,1,nchar(RESOLUTION)-2))*1000


### LOAD INPUT DATA --> raw data will be available at publication
# Profiles
sc.files.dlp78=dir(paste0(INPUT_DIR,"/KOs/dlp/DLP78"),full=T,pattern="absoluteCN")
sc.files.dlp84=dir(paste0(INPUT_DIR,"/KOs/dlp/DLP84"),full=T,pattern="absoluteCN")

# Feature components
MIXMODELS=readRDS(paste0(BASE,"/data/SC_MIXMODELS.rds"))


################################################################################
#### 0. PROCESS SINGLE CELL PROFILES FROM CELL LINE MODELS DOX TREATED #########

## Select cells that pass QC filtering
# these are high-quality cells! manually checked for another project
cells=c(dir(paste0(INPUT_DIR,"/KOs/dlp/DLP78"),full=F,pattern="absolute.png"),
        dir(paste0(INPUT_DIR,"/KOs/dlp/DLP84"),full=F,pattern="absolute.png"))
cells=gsub("_absolute.png","",cells)
cells=cells[cells!="dlp84_A1.markdup"]

## Get usable profiles
# Read input file(s)
cnobject=c(readRDS(sc.files.dlp78),readRDS(sc.files.dlp84))
# Select ok cells
cnobject=cnobject[cells]

# Get segment tables
SCELLS.CN=as.data.frame(do.call(rbind,lapply(1:length(cnobject),function(c){
    cn=cnobject[[c]]
    copynumber=assayDataElement(cn, "copynumber")
    thisCell=colnames(copynumber)
    segTab=CNpare::getSegTable(cn[,thisCell])
    segTab$sample=thisCell
    segTab[,2:4]=sapply(2:4,function(x)as.numeric(segTab[,x]))
    return(segTab)
})))
SCELLS.CN$sample=gsub(".markdup","",SCELLS.CN$sample)
info=fread(paste0(INPUT_DIR,"/KOs/sc_info_dlp84_20250629.txt"))
SCELLS.CN$model=info$model[match(SCELLS.CN$sample,info$sample)]
SCELLS.CN$model[SCELLS.CN$sample=="dlp78-F8"]="Doxo"
SCELLS.CN$model[is.na(SCELLS.CN$model)]="DMSO"

# round copy number values
SCELLS.CN$segVal=round(SCELLS.CN$segVal)
saveRDS(SCELLS.CN,paste0(INTER_DIR, "/5.1.0_Doxo_Segments.rds")) #input for next step

# how many cells per sample?
stats_sample = as.data.frame(do.call(rbind,lapply(unique(SCELLS.CN$model), function(thisMod){
    cn=SCELLS.CN[SCELLS.CN$model==thisMod,]
    n=length(unique(cn$sample))
    return(c(thisMod,n))
})))
colnames(stats_sample)=c("sample","ncells")


## Extract copy number features
print("Extract copy number features...")
cnprofiles=SCELLS.CN[,1:5]
lProfiles=split(cnprofiles, as.factor(cnprofiles$sample))
ECNF=extractThreeCopynumberFeatures.sc(lProfiles, cores=CORES, prePath=file.path(BASE,"data"), rmNorm=RMNORM, size=SIZE, res=BIN.SIZE)
saveRDS(ECNF, paste0(INTER_DIR, "/5.1.0_Doxo_ECNF.rds"))


## Get ExC matrix
print("Generate event-per-components matrix...")
cells=unique(cnprofiles$sample)

FEATS=names(MIXMODELS)
lMats=lapply(FEATS, function(thisFeature) {
    print(thisFeature)
    thisEcnf=ECNF[[thisFeature]]
    thisModel=MIXMODELS[[thisFeature]]
    
    # Select cells
    thisEcnf=thisEcnf[thisEcnf$ID%in%cells,]
    
    # We want a posterior, hence likelihood (density) times prior (weight)
    dat=as.numeric(thisEcnf[,2])
    print("Calculating probabilities...")
    if(thisFeature!="segsize"){
        # Poisson model
        print("Poisson-based posterior")
        postDatUnscaled = sapply(1:nrow(thisModel), function(x) dpois(x = dat, lambda = thisModel[[x,"Lambda"]]))
    }
    else{
        # Gaussian model
        print("Gauss-based posterior")
        postDatUnscaled = sapply(1:nrow(thisModel), function(x) dnorm(x = dat, mean = thisModel[[x,"Mean"]], sd = thisModel[[x,"SD"]]))
    }
    # Normalize densities to probabilities
    postDatScaled = data.frame(postDatUnscaled/rowSums(postDatUnscaled))
    # Make sure that focal amplifications with huge number of copies have weight in the chp=8 
    # Make sure that windows with a lot of segments have weight in the bpw=8
    id = which(rowSums(is.na(postDatScaled)) == ncol(postDatScaled))
    if (thisFeature == "changepoint"){
        if(length(id)!=0){postDatScaled[id,]=matrix(rep(c(rep(0,7),1), length(id)), ncol = 8, nrow =length(id), byrow = T)}
    }
    if (thisFeature == "bpwindow"){
        if(length(id)!=0){postDatScaled[id,]=matrix(rep(c(rep(0,7),1), length(id)), ncol = 8, nrow =length(id), byrow = T)}
    }
    # Create component matrix
    colnames(postDatScaled) = paste0(thisFeature, 1:ncol(postDatScaled))
    postDatScaled=as.matrix(postDatScaled)
    
    return(postDatScaled)
})
names(lMats)=FEATS
ExC=do.call(cbind, lMats)
saveRDS(ExC,paste0(INTER_DIR, "/5.1.1_Doxo_ExC.rds"))


#### 1. GET BIN-LEVEL DATA FOR CREATE NULL DISTRIBUTION ########################

### Split copy number profiles per sample
SAMPS=unique(SCELLS.CN$model)
lSCELLS.CN=lapply(SAMPS,function(thisSamp){
    cn=SCELLS.CN[SCELLS.CN$model==thisSamp,1:5]
    return(cn)
})
names(lSCELLS.CN)=SAMPS

### Get copy number values per bin
lSCELLS.BIN=lapply(SAMPS,function(thisSamp){
    print(thisSamp)
    cells=paste0(unique(SCELLS.CN$sample[SCELLS.CN$model==thisSamp]),".markdup")
    cnobject=cnobject[cells]
    
    # get raw values per bin
    bins=as.data.frame(do.call(rbind,lapply(1:length(cnobject),function(c){
        cn=cnobject[[c]]
        binnumber=assayDataElement(cn,"copynumber")
        copynumber=assayDataElement(cn, "segmented")
        thisID=colnames(cn)
        #locations
        chromosome=gsub(":.+$", "", row.names(cn))
        pos=gsub("^.+:", "", row.names(cn))
        start=as.numeric(gsub("-.+$", "", pos))
        end=as.numeric(gsub("^.+-", "", pos))
        
        #copy number values
        segVal=copynumber[,thisID]
        binVal=binnumber[,thisID]
        
        #outputs
        out=cbind(chromosome,start,end,segVal,binVal,sample=thisID)
        return(out)
    })))
    colnames(bins)=c("chromosome","start","end","segVal","binVal","sample")
    bins=bins[!is.na(bins$segVal),]
    bins[,c(2:5)]=sapply(2:5,function(x)as.numeric(bins[,x]))
    bins$sample=gsub(".markdup","",bins$sample)
    row.names(bins)=NULL
    return(bins)
})
names(lSCELLS.BIN)=SAMPS
saveRDS(lSCELLS.BIN,paste0(INTER_DIR, "/5.1.1_Doxo_Bins.rds")) #input for next step
saveRDS(lSCELLS.CN,paste0(INTER_DIR, "/5.1.1_Doxo_Segments.rds")) #input for next step


## Next step! 
# Run in the cluster the "05_1b_sc_confidence_scores.R" script
# The script outputs a file per sample with segTables of all cells with confidence
# score respect to other cells and the minimal changepoint to other cells
# these files are saved in /confidence folder
