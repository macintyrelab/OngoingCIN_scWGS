################################################################################
## This script is to dissect ongoing CIN in pancreatic cancer organoids
## Started: 06/04/25
## Updated: 22/04/25
##
## I derived copy number profiles from single-cell sequencing data generated in 
## the lab using the DLP/PicoPlex protocol. I selected those cells that have >30nrpbpcc, 
## and I excluded cells with noisy and/or overfragmented profiles
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
NEWSIGS_DIR=file.path(BASE,"02_Derive_New_scCINSignatures")
INPUT_DIR=file.path(BASE,"data")
INTER_DIR=file.path(SCELL_DIR,"intermediate")
dir.create(INTER_DIR, showWarnings = FALSE, recursive = TRUE)
# Create directories for outputs
OUTPUT_DIR=file.path(SCELL_DIR,"outputs")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
PLOTS_DIR=file.path(SCELL_DIR,"plots/pdac_pdos")
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
sc.files=dir(paste0(INPUT_DIR,"/pdos_pdac"),full=T,pattern="absoluteCN.rds") 

# Feature components
MIXMODELS=readRDS(paste0(BASE,"/data/SC_MIXMODELS.rds"))


################################################################################
#### 0. PROCESS SINGLE CELL PROFILES FROM PDAC ORGANOIDS #######################

## Select cells that are used in the phylogenetic trees
# these are high-quality cells! manually checked for another project
cells=c(dir(paste0(INPUT_DIR,"/pdos_pdac/DLP62-64"),full=F,pattern="absolute.png"),
        dir(paste0(INPUT_DIR,"/pdos_pdac/DLP73"),full=F,pattern="absolute.png"),
        dir(paste0(INPUT_DIR,"/pdos_pdac/scMUR72"),full=F,pattern="absolute.png"))
cells=gsub("_absolute.png","",cells)

## Get usable profiles
# Read input file(s)
cnobject=c(readRDS(sc.files[1]),readRDS(sc.files[2]),readRDS(sc.files[3]),readRDS(sc.files[4]))
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
  thisModel=ifelse(grepl("DLP62",thisCell),"MUR54",
                   ifelse(grepl("DLP63",thisCell),"MUR59",
                          ifelse(grepl("DLP64",thisCell),"MUR67",
                                 ifelse(grepl("DLP73",thisCell),"MUR79","MUR72"))))
  segTab$model=thisModel
  return(segTab)
})))

# Avoiding artefacts
SCELLS.CN$segVal[SCELLS.CN$segVal<0]=0
SCELLS.CN=SCELLS.CN[(SCELLS.CN$end-SCELLS.CN$start)!=0,]

# round copy number values
SCELLS.CN$segVal=round(SCELLS.CN$segVal)

# Merge segments only when segments with same copy number value follow each other 
df=idSmoothingTargets(SCELLS.CN, WIGGLE = 0, colNameSegVal = "segVal", colNameChr = "chromosome", IGNOREDELS = IGNOREDELS)
# Split by sample name
lRaw=split(df, df$sample)
# Smooth segments by taking the weighted average of the segVal and their lengths
lSmooth=smoothSegments(lRaw, CORES = 1, SMOOTHINGFACTOR = 0, colNameMerge = "segVal", colNameChr = "chromosome",
                       colNameStart = "start", colNameEnd = "end", IGNOREDELS = IGNOREDELS, asDf = FALSE)

SCELLS.CN=data.frame(rbindlist(lSmooth))
SCELLS.CN$sample=gsub(".markdup","",SCELLS.CN$sample)

# how many cells per sample?
stats_sample = as.data.frame(do.call(rbind,lapply(unique(SCELLS.CN$model), function(thisMod){
  cn=SCELLS.CN[SCELLS.CN$model==thisMod,]
  n=length(unique(cn$sample))
  return(c(thisMod,n))
})))
colnames(stats_sample)=c("sample","ncells")

# I select samples with at least >10 cells
SAMPS=unique(stats_sample$sample[stats_sample$ncells>=10])
SCELLS.CN=SCELLS.CN[SCELLS.CN$model%in%SAMPS,]
SCELLS.CN$sample=factor(SCELLS.CN$sample,levels=unique(SCELLS.CN$sample)) # this is needed to have the correct order in the ExC matrix
saveRDS(SCELLS.CN, paste0(INTER_DIR, "/4.1.0_PDAC-PDOs_Segments.rds")) # 122 cells from 5 PDOs


## Extract copy number features
print("Extract copy number features...")
cnprofiles=SCELLS.CN[,1:5]
lProfiles=split(cnprofiles, as.factor(cnprofiles$sample))
ECNF=extractThreeCopynumberFeatures.sc(lProfiles, cores=CORES, prePath=file.path(BASE,"data"), rmNorm=RMNORM, size=SIZE, res=BIN.SIZE)
saveRDS(ECNF, paste0(INTER_DIR, "/4.1.0_PDAC-PDOs_ECNF.rds"))


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
  # Make sure that windows with a lot of segments have weight in the bpw=9
  id = which(rowSums(is.na(postDatScaled)) == ncol(postDatScaled))
  if (thisFeature == "changepoint"){
    if(length(id)!=0){postDatScaled[id,]=matrix(rep(c(rep(0,7),1), length(id)), ncol = 8, nrow =length(id), byrow = T)}
  }
  if (thisFeature == "bpwindow"){
    if(length(id)!=0){postDatScaled[id,]=matrix(rep(c(rep(0,8),1), length(id)), ncol = 9, nrow =length(id), byrow = T)}
  }
  # Create component matrix
  colnames(postDatScaled) = paste0(thisFeature, 1:ncol(postDatScaled))
  postDatScaled=as.matrix(postDatScaled)
  
  return(postDatScaled)
})
names(lMats)=FEATS
ExC=do.call(cbind, lMats)
saveRDS(ExC,paste0(INTER_DIR, "/4.1.1_PDAC-PDOs_ExC.rds"))


#### 1. GET BIN-LEVEL DATA FOR CREATE NULL DISTRIBUTION ########################

### Split copy number profiles per sample
# Split the MUT79 PDO in two clones & 3 subclones => one is diploid & another is tetraploid (with 2 subclones)
tetr.A=c("DLP73-D3","DLP73-F10","DLP73-H5","DLP73-E7","DLP73-B9","DLP73-G8","DLP73-G5","DLP73-E4","DLP73-D1","DLP73-H7")
tetr.B=c("DLP73-C7","DLP73-B7","DLP73-B1","DLP73-E9","DLP73-B3","DLP73-G9","DLP73-B5")
tetr=c(tetr.A,tetr.B)
SCELLS.CN$model[SCELLS.CN$sample%in%tetr]="MUR79-tetr"

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
saveRDS(lSCELLS.BIN,paste0(INTER_DIR, "/4.1.1_PDAC-PDOs_Bins.rds")) #input for next step
saveRDS(lSCELLS.CN,paste0(INTER_DIR, "/4.1.1_PDAC-PDOs_Segments.rds")) #input for next step


## Next step! 
# Run in the cluster the "04_1b_sc_confidence_scores.R" script
# The script outputs a file per sample with segTables of all cells with confidence
# score respect to other cells and the minimal changepoint to other cells
# these files are saved in /confidence folder
