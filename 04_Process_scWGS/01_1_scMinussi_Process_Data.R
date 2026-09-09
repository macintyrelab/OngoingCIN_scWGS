################################################################################
## This script is to process single-cell sequencing data from Minussi et al. 2022
## Started: 07/06/24
## Updated: 09/02/25
##
## I derived copy number profiles from scWGS data downloaded from Minussi et al. 2022
## To do so, I selected those cells that have >15nrpbpcc. I fitted profiles using
## a range of ploidies close to the one reported in the paper. 
## Once I have the profiles, I explore the noise of them and exclude those cells
## that have high noise becase they may be cycling. 
##
## This data is used for generating the null distribution for calling unique events
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
PLOTS_DIR=file.path(SCELL_DIR,"plots/minussi")
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
REMOVEQUIET=TRUE # Just to check that all samples have at least 20 CNAs. If not remove.
WINDOW=100e3 #resolution of single-cell data x2 for calling a segment
NOISE=0.18 #obtained from the "01_1a_scMinussi_Explore_Noise.R" script
RMNORM=TRUE
CORES=1
WIGGLE=0.1
IGNOREDELS=FALSE
SIZE=5000000
UNINFPRIOR=TRUE
RESOLUTION="100kb"
BIN.SIZE=as.numeric(substr(RESOLUTION,1,nchar(RESOLUTION)-2))*1000


### LOAD INPUT DATA
# Raw data can be downloaded from https://www.ncbi.nlm.nih.gov/sra?term=PRJNA629885
# Raw data is processed as indicated in the Methods section
# Profiles powered
sc.files=dir(paste0(INPUT_DIR,"/minussi_sc"),full=T,pattern="absoluteCN")
# Threshold for excluding dividing cells
thresh.div=readRDS(paste0(INTER_DIR, "/1.0_scMinussi_Powered_Sphase_THRESH.rds"))

# Feature components
MIXMODELS=readRDS(paste0(BASE,"/data/SC_MIXMODELS.rds"))


################################################################################
#### 0. GET SINGLE-CELL PROFILES FROM TNBC TUMOURS & CELL LINES ################

### Get usable profiles
# Remove samples with high noise & smooth diploid segments
SCELLS.CN=as.data.frame(do.call(rbind,lapply(sc.files,function(thisFile){
    print(basename(thisFile))
    
    ## Step 1: Process data
    # read file
    cnobject=readRDS(thisFile)
    # tresh for dividing cells
    samp=thresh.div$sample[grepl(gsub("_absoluteCN_100kb.rds","",basename(thisFile)),thresh.div$sample)]
    thresh=thresh.div$thresh[thresh.div$sample==samp]
    
    # filter out noisy profiles
    if(length(cnobject)>0){
        segTabs=c()
        for(c in 1:length(cnobject)){
            cn=cnobject[[c]]
            copynumber=assayDataElement(cn, "copynumber")
            noise=apply(copynumber, 2, sdFUN, na.rm=TRUE)
            if(noise<=NOISE){
                thisClon=colnames(copynumber)
                segTab=CNpare::getSegTable(cn[,thisClon])
                segTab$sample=thisClon
                segTab[,2:4]=sapply(2:4,function(x)as.numeric(segTab[,x]))
                if(nrow(segTab)<thresh){
                    segTabs=rbind(segTabs,segTab)
                }
            }
        }
        
        # for those samples with powered & non-noisy cells
        if(!is.null(segTabs)){
            # avoiding artefacts
            segTabs$segVal[segTabs$segVal<0]=0
            segTabs=segTabs[(segTabs$end-segTabs$start)!=0,]
            
            # round copy number values
            segTabs$segVal=round(segTabs$segVal)
            
            # Merge segments only when segments with same copy number value follow each other 
            df=idSmoothingTargets(segTabs, WIGGLE = 0, colNameSegVal = "segVal", colNameChr = "chromosome", IGNOREDELS = IGNOREDELS)
            # Split by sample name
            lRaw=split(df, df$sample)
            # Smooth segments by taking the weighted average of the segVal and their lengths
            lSmooth=smoothSegments(lRaw, CORES = 1, SMOOTHINGFACTOR = 0, colNameMerge = "segVal", colNameChr = "chromosome",
                                   colNameStart = "start", colNameEnd = "end", IGNOREDELS = IGNOREDELS, asDf = FALSE)
            
            segTabs=data.frame(rbindlist(lSmooth))
            segTabs$model=gsub("_absoluteCN_100kb.rds","",basename(thisFile))
            return(segTabs) 
        }
    }
})))
SCELLS.CN$model=gsub("sc_","",SCELLS.CN$model) # 7 samples & 1125 cells

# how many cells per sample?
stats_sample = as.data.frame(do.call(rbind,lapply(unique(SCELLS.CN$model), function(thisMod){
    cn=SCELLS.CN[SCELLS.CN$model==thisMod,]
    n=length(unique(cn$sample))
    return(c(thisMod,n))
})))
colnames(stats_sample)=c("sample","ncells") 

# I am going to use this dataset to create the null distribution
# I select samples with at least >10 cells
SAMPS=unique(stats_sample$sample[as.numeric(stats_sample$ncells)>=10])
SCELLS.CN=SCELLS.CN[SCELLS.CN$model%in%SAMPS,]
saveRDS(SCELLS.CN, paste0(INTER_DIR, "/1.0_Minussi_Powered_Segments.rds")) 

## Extract copy number features
print("Extract copy number features...")
cnprofiles=SCELLS.CN[,1:5]
lProfiles=split(cnprofiles, as.factor(cnprofiles$sample))
ECNF=extractThreeCopynumberFeatures.sc(lProfiles, cores=CORES, prePath=file.path(BASE,"data"), rmNorm=RMNORM, size=SIZE, res=BIN.SIZE)
saveRDS(ECNF, paste0(INTER_DIR, "/1.1.0_Minussi_ECNF.rds"))


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
        if(length(id)!=0){postDatScaled[id,]=matrix(rep(c(rep(0,8),1), length(id)), ncol = 1, nrow =length(id), byrow = T)}
    }
    # Create component matrix
    colnames(postDatScaled) = paste0(thisFeature, 1:ncol(postDatScaled))
    postDatScaled=as.matrix(postDatScaled)
    
    return(postDatScaled)
})
names(lMats)=FEATS
ExC=do.call(cbind, lMats)
saveRDS(ExC,paste0(INTER_DIR, "/1.1.1_Minussi_ExC.rds"))


#### 1. GET BIN-LEVEL DATA FOR CREATE NULL DISTRIBUTION ########################

### Split copy number profiles per sample
lSCELLS.CN=lapply(SAMPS,function(thisSamp){
    cn=SCELLS.CN[SCELLS.CN$model==thisSamp,1:5]
    return(cn)
})
names(lSCELLS.CN)=SAMPS

### Get copy number values per bin
lSCELLS.BIN=lapply(SAMPS,function(thisSamp){
    print(thisSamp)
    cnobject=readRDS(sc.files[grepl(thisSamp,sc.files)])

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
    row.names(bins)=NULL
    #remove cells that do not pass filters
    samps=lSCELLS.CN[[thisSamp]]$sample
    bins=bins[bins$sample%in%samps,]
    #output
    return(bins)
})
names(lSCELLS.BIN)=SAMPS
saveRDS(lSCELLS.BIN,paste0(INTER_DIR, "/1.1.1_Minussi_Bins.rds")) #input for next step
saveRDS(lSCELLS.CN,paste0(INTER_DIR, "/1.1.1_Minussi_Segments.rds")) #input for next step

## Next step! 
# Run in the cluster the "01_1b_sc_confidence_scores.R" script
# The script outputs a file per sample with segTables of all cells with confidence
# score respect to other cells and the minimal changepoint to other cells
# these files are saved in /confidence folder


############ END SCRIPT ! ######################################################