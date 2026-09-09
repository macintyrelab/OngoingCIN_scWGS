################################################################################
## This script is to process single-cell sequencing data from BRCA1/2 KOs
## Started: 01/08/24
## Updated: 02/06/26
##
## I derived copy number profiles from scWGS data generated in the lab using 
## PicoPlex/DLP protocol. These KO models may have copy number alterations induced 
## by deficiencies in the HR repair mechanism. These models have been done in a 
## TP53 null background which allows the accumulation of CNAs.
## 
## I selected those cells that have >15nrpbpcc. Once I have the profiles, I explore 
## the noise of them, exclude those cells that have high noise and also those that
## have overfragmented profiles, because they may be cycling.
##
## We also sequenced BRCA1/2 cells prior to treatment to PARPi and platins
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
PLOTS_DIR=file.path(SCELL_DIR,"plots/brca_ko")
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
sc.files.picoplex=dir(paste0(INPUT_DIR,"/KOs/picoplex"),full=T,pattern="absoluteCN")
sc.files.dlp=dir(paste0(INPUT_DIR,"/KOs/dlp"),full=T,pattern="absoluteCN")
# Profiles drug screening
sc.files.dlp367=dir(paste0(INPUT_DIR,"/KOs/dlp/DLP367"),full=T,pattern="absoluteCN") #TP53
sc.files.dlp369=dir(paste0(INPUT_DIR,"/KOs/dlp/DLP369"),full=T,pattern="absoluteCN") #BRCA1
sc.files.dlp371=dir(paste0(INPUT_DIR,"/KOs/dlp/DLP371"),full=T,pattern="absoluteCN") #BRCA2


# Threshold for excluding dividing cells
thresh.div=readRDS(paste0(INTER_DIR, "/3.0_BRCA_KO_Sphase_THRESH.rds"))

# Feature components
MIXMODELS=readRDS(paste0(BASE,"/data/SC_MIXMODELS.rds"))


################################################################################
#### 0. PROCESS SINGLE CELL PROFILES FROM CELL LINE MODELS ###############

### Get usable profiles 
# Remove samples with high noise & smooth diploid segments

# Picoplex data
SCELLS.CN.PPLEX=as.data.frame(do.call(rbind,lapply(sc.files.picoplex,function(thisFile){
    print(basename(thisFile))
    
    ## Step 1: Process data
    # read file
    cnobject=readRDS(thisFile)
    # tresh for dividing cells
    thresh=thresh.div[thresh.div$tech=="picoplex",]
    samp=thresh$sample[grepl(gsub("_absoluteCN_100kb.rds","",basename(thisFile)),thresh$sample)]
    if(length(samp)>1){samp=samp[1]}
    thresh=thresh$mean[thresh$sample==samp]+thresh$sd[thresh$sample==samp]
    
    # filter out noisy profiles
    if(length(cnobject)>0){
        segTabs=c()
        for(c in 1:length(cnobject)){
            cn=cnobject[[c]]
            copynumber=assayDataElement(cn, "copynumber")
            noise=apply(copynumber, 2, sdFUN, na.rm=TRUE)
            if(noise<=NOISE.PPLEX){
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
# IDs of single-cells
SCELLS.CN.PPLEX$sample[1:nrow(SCELLS.CN.PPLEX)]=sapply(1:nrow(SCELLS.CN.PPLEX),function(x)strsplit(SCELLS.CN.PPLEX$sample, '_')[[x]][1])
SCELLS.CN.PPLEX=SCELLS.CN.PPLEX[SCELLS.CN.PPLEX$sample!="brca2-H1",] #ploidy?

#add model info in cell lines
info=rbind(fread(paste0(INPUT_DIR,"/KOs/sc_brca2_picoplex_15082024.txt")),
           fread(paste0(INPUT_DIR,"/KOs/sc_brca1_screen_picoplex_15082024.txt")),
           fread(paste0(INPUT_DIR,"/KOs/sc_brca1_picoplex_12062024.txt")))
cn=SCELLS.CN.PPLEX[SCELLS.CN.PPLEX$sample%in%info$sample,1:5]
SCELLS.CN.PPLEX=left_join(cn,info,by="sample")
SCELLS.CN.PPLEX$tech="picoplex"

# DLP data
sc.files.dlp = c(sc.files.dlp, sc.files.dlp367, sc.files.dlp369, sc.files.dlp371)
SCELLS.CN.DLP=as.data.frame(do.call(rbind,lapply(sc.files.dlp,function(thisFile){
    print(basename(thisFile))
    
    ## Step 1: Process data
    # read file
    cnobject=readRDS(thisFile)
    # tresh for dividing cells
    thresh=thresh.div[thresh.div$tech=="dlp",]
    samp=thresh$sample[grepl(gsub("_absoluteCN_100kb.rds","",basename(thisFile)),thresh$sample)]
    thresh=thresh$mean[thresh$sample==samp]#+thresh$sd[thresh$sample==samp]
    
    # filter out noisy profiles
    if(length(cnobject)>0){
        segTabs=c()
        for(c in 1:length(cnobject)){
            cn=cnobject[[c]]
            copynumber=assayDataElement(cn, "copynumber")
            noise=apply(copynumber, 2, sdFUN, na.rm=TRUE)
            if(noise<=NOISE.DLP){
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
SCELLS.CN.DLP$tech="dlp"
SCELLS.CN.DLP$model=gsub("_absoluteCN.rds","",SCELLS.CN.DLP$model)

# Exclude early cycling samples => manual check
rm=c("DLP22_E5_S469_L005.markdup","DLP22_B6_S474_L005.markdup", "DLP371_D06.markdup", "DLP371_L04.markdup")
SCELLS.CN.DLP=SCELLS.CN.DLP[!SCELLS.CN.DLP$sample%in%rm,]

# Add info
info=rbind(fread(paste0(INPUT_DIR,"/KOs/sc_info_dlp75-76_20250519.txt")),
           fread(paste0(INPUT_DIR,"/KOs/sc_info_dlp367-372_drug_screenings.txt")))
colnames(info)[2]="cell"
info$cell=paste0(info$cell,".markdup")
SCELLS.CN.DLP$model[SCELLS.CN.DLP$model=="scWT&TP53_DLP75-76"]=info$model[match(SCELLS.CN.DLP$sample[SCELLS.CN.DLP$model=="scWT&TP53_DLP75-76"],info$cell)]

# Join datasets
SCELLS.CN=rbind(SCELLS.CN.PPLEX,SCELLS.CN.DLP)

# how many cells per sample?
stats_sample = as.data.frame(do.call(rbind,lapply(unique(SCELLS.CN$model), function(thisMod){
    cn=SCELLS.CN[SCELLS.CN$model==thisMod,]
    n=length(unique(cn$sample))
    return(c(thisMod,n))
})))
colnames(stats_sample)=c("sample","ncells")


# Save
SCELLS.CN$sample=factor(SCELLS.CN$sample, levels=unique(SCELLS.CN$sample)) #this is needed to have the correct order in the ExC matrix
saveRDS(SCELLS.CN, paste0(INTER_DIR, "/3.1.0_BRCA_KO_Segments.rds")) #141 cells => 138 if limiting to >5 cells per clone + 204 cells of drug screening


## Extract copy number features
print("Extract copy number features...")
cnprofiles=SCELLS.CN[,1:5]
lProfiles=split(cnprofiles, cnprofiles$sample)
ECNF=extractThreeCopynumberFeatures.sc(lProfiles, cores=CORES, prePath=file.path(BASE,"data"), rmNorm=RMNORM, size=SIZE, res=BIN.SIZE)
saveRDS(ECNF, paste0(INTER_DIR, "/3.1.0_BRCA_KO_ECNF.rds"))


## Get ExC matrix
print("Generate event-per-components matrix...")

FEATS=names(MIXMODELS)
lMats=lapply(FEATS, function(thisFeature) {
    print(thisFeature)
    thisEcnf=ECNF[[thisFeature]]
    thisModel=MIXMODELS[[thisFeature]]
    
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
        if(length(id)!=0){postDatScaled[id,]=matrix(rep(c(rep(0,8),1), length(id)), ncol = 9, nrow =length(id), byrow = T)}
    }
    # Create component matrix
    colnames(postDatScaled) = paste0(thisFeature, 1:ncol(postDatScaled))
    postDatScaled=as.matrix(postDatScaled)
    
    return(postDatScaled)
})
names(lMats)=FEATS
ExC=do.call(cbind, lMats)
saveRDS(ExC,paste0(INTER_DIR, "/3.1.1_BRCA_KO_ExC.rds"))


#### 1. GET BIN-LEVEL DATA FOR CREATE NULL DISTRIBUTION ########################
SAMPS=unique(SCELLS.CN$model)

### Get copy number values per bin
lSCELLS.CN=lapply(SAMPS,function(thisSamp){
    cn=SCELLS.CN[SCELLS.CN$model==thisSamp,1:5]
    return(cn)
})
names(lSCELLS.CN)=SAMPS

### Get copy number values per bin
lSCELLS.BIN.PPLEX=as.data.frame(do.call(rbind,lapply(sc.files.picoplex,function(thisFile){
    print(basename(thisFile))
    
    # read file
    cnobject=readRDS(thisFile)
    cells=sapply(1:length(cnobject),function(x) strsplit(names(cnobject), '_')[[x]][1])
    idx=which(cells%in%SCELLS.CN$sample)
    cnobject=cnobject[idx]
    
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
        out=cbind(chromosome,start,end,segVal,binVal,sample=strsplit(thisID, '_')[[1]][1])
        return(out)
    })))
    colnames(bins)=c("chromosome","start","end","segVal","binVal","sample")
    bins=bins[!is.na(bins$segVal),]
    bins[,c(2:5)]=sapply(2:5,function(x)as.numeric(bins[,x]))
    row.names(bins)=NULL
    return(bins)
})))

lSCELLS.BIN.DLP=as.data.frame(do.call(rbind,lapply(sc.files.dlp,function(thisFile){
    print(basename(thisFile))
    
    # read file
    cnobject=readRDS(thisFile)
    cells=names(cnobject)
    idx=which(cells%in%SCELLS.CN$sample)
    cnobject=cnobject[idx]
    
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
    return(bins)
})))

# split by sample/model
lSCELLS.BIN=lapply(SAMPS,function(thisSamp){
    samps=SCELLS.CN$sample[SCELLS.CN$model==thisSamp]
    tech=unique(SCELLS.CN$tech[SCELLS.CN$model==thisSamp])
    if(tech=="picoplex"){bins=lSCELLS.BIN.PPLEX[lSCELLS.BIN.PPLEX$sample%in%samps,]}
    if(tech=="dlp"){bins=lSCELLS.BIN.DLP[lSCELLS.BIN.DLP$sample%in%samps,]}
    return(bins)
})
names(lSCELLS.BIN)=SAMPS

# Save both lists
saveRDS(lSCELLS.BIN,paste0(INTER_DIR, "/3.1.1_BRCA_KO_Bins.rds")) #input for next step
saveRDS(lSCELLS.CN,paste0(INTER_DIR, "/3.1.1_BRCA_KO_Segments.rds")) #input for next step


## Next step! 
# Run in the cluster the "03_1b_sc_confidence_scores.R" script
# The script outputs a file per sample with segTables of all cells with confidence
# score respect to other cells and the minimal changepoint to other cells
# these files are saved in /confidence folder


############ END SCRIPT ! ######################################################