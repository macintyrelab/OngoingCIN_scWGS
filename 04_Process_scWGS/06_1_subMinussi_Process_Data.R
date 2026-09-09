################################################################################
## This script is to process subclon-pseudobulk sequencing data from Minussi et al. 2022
## Started: 16/10/25
## Updated: 16/10/25
##
## I derived copy number profiles of the subclones found in each sample by
## concatenating fasta files of all cells included in each subclone. I fitted profiles using
## a range of ploidies close to the one reported in the paper. I used clonal structure
## reported in the paper for generating the subclon profiles: https://navinlabcode.github.io/ACT_paper/instructions.html
##
## In this case, all profiles have low noise levels and are not cycling. Therefore, 
## I don't need to clean up the input data
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
PLOTS_DIR=file.path(SCELL_DIR,"plots/minussi_subc")
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
WINDOW=100e3 #resolution of single-cell data x2 for calling a segment
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
# Subclonal structure and cells included in each subclone can be found here https://github.com/navinlabcode/ACT_paper
# Profiles powered
sub.files=dir(paste0(INPUT_DIR,"/minussi_subc"),full=T,pattern="absoluteCN")

# Feature components
MIXMODELS=readRDS(paste0(BASE,"/data/SC_MIXMODELS.rds"))


################################################################################
#### 0. GET SINGLE-CELL PROFILES FROM TNBC TUMOURS & CELL LINES ################

### Get usable profiles
# Remove samples with high noise & smooth diploid segments
SCLONES.CN=as.data.frame(do.call(rbind,lapply(sub.files,function(thisFile){
    print(basename(thisFile))
    
    ## Step 1: Process data
    # read file
    cnobject=readRDS(thisFile)
    
    # filter out noisy profiles
    if(length(cnobject)>0){
        segTabs=c()
        for(c in 1:length(cnobject)){
            cn=cnobject[[c]]
            copynumber=assayDataElement(cn, "copynumber")
            # noise=apply(copynumber, 2, sdFUN, na.rm=TRUE)
            thisClon=colnames(copynumber)
            segTab=CNpare::getSegTable(cn[,thisClon])
            segTab$sample=thisClon
            segTab[,2:4]=sapply(2:4,function(x)as.numeric(segTab[,x]))
            segTabs=rbind(segTabs,segTab)
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
            segTabs$model=gsub("_absoluteCN.rds","",basename(thisFile))
            return(segTabs) 
        }
    }
})))
rm=unique(SCLONES.CN$sample[grep("_s",SCLONES.CN$sample)]) #remove some superclone profiles
SCLONES.CN=SCLONES.CN[!SCLONES.CN$sample%in%rm,]
saveRDS(SCLONES.CN, paste0(INTER_DIR, "/6.0_subBRCA_Minussi_Powered_Segments.rds"))

# how many clones per sample?
stats_sample = as.data.frame(do.call(rbind,lapply(unique(SCLONES.CN$model), function(thisMod){
    cn=SCLONES.CN[SCLONES.CN$model==thisMod,]
    n=length(unique(cn$sample))
    return(c(thisMod,n))
})))
colnames(stats_sample)=c("sample","nclones") 


## Extract copy number features
print("Extract copy number features...")
cnprofiles=SCLONES.CN[,1:5]
lProfiles=split(cnprofiles, as.factor(cnprofiles$sample))
ECNF=extractThreeCopynumberFeatures.sc(lProfiles, cores=CORES, prePath=file.path(BASE,"data"), rmNorm=RMNORM, size=SIZE, res=BIN.SIZE)
saveRDS(ECNF, paste0(INTER_DIR, "/6.1.0_subMinussi_ECNF.rds"))


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
saveRDS(ExC,paste0(INTER_DIR, "/6.1.1_subMinussi_ExC.rds"))


#### 1. GET BIN-LEVEL DATA FOR TESTING CONFIDENCE OF UNIQUE CALLING ############
# This data is needed for computing segVal differences across profiles
# of the same sample, and it's then used for testing if unique events are real
# We can do the unique calling by sample, or by branch

SAMPS=unique(SCLONES.CN$model)


### 1.1. Analysis at sample level

## Split copy number profiles per sample
lSCLONES.CN=lapply(SAMPS,function(thisSamp){
    cn=SCLONES.CN[SCLONES.CN$model==thisSamp,1:5]
    return(cn)
})
names(lSCLONES.CN)=SAMPS

## Get copy number values per bin
lSCLONES.BIN=lapply(SAMPS,function(thisSamp){
    print(thisSamp)
    cnobject=readRDS(sub.files[grepl(thisSamp,sub.files)])

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
})
names(lSCLONES.BIN)=SAMPS
saveRDS(lSCLONES.BIN,paste0(INTER_DIR, "/6.1.1_Minussi_Bins_subpersample.rds")) #input for next step
saveRDS(lSCLONES.CN,paste0(INTER_DIR, "/6.1.1_Minussi_Segments_subpersample.rds")) #input for next step


### 1.2. Analysis at branch level
# In this case, I only apply it for TN1 and TN2 which have two different branches
SAMPS=c("TN1","TN2")
trees=dir(paste0(INPUT_DIR,"/minussi_subc/clonal_structure"),full=T,pattern="clustering")

## Split copy number profiles per sample
lSCLONES.CN=do.call(c,lapply(SAMPS,function(thisSamp){
    cn=SCLONES.CN[SCLONES.CN$model==thisSamp,1:5]
    tree=readRDS(trees[grep(thisSamp,trees)])
    branches=unique(tree$branch)
    
    cn.branch=lapply(branches,function(thisBranch){
        clones=paste0(unique(tree$subclones[tree$branch==thisBranch]),"_")
        cn.clones=as.data.frame(do.call(rbind,lapply(clones, function(thisClon){
            dt=cn[cn$sample==unique(cn$sample[grep(thisClon,cn$sample)]),]
            return(dt)
        })))
    })
    names(cn.branch)=paste0(thisSamp,"-",branches)
    return(cn.branch)
}))

## Get copy number values per bin
lSCLONES.BIN=do.call(c,lapply(SAMPS,function(thisSamp){
    print(thisSamp)
    cnobject=readRDS(sub.files[grepl(thisSamp,sub.files)])
    
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
    
    # separate per branches
    branches = names(lSCLONES.CN[grep(thisSamp,names(lSCLONES.CN))])
    bins.branch=lapply(branches,function(thisBranch){
        clones=unique(lSCLONES.CN[[thisBranch]]$sample)
        bins.clones=as.data.frame(do.call(rbind,lapply(clones, function(thisClon){
            dt=bins[bins$sample==unique(bins$sample[grep(thisClon,bins$sample)]),]
            return(dt)
        })))
    })
    names(bins.branch)=branches
    return(bins.branch)
}))
saveRDS(lSCLONES.BIN,paste0(INTER_DIR, "/6.1.1_Minussi_Bins_subperbranch.rds")) #input for next step
saveRDS(lSCLONES.CN,paste0(INTER_DIR, "/6.1.1_Minussi_Segments_subperbranch.rds")) #input for next step

## Next step! 
# Run in the cluster the "06_1b_sc_confidence_scores.R" script
# The script outputs a file per sample with segTables of all cells with confidence
# score respect to other cells and the minimal changepoint to other cells
# these files are saved in /confidence folder


############ END SCRIPT ! ######################################################