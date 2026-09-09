################################################################################
## This script is to process single clone sequencing data from Watson 2024
## Started: 08/08/24
## Updated: 12/03/25
##
## I derived copy number profiles from sWGS data generated in Watson et al 2024
## In this paper, normal mammary cell lines are used. They used reversine to induce
## missegregation errors, and then explore evolution in vitro. 
##
## I selected those clones that have >15nrpbpcc. Once I have the profiles, 
## I explore the noise of them and exclude those clones that have high noise or are 
## overfragmented because they may be cycling.
## Untreated clone has a perfect normal profile (2 copies in all chromosomes)
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
PLOTS_DIR=file.path(SCELL_DIR,"plots/watson")
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
source(file.path(BASE, "scripts/plotSCHeatmap.R"))

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
NOISE=0.18 #close to the median value of noise => run "01_1a_Watson_Explore_Noise.R"
RMNORM=TRUE
CORES=1
WIGGLE=0.1
IGNOREDELS=FALSE
SIZE=5000000
UNINFPRIOR=TRUE
RESOLUTION="100kb"
BIN.SIZE=as.numeric(substr(RESOLUTION,1,nchar(RESOLUTION)-2))*1000

### LOAD INPUT DATA
# Raw data can be downloaded from https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA634423
# Raw data is processed as indicated in the Methods section
# Profiles powered
sc.files=dir(paste0(INPUT_DIR,"/watson"),full=T,pattern="absoluteCN")
# Threshold for excluding dividing cells => run first the "02_1a_Watson_Explore_Noise.R"
thresh.div=readRDS(paste0(INTER_DIR, "/2.0_Watson_Sphase_THRESH.rds"))

# Feature components
MIXMODELS=readRDS(paste0(BASE,"/data/SC_MIXMODELS.rds"))


################################################################################
#### 0. PROCESS SINGLE CLONE PROFILES FROM hTERT–HMEC CELL LINE ################

## Get usable profiles
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

#add model info in cell lines
info=fread(paste0(INPUT_DIR,"/watson/Watson_SraRunTable_Metadata.txt"))
info=info[info$`Assay Type`=="WGS",]
colnames(info)[1]="sample"
info$sample=paste0(info$sample,".markdup")
SCELLS.CN=left_join(SCELLS.CN,info[,c(1,29,33,34,38)],by="sample")
SCELLS.CN=SCELLS.CN[!SCELLS.CN$Figure%in%c("Fig_1c","Fig_1b"),] #exclude single-cell data
SCELLS.CN=SCELLS.CN[,c(1:5,7,9,10)] #only mammary epithelial cells 
colnames(SCELLS.CN)=c("chromosome","start","end","segVal","sample","clone","model","experiment")
SCELLS.CN$experiment=ifelse(SCELLS.CN$experiment=="Fig_S5a" | SCELLS.CN$experiment=="Fig_S5b" | SCELLS.CN$experiment=="Fig_S5c","evolution","reversine")
SCELLS.CN$model=gsub("-aneuploid","",SCELLS.CN$model)
SCELLS.CN$sample=factor(SCELLS.CN$sample, levels=unique(SCELLS.CN$sample)) #this is needed to have the correct order in the ExC matrix
SCELLS.CN=SCELLS.CN[SCELLS.CN$experiment=="reversine",]
saveRDS(SCELLS.CN, paste0(INTER_DIR, "/2.1.0_Watson_Segments.rds")) # 41 clones => all 4N


## Extract copy number features
print("Extract copy number features...")
cnprofiles=SCELLS.CN[,1:5]
lProfiles=split(cnprofiles, as.factor(cnprofiles$sample))
ECNF=extractThreeCopynumberFeatures.sc(lProfiles, cores=CORES, prePath=file.path(BASE,"data"), rmNorm=RMNORM, size=SIZE, res=BIN.SIZE)
saveRDS(ECNF, paste0(INTER_DIR, "/2.1.0_Watson_ECNF.rds"))


## Get ExC matrix
print("Generate event-per-components matrix...")
ecnf=readRDS(paste0(INTER_DIR, "/2.1.0_Watson_ECNF.rds"))
cells=unique(cnprofiles$sample)

FEATS=names(MIXMODELS)
lMats=lapply(FEATS, function(thisFeature) {
    print(thisFeature)
    thisEcnf=ecnf[[thisFeature]]
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
saveRDS(ExC,paste0(INTER_DIR, "/2.1.1_Watson_Reversine_ExC.rds"))


############ END SCRIPT ! ######################################################






