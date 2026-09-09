################################################################################
## This script is to call unique events in Minussi samples
## Started: 14/11/25
## Updated: 14/11/25
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
SCELL_DIR=file.path(BASE,"04_Process_scWGS")
NEWSIGS_DIR=file.path(BASE,"02_Derive_New_scCINSignatures")
INPUT_DIR=file.path(SCELL_DIR,"intermediate")
MISC_DIR=dirname(this.path())
# Create directories for outputs
OUTPUT_DIR=file.path(MISC_DIR,"outputs")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)


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
source(file.path(BASE, "scripts/uniqueEvents_functions.R"))


### PARAMETERS
RESOLUTION="100kb"
BIN.SIZE=as.numeric(substr(RESOLUTION,1,nchar(RESOLUTION)-2))*1000
WINDOW=as.numeric(BIN.SIZE*2)
THRESH=1.00 #% for defining a shared event as clonal/subclonal
SIGN.THRESH=0.01 #for defining significant unique segments


### LOAD INPUT DATA
# Profiles 
SCELLS.CN=readRDS(paste0(INPUT_DIR, "/1.0_Minussi_Powered_Segments.rds")) 
# ECNF 
ECNF=readRDS(paste0(INPUT_DIR, "/1.1.0_Minussi_ECNF.rds"))
# ExC matrix
ExC=readRDS(paste0(INPUT_DIR, "/1.1.1_Minussi_ExC.rds"))

# Feature components
MIXMODELS=readRDS(paste0(BASE,"/data/SC_MIXMODELS.rds"))
# Signature definitions
DEFS=readRDS(file=paste0(NEWSIGS_DIR,"/outputs/2_New_scCINSignatures_Definitions_Corrected.rds"))

# Null distribution matrix => for unique calling
lDistributions=readRDS(paste0(INPUT_DIR, "/confidence/NullDistribution_Confidence_Segsize.rds"))


################################################################################
###### CLASSIFY EVENTS AS CLONAL/SHARED/UNIQUE #################################

# Prepare events table
CN.EVENTS=SCELLS.CN[SCELLS.CN$segVal!=2,]
CN.EVENTS$segVal[CN.EVENTS$segVal<0]=0
CN.EVENTS$length=CN.EVENTS$end-CN.EVENTS$start
CN.EVENTS=CN.EVENTS[CN.EVENTS$length!=0,]
CN.EVENTS=cbind(CN.EVENTS,ExC)

# Classify CNAs as clonal, subclonal & unique
SAMPS=unique(SCELLS.CN$model)
lcntables=lapply(SAMPS,function(thisSamp){
    cn=SCELLS.CN[SCELLS.CN$model==thisSamp,1:5]
    return(cn)
})
names(lcntables)=SAMPS

# Get confidence matrix per segment
confidence.files=dir(paste0(INPUT_DIR,"/confidence/minussi"),full=T,pattern="Confidence")
lConfidence=lapply(confidence.files,function(thisFile){
    segs=readRDS(thisFile)
    return(segs)
})
names(lConfidence)=gsub("_Segments_Confidence.rds","",basename(confidence.files))

# Unify boundaries
print("Unifying matrix to classify events...")
lEvents=lapply(1:length(lcntables),function(thisTab){
    print(paste0("Sample ID: ",names(lcntables)[[thisTab]]))
    cntable=lcntables[[thisTab]]
    # unify segments
    unify=unifySegments.multiple(profiles=cntable,window=WINDOW,bin.size=BIN.SIZE)
    # round copy number values
    unify[,4:ncol(unify)]=sapply(4:ncol(unify),function(x)round(unify[,x]))
    # get unique & shared events
    events=rbindlist(classEvents.sc(unify,threshold=THRESH),idcol = T)
    events=as.data.frame(cbind(events[,c(2:5,1,6)]))
    colnames(events)[5]="sample"
    # remove normal segments
    events=events[events$segVal!=2,]
    return(events)
})
names(lEvents)=SAMPS
# Save list of events per sample
saveRDS(lEvents,file.path(OUTPUT_DIR,"all_events_minussi.rds"))



# Get unique segments with high confidence
print("Getting segments unique to cells...")
lUnique=lapply(1:length(lcntables), function(thisTab){
    print(paste0("Sample ID: ",names(lcntables)[[thisTab]]))
    #genomic locations of events
    events=lEvents[[thisTab]]
    unique=events[events$class=="unique",]
    unique.gr <- GRanges(seqnames = unique$chromosome, IRanges(start = unique$start, end = unique$end))
    #genomic locations of segments
    cntable=lcntables[[thisTab]]
    segments.gr = GRanges(seqnames = cntable$chromosome, IRanges(start = cntable$start, end = cntable$end))
    #get overlaps
    overlaps <- findOverlaps(unique.gr, segments.gr)
    # Get segments that map with the unique event in all samples
    seg.unique <- do.call(rbind,lapply(1:nrow(unique), function(thisUnq){
        # get the index of the events of all clones in the same location of unique event 
        idxs=overlaps@to[overlaps@from == thisUnq]
        segs=cntable[idxs,]
        samp=unique[thisUnq,5]
        usegs=segs[segs$sample==samp,]
        return(usegs)
    }))
    #remove duplicates => unique events can map in the same segment
    seg.unique=seg.unique[!duplicated(seg.unique),]
    seg.unique=seg.unique[order(seg.unique$chromosome,seg.unique$start,seg.unique$end),]
    
    if(nrow(seg.unique)>0){
        # Check that unique segments are really unique
        #confidence of event calls
        confis=lConfidence[[names(lcntables)[[thisTab]]]]
        seg.unique=left_join(seg.unique,confis,by=c("chromosome","start","end","segVal","sample"))
        
        #know segsize category to compare to null distribution
        length=seg.unique$end-seg.unique$start
        thisModel=MIXMODELS[["segsize"]]
        exc=sapply(1:nrow(thisModel), function(x) dnorm(x = length, mean = thisModel[[x,"Mean"]], sd = thisModel[[x,"SD"]]))
        if(is.null(nrow(exc))){exc=t(data.frame(exc))}
        # Normalize densities to probabilities
        exc=data.frame(exc/rowSums(exc))
        # Create component matrix
        colnames(exc) = paste0("segsize", 1:ncol(exc))
        exc=as.matrix(exc)
        max=colnames(exc)[apply(exc,1,which.max)]
        seg.unique$segsize=max
        # Compare
        stats <- as.data.frame(do.call(rbind,lapply(1:nrow(seg.unique), function(thisSeg){
            # compare to null distribution
            comp=seg.unique$segsize[thisSeg]
            min=seg.unique$minchange[thisSeg]
            dist=lDistributions[[comp]]
            if(length(dist)>1){t.pval=(t.test(min,dist,var.equal = T, alternative = "greater"))$p.value}else{t.pval=1}
            return(out=c(t.pval,mean(dist,na.rm=T),length(dist)))
        })))
        colnames(stats)=c("t.pvalue","mean.null","n")
        seg.unique=cbind(seg.unique,stats)
        return(seg.unique)
    }
})
names(lUnique)=SAMPS
dfUnique=rbindlist(lUnique,idcol = T)
colnames(dfUnique)[1]="model"
# Significant unique segments
dfUnique$t.qvalue=p.adjust(as.numeric(dfUnique$t.pvalue), method = "BH")
dfUnique$significant=ifelse(dfUnique$t.qvalue<SIGN.THRESH,"yes","no")
dfUnique$significant[dfUnique$minchange<0.7]="no" #minimal difference with other cell may be >0.7 to be considered as real
# Retain only high confident unique
dfUnique.filt=dfUnique[dfUnique$significant=="yes",]

# Stats
nUnique.model=dfUnique.filt[dfUnique.filt$significant=="yes",]%>%group_by(model)%>%summarise(n=n())

# Save unique table as ground truth
saveRDS(dfUnique.filt,file.path(OUTPUT_DIR,"real_unique_minussi.rds"))


############ END SCRIPT ! ######################################################