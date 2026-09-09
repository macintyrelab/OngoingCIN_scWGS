################################################################################
## This script is for evaluating different parameters used in our single-cell framework
## Started: 16/04/2026
## Updated: 20/04/2026
##
## Here, I evaluate the different thresholds applied to select high-confidence cells
## for quantifying ongoing CIN. 
##  1) 100kb bin resolution (already solved) --> we have tested different bin resolutions 
##     for quantifying ongoingCIN and show that bin sizes >100kb precludes signature detection
##  2) Sequencing depth (already solved) --> we have explored the stability of noise, number of powered cells,
##     number of unique events, etc. at different nrpbpcc to justify the threshold applied for retaining cells
##  3) Technical noise --> we removed noisy cells based on bin-level reads overdispersion (as in Minussi et al. 2021)
##     Here we explored the impact of this threshold
##  4) Overfragmentation --> we excluded overfragmented cells (as in Minussi et a. 2021). Here we explore the impact
##     of this threshold
##  5) Null distribution for excluding FP unique events --> we also explore here the impact of not excluding
##     FP unique events in detecting ongoing CIN.
##
## I explore threshold for (3), (4) and (5) using our DLP data from the BRCA2-/- model,
## with profiles generated at 100kb and 30nrpbpcc.
##
#################################################################################


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
MISC_DIR=dirname(this.path())
INPUT_DIR=file.path(BASE,"data")
# Create directories for outputs
OUTPUT_DIR=file.path(MISC_DIR,"outputs")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
PLOTS_DIR=file.path(MISC_DIR,"plots")
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(PLOTS_DIR,"nrpbpcc"), showWarnings = FALSE, recursive = TRUE)

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
source(file.path(BASE, "scripts/uniqueEvents_functions.R"))

# For evaluation of observed noise
sdDiffTrim <- function(x, ..., trim=0.001, scale=TRUE) {
    if (scale)
        x <- x / mean(x, na.rm=TRUE)
    matrixStats::sdDiff(x, ..., trim=trim)
}
sdFUN <- match.fun(sdDiffTrim)

### PARAMETERS
IGNOREDELS=FALSE
noise.thresh=c(0.20, 0.22, 0.28, 0.50) #from "03_1a_BRCA_KO_Explore_Noise.R" --> median, mean+2sd, q75, >max
breaks.thresh=c(45, 64, 110, 250) #from "03_1a_BRCA_KO_Explore_Noise.R" --> median, mean, mean+2sd, >max

### LOAD INPUT DATA
# Profiles at 30 nrpbpcc
sc.files=dir(paste0(INPUT_DIR,"/KOs/nrpbpcc/30nrpbpcc"),full=T,pattern="absoluteCN")


################################################################################
#### 1. APPLY NOISE & BREAKS THRESHOLDS ########################################
# Get single-cell profiles for different combinations of noise & breaks thresholds

SCELLS.CN=as.data.frame(do.call(rbind,lapply(sc.files,function(thisFile){
    print(basename(thisFile))
    
    ## Step 1: Process data
    # read file
    cnobject=readRDS(thisFile)
    
    # combination of noise thresholds
    sc.segsTabs=c()
    for(n in noise.thresh){
        NOISE=n
        print(NOISE)
        
        # tresh for dividing cells
        sc.segsTab=c()
        for(b in breaks.thresh){
            BREAKS=b
            print(BREAKS)
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
                        if(nrow(segTab)<BREAKS){
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
                    segTabs$breaks=BREAKS
                    segTabs$noise=NOISE
                    sc.segsTab=rbind(sc.segsTab,segTabs)
                }
            }
        }
        sc.segsTabs=rbind(sc.segsTabs,sc.segsTab)
    }
    return(sc.segsTabs)
})))

SCELLS.CN$noise=factor(SCELLS.CN$noise, levels = noise.thresh, labels = c("median","mean+2sd","q75","none"))
SCELLS.CN$breaks=factor(SCELLS.CN$breaks, levels = breaks.thresh, labels = c("median","mean","mean+2sd","none"))
SCELLS.CN$combination=paste0(SCELLS.CN$noise,"&",SCELLS.CN$breaks)

# how many cells per combination?
stats = as.data.frame(do.call(rbind,lapply(unique(SCELLS.CN$combination), function(thisComb){
    cn=SCELLS.CN[SCELLS.CN$combination==thisComb,]
    n=length(unique(cn$sample))
    return(c(thisComb,n,as.character(unique(cn$noise)),as.character(unique(cn$breaks))))
})))
colnames(stats)=c("combination","ncells","noise","breaks")
stats$combination=factor(stats$combination,levels = unique(SCELLS.CN$combination))
stats$noise=factor(stats$noise, levels = c("median","mean+2sd","q75","none"))
stats$breaks=factor(stats$breaks, levels = c("median","mean","mean+2sd","none"))
stats$ncells=as.numeric(stats$ncells)

pdf(file=paste0(PLOTS_DIR, "/nrpbpcc/BoxPlot_thresholds.pdf"), width = 100/25.4, height = 100/25.4)
p1=ggplot(stats,aes(x=noise,y=ncells))+
    facet_wrap(~breaks, ncol = 1) +
    geom_bar(stat = "identity")+
    geom_text(aes(label=ncells), vjust=0, size = 1.8) +
    labs(subtitle = "Number of PASS cells per breaks threshold",
         x = "Noise threshold", y = "Number of cells") +
    theme(axis.text.x = element_text(size = 5, angle = 0, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 6),
          legend.position = "none", legend.title = element_blank())
p2=ggplot(stats,aes(x=breaks,y=ncells))+
    facet_wrap(~noise, ncol = 1) +
    geom_bar(stat = "identity")+
    geom_text(aes(label=ncells), vjust=0, size = 1.8) +
    labs(subtitle = "Number of PASS cells per noise threshold",
         x = "Breaks threshold", y = "Number of cells") +
    theme(axis.text.x = element_text(size = 5, angle = 0, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 6),
          legend.position = "none", legend.title = element_blank())
p=ggarrange(p1,p2,ncol=2)
print(p)
dev.off()


################################################################################
#### 2. GET UNIQUE EVENTS IN CELLS ACROSS THRESHOLD COMBINATIONS ###############

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
THRESH=1.00 #% for defining a shared event as recent
SIGN.THRESH=0.05 #for defining significant unique segments

# Mixmodels
MIXMODELS=readRDS(paste0(BASE,"/02_Derive_New_scCINSignatures/data/SC_MIXMODELS.rds"))
# Null distribution matrix => for unique calling
lDistributions=readRDS(paste0(BASE, "/04_Process_scWGS/intermediate/confidence/NullDistribution_Confidence_Segsize.rds"))


### Process for each sequencing depth
dfUnique=c()
dfUnique_mb=c()
CN.EVENTS=c()

for(thisThresh in unique(SCELLS.CN$combination)){
    print(paste0("Thresh combination: ",thisThresh))
    
    # 1. Get profiles of selected cells
    cnprofiles=SCELLS.CN[SCELLS.CN$combination==thisThresh,1:5]
    cnprofiles$sample=factor(cnprofiles$sample, levels=unique(cnprofiles$sample))
    cells=as.character(unique(cnprofiles$sample))
    
    # 2. Get raw values per bin
    cnobject=readRDS(sc.files)
    cnobject=cnobject[cells]
    
    bintable=as.data.frame(do.call(rbind,lapply(1:length(cnobject),function(c){
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
    colnames(bintable)=c("chromosome","start","end","segVal","binVal","sample")
    bintable=bintable[!is.na(bintable$segVal),]
    bintable[,c(2:5)]=sapply(2:5,function(x)as.numeric(bintable[,x]))
    row.names(bintable)=NULL
    
    # 3. Get confidence matrix
    Confidence=as.data.frame(do.call(rbind,lapply(cells, function(thisCell){
        print(paste0("Cell ID: ",thisCell))
        tab=cnprofiles[cnprofiles$sample==thisCell,]
        
        # get confidence scores
        confis = as.data.frame(do.call(rbind,lapply(1:nrow(tab), function(thisSeg){
            # get genomic location of segment
            segs=tab[thisSeg,]
            segments.gr = GRanges(seqnames = segs$chromosome, IRanges(start = segs$start, end = segs$end))
            
            # get raw value in specific cell
            bin=bintable[bintable$sample%in%segs$sample & bintable$chromosome%in%segs$chromosome,]
            bins.gr=GRanges(seqnames = bin$chromosome, IRanges(start = bin$start, end = bin$end))
            overlaps <- findOverlaps(segments.gr, bins.gr)
            rawVal=mean(bin[overlaps@to,4])
            
            # get raw values in the other cells
            bin=bintable[!bintable$sample%in%segs$sample & bintable$chromosome%in%segs$chromosome,]
            bins.gr=GRanges(seqnames = bin$chromosome, IRanges(start = bin$start, end = bin$end))
            overlaps <- findOverlaps(segments.gr, bins.gr)
            oths=bin[overlaps@to,]
            oths <- do.call(c,lapply(unique(oths$sample), function(thisOth){
                bins=mean(oths[oths$sample==thisOth,4])
                return(bins)
            }))
            
            # confidence score => mode of difference with value in other cells
            confi=abs(rawVal-round(oths))
            min=min(confi)
            confi=as.numeric(names(sort(-table(confi))[1]))
            
            # output
            segs$confidence=confi
            segs$minchange=min
            return(segs)
        })))
        colnames(confis)=c("chromosome","start","end","segVal","sample","confidence","minchange")
        return(confis)
    })))
    
    # 4. Extract copy number features
    cnprofiles=cnprofiles[,1:5]
    lProfiles=split(cnprofiles, cnprofiles$sample)
    ECNF=extractThreeCopynumberFeatures.sc(lProfiles, cores=CORES, prePath=file.path(BASE,"data"), rmNorm=RMNORM, size=SIZE, res=BIN.SIZE)
    
    # 5. Get ExC matrix
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
    
    # 6. Prepare events table --> for later!
    cnevents=SCELLS.CN[SCELLS.CN$combination==thisThresh & SCELLS.CN$sample%in%cells & SCELLS.CN$segVal!=2,]
    cnevents$segVal[cnevents$segVal<0]=0
    cnevents$length=cnevents$end-cnevents$start
    cnevents=cnevents[cnevents$length!=0,]
    cnevents$sample=factor(cnevents$sample, levels=unique(cnevents$sample))
    cnevents=cbind(cnevents,ExC)
    CN.EVENTS=rbind(CN.EVENTS,cnevents)
    
    ### 100kb resolution    
    # 7. Classify events
    unify=unifySegments.multiple(profiles=cnprofiles,window=WINDOW,bin.size=BIN.SIZE)
    unify[,4:ncol(unify)]=sapply(4:ncol(unify),function(x)round(unify[,x]))
    events=rbindlist(classEvents.sc(unify,threshold=THRESH),idcol = T) # get unique & shared events
    events=as.data.frame(cbind(events[,c(2:5,1,6)]))
    colnames(events)[5]="sample"
    events=events[events$segVal!=2,] # remove normal segments
    
    # 8. Get unique events
    unique=events[events$class=="unique",]
    unique.gr <- GRanges(seqnames = unique$chromosome, IRanges(start = unique$start, end = unique$end))
    segments.gr = GRanges(seqnames = cnprofiles$chromosome, IRanges(start = cnprofiles$start, end = cnprofiles$end))
    overlaps <- findOverlaps(unique.gr, segments.gr)
    # Get segments that map with the unique event in all samples
    seg.unique <- do.call(rbind,lapply(1:nrow(unique), function(thisUnq){
        # get the index of the events of all clones in the same location of unique event 
        idxs=overlaps@to[overlaps@from == thisUnq]
        segs=cnprofiles[idxs,]
        samp=unique[thisUnq,5]
        usegs=segs[segs$sample==samp,]
        return(usegs)
    }))
    #remove duplicates => unique events can map in the same segment
    seg.unique=seg.unique[!duplicated(seg.unique),]
    seg.unique=seg.unique[order(seg.unique$chromosome,seg.unique$start,seg.unique$end),]
    
    # Check that unique segments are really unique
    confis=Confidence
    seg.unique=left_join(seg.unique,confis,by=c("chromosome","start","end","segVal","sample"))
    
    # know segsize category to compare to null distribution
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
    
    # Significant unique segments
    seg.unique$t.qvalue=p.adjust(as.numeric(seg.unique$t.pvalue), method = "BH")
    seg.unique$significant=ifelse(seg.unique$t.qvalue<SIGN.THRESH,"yes","no")
    seg.unique$significant[seg.unique$minchange<0.7]="no" #minimal difference with other cell may be >0.5 to be considered as real
    seg.unique$combination=thisThresh
    
    # output
    dfUnique=rbind(dfUnique,seg.unique)
    
    ### 1Mb resolution    
    # 7. Classify events
    unify=unifySegments.multiple(profiles=cnprofiles,window=1000e3,bin.size=1000e3)
    unify[,4:ncol(unify)]=sapply(4:ncol(unify),function(x)round(unify[,x]))
    events=rbindlist(classEvents.sc(unify,threshold=THRESH),idcol = T) # get unique & shared events
    events=as.data.frame(cbind(events[,c(2:5,1,6)]))
    colnames(events)[5]="sample"
    events=events[events$segVal!=2,] # remove normal segments
    
    # 8. Get unique events
    unique=events[events$class=="unique",]
    unique.gr <- GRanges(seqnames = unique$chromosome, IRanges(start = unique$start, end = unique$end))
    segments.gr = GRanges(seqnames = cnprofiles$chromosome, IRanges(start = cnprofiles$start, end = cnprofiles$end))
    overlaps <- findOverlaps(unique.gr, segments.gr)
    # Get segments that map with the unique event in all samples
    seg.unique <- do.call(rbind,lapply(1:nrow(unique), function(thisUnq){
        # get the index of the events of all clones in the same location of unique event 
        idxs=overlaps@to[overlaps@from == thisUnq]
        segs=cnprofiles[idxs,]
        samp=unique[thisUnq,5]
        usegs=segs[segs$sample==samp,]
        return(usegs)
    }))
    #remove duplicates => unique events can map in the same segment
    seg.unique=seg.unique[!duplicated(seg.unique),]
    seg.unique=seg.unique[order(seg.unique$chromosome,seg.unique$start,seg.unique$end),]
    
    # Check that unique segments are really unique
    confis=Confidence
    seg.unique=left_join(seg.unique,confis,by=c("chromosome","start","end","segVal","sample"))
    
    # know segsize category to compare to null distribution
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
    
    # Significant unique segments
    seg.unique$t.qvalue=p.adjust(as.numeric(seg.unique$t.pvalue), method = "BH")
    seg.unique$significant=ifelse(seg.unique$t.qvalue<SIGN.THRESH,"yes","no")
    seg.unique$significant[seg.unique$minchange<0.7]="no" #minimal difference with other cell may be >0.5 to be considered as real
    seg.unique$combination=thisThresh
    
    # output
    dfUnique_mb=rbind(dfUnique_mb,seg.unique)
}

# Stats
nUnique=dfUnique%>%group_by(combination, significant)%>%summarise(n=n())
nUnique$freq=NA
for(thresh in unique(nUnique$combination)){nUnique$freq[nUnique$combination==thresh]=nUnique$n[nUnique$combination==thresh]/sum(nUnique$n[nUnique$combination==thresh])}
nUnique$noise=NA
nUnique$noise[1:nrow(nUnique)]=sapply(1:nrow(nUnique),function(x)strsplit(nUnique$combination, '&')[[x]][1])
nUnique$breaks=NA
nUnique$breaks[1:nrow(nUnique)]=sapply(1:nrow(nUnique),function(x)strsplit(nUnique$combination, '&')[[x]][2])

# Plotting unique events
nUnique$noise=factor(nUnique$noise, levels = c("median","mean+2sd","q75","none"))
nUnique$breaks=factor(nUnique$breaks, levels = c("median","mean","mean+2sd","none"))

pdf(file=paste0(PLOTS_DIR, "/nrpbpcc/BoxPlot_thresholds_nUnique.pdf"), width = 100/25.4, height = 100/25.4)
p1=ggplot(nUnique,aes(x=noise,y=n)) +
    facet_wrap(~breaks, ncol = 1) +
    geom_bar(aes(fill=significant), stat = "identity") +
    scale_fill_manual(values=c("red","seagreen3")) +
    geom_text(aes(label=n), vjust=0, size = 1.8) +
    # scale_y_continuous(breaks = c(seq(0,140,25))) +
    labs(subtitle = "Number of unique from PASS cells per break thresholds",
         x = "Noise threshold", y = "Number of unique events") +
    theme(axis.text.x = element_text(size = 5, angle = 0, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 6),
          legend.position = "none", legend.title = element_blank())
p2=ggplot(nUnique,aes(x=breaks,y=n)) +
    facet_wrap(~noise, ncol = 1) +
    geom_bar(aes(fill=significant), stat = "identity") +
    scale_fill_manual(values=c("red","seagreen3")) +
    geom_text(aes(label=n), vjust=0, size = 1.8) +
    # scale_y_continuous(breaks = c(seq(0,140,25))) +
    labs(subtitle = "Number of unique from PASS cells per noise thresholds",
         x = "Breaks threshold", y = "Number of unique events") +
    theme(axis.text.x = element_text(size = 5, angle = 0, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 6),
          legend.position = "none", legend.title = element_blank())
p=ggarrange(p1,p2,ncol=2,hjust=T)
print(p)
dev.off()


################################################################################
#### 3. MAP UNIQUE EVENTS IN CELLS ACROSS SEQUENCING DEPTHS ####################

# Signature definitions
DEFS=readRDS(file=paste0(BASE,"/02_Derive_New_scCINSignatures/outputs/2_New_scCINSignatures_Definitions_Corrected.rds"))
DEFS=DEFS[row.names(DEFS)!="CX7",]

## Removing false-positive unique events
# Map signatures to unique segments 
dfUnique.filt=dfUnique[dfUnique$significant=="yes",]
UNIQUE=left_join(dfUnique.filt,CN.EVENTS,by=c("chromosome","start","end","segVal","sample","combination"))
SIG.EXPOSURES=NULL
SIG.DEFINITIONS=DEFS

SCORES=as.data.frame(do.call(rbind,lapply(1:nrow(UNIQUE),function(thisEvent){
    probs=UNIQUE[thisEvent,18:54]
    if(is.null(SIG.EXPOSURES)){scores=sapply(1:nrow(SIG.DEFINITIONS),FUN=function(x) sum(SIG.DEFINITIONS[x,]*probs))}
    if(!is.null(SIG.EXPOSURES)){
        scores=sapply(1:nrow(SIG.DEFINITIONS),FUN=function(x) sum(SIG.DEFINITIONS[x,]*probs)*SIG.EXPOSURES[x])}
    scores=scores/sum(scores)
    names(scores)=row.names(SIG.DEFINITIONS)
    max=names(which(scores==max(scores)))
    scores=c(scores,max)
    return(scores)
})))
colnames(SCORES)=c(row.names(SIG.DEFINITIONS),"map.sig")
SCORES[,1:16]=sapply(1:16,function(x)as.numeric(SCORES[,x]))
UNIQUE.SCORES=cbind(UNIQUE,SCORES)

# Stats
UNIQUE.SCORES$map.sig=factor(UNIQUE.SCORES$map.sig,levels = row.names(DEFS))
dfUniqueAssign = UNIQUE.SCORES %>% group_by(combination,map.sig,.drop=F) %>% summarise(n=n())
dfUniqueAssign$freq=NA
for(thresh in unique(dfUniqueAssign$combination)){dfUniqueAssign$freq[dfUniqueAssign$combination==thresh]=dfUniqueAssign$n[dfUniqueAssign$combination==thresh]/sum(dfUniqueAssign$n[dfUniqueAssign$combination==thresh])}
dfUniqueAssign$noise=NA
dfUniqueAssign$noise[1:nrow(dfUniqueAssign)]=sapply(1:nrow(dfUniqueAssign),function(x)strsplit(dfUniqueAssign$combination, '&')[[x]][1])
dfUniqueAssign$breaks=NA
dfUniqueAssign$breaks[1:nrow(dfUniqueAssign)]=sapply(1:nrow(dfUniqueAssign),function(x)strsplit(dfUniqueAssign$combination, '&')[[x]][2])

# Plot signature activities across threshold combinations
dfUniqueAssign$noise=factor(dfUniqueAssign$noise, levels = c("median","mean+2sd","q75","none"))
dfUniqueAssign$breaks=factor(dfUniqueAssign$breaks, levels = c("median","mean","mean+2sd","none"))

sigs_palette=c("#4cbcd1", "#db9540","#b16f45", "#4d75b2","#73792f", "#d35138","#ad5fd0", "#6cb543",
               "#ca48a0", "#3a814f","#6861c4", "#b9b04b","#8b99e4", "#ce82bc","#d95070", "#a34e6a", "#5bbf8a")
names(sigs_palette)=paste0("CX", rep(1:17))

pdf(file=paste0(PLOTS_DIR, "/nrpbpcc/BarPlot_mapUnique_TPUnique_threshs_color.pdf"), width = 100/25.4, height = 100/25.4)
p1 = ggplot(dfUniqueAssign, aes(x=noise, y=freq, fill=map.sig))+
    facet_wrap(~breaks,ncol=1) +
    geom_bar(stat = "identity")+
    ylim(0,1) +
    scale_fill_manual(values = sigs_palette) +
    # coord_flip() +
    labs(subtitle = "TP unique events",
         x = "Noise threshold", y = "Relative activity") +
    theme(axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          legend.position = "none")
p2 = ggplot(dfUniqueAssign, aes(x=breaks, y=freq, fill=map.sig))+
    facet_wrap(~noise,ncol=1) +
    geom_bar(stat = "identity")+
    ylim(0,1) +
    scale_fill_manual(values = sigs_palette) +
    # coord_flip() +
    labs(subtitle = "TP unique events",
         x = "Breaks threshold", y = "Relative activity") +
    theme(axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          legend.position = "none")
p=ggarrange(p1,p2,ncol=2,hjust=T)
print(p)
dev.off()


# Plot similarity across definitions
set.seed(433)
pal = circlize::colorRamp2(c(0, 0.8, 0.85, 1), c("white", "#FFE5E0", "#FFBEB2", "#E44244"))(seq(0,1,by=0.01))

tp53.sig=readRDS(file.path(BASE,"05_Induced_CIN/outputs/2.2.1_TP53_signature_definition_unique.rds"))

mCosine=c()
for(thisThresh in unique(UNIQUE$combination)){
    print(paste0("Thres combination: ",thisThresh))
    
    events=UNIQUE[UNIQUE$combination==thisThresh,]
    sxc=colMeans(events[,18:54])
    
    # Subtract TP53 background 
    brca.sig=(as.numeric(sxc)*0.7198475)-(tp53.sig*(1-0.7198475))
    brca.sig[brca.sig<0]=0
    brca.sig[1:20]=brca.sig[1:20]/sum(brca.sig[1:20])
    brca.sig[21:28]=brca.sig[21:28]/sum(brca.sig[21:28])
    brca.sig[29:37]=brca.sig[29:37]/sum(brca.sig[29:37])
    
    # Cosine similarities
    cosine<-do.call(rbind,lapply(1:nrow(DEFS), function(thisSig){
        cos=lsa::cosine(DEFS[thisSig,1:28], brca.sig[1:28])}))
    row.names(cosine)<-row.names(DEFS)
    mCosine=cbind(mCosine,cosine)
}
colnames(mCosine)=unique(UNIQUE$combination)

pdf(file=paste0(PLOTS_DIR, "/nrpbpcc/Heatmap_CXSigs_thresholds_cosine.pdf"), width = 150/25.4, height = 100/25.4)
p=Heatmap(t(mCosine), col = pal, column_title = "RPE1-TP53+BRCA2 (-TP53)", cluster_rows = TRUE, cluster_columns = FALSE,
          heatmap_legend_param = list(
              title = "Cosine"
          ))
print(p)
dev.off()


################################################################################
#### 4. ONGOING CIN BASED ON "REAL" UNIQUE EVENTS ##############################
# Only for good thresholds

## Removing false-positive unique events --> p-val=0.05
# Map signatures to unique segments 
dfUnique.filt=dfUnique[dfUnique$combination=="mean+2sd&mean+2sd",]
dfUnique.filt=dfUnique.filt[dfUnique.filt$significant=="yes",]
UNIQUE=left_join(dfUnique.filt,CN.EVENTS,by=c("chromosome","start","end","segVal","sample","combination"))
SIG.EXPOSURES=NULL
SIG.DEFINITIONS=DEFS

SCORES=as.data.frame(do.call(rbind,lapply(1:nrow(UNIQUE),function(thisEvent){
    probs=UNIQUE[thisEvent,18:54]
    if(is.null(SIG.EXPOSURES)){scores=sapply(1:nrow(SIG.DEFINITIONS),FUN=function(x) sum(SIG.DEFINITIONS[x,]*probs))}
    if(!is.null(SIG.EXPOSURES)){
        scores=sapply(1:nrow(SIG.DEFINITIONS),FUN=function(x) sum(SIG.DEFINITIONS[x,]*probs)*SIG.EXPOSURES[x])}
    scores=scores/sum(scores)
    names(scores)=row.names(SIG.DEFINITIONS)
    max=names(which(scores==max(scores)))
    scores=c(scores,max)
    return(scores)
})))
colnames(SCORES)=c(row.names(SIG.DEFINITIONS),"map.sig")
SCORES[,1:16]=sapply(1:16,function(x)as.numeric(SCORES[,x]))
UNIQUE.SCORES=cbind(UNIQUE,SCORES)


# Plot similarity across definitions
UNIQUE.SCORES$map.sig=factor(UNIQUE.SCORES$map.sig,levels = row.names(DEFS))
dfUniqueAssign = UNIQUE.SCORES %>% group_by(map.sig,.drop=F) %>% summarise(n=n())
dfUniqueAssign$thresh="P<0.05"
dfUniqueAssign$freq=dfUniqueAssign$n/sum(dfUniqueAssign$n)

pdf(file=paste0(PLOTS_DIR, "/nrpbpcc/BarPlot_mapUnique_TPUnique.pdf"), width = 90/25.4, height = 45/25.4)
p = ggplot(dfUniqueAssign,aes(x=map.sig,y=freq))+
    # facet_wrap(~map.sig, ncol=4) +
    ylim(0, 0.4) +
    geom_bar(stat = "identity", position = position_dodge())+
    geom_text(aes(label=n), vjust=0, size = 1.8) +
    labs(subtitle = "Signatures assigned to unique segments",
         x = "CIN signatures", y = "Frequency of segments (%)") +
    theme(axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          legend.position = "none")
print(p)
dev.off()

pdf(file=paste0(PLOTS_DIR, "/nrpbpcc/BarPlot_mapUnique_TPUnique_color.pdf"), width = 20/25.4, height = 45/25.4)
p = ggplot(dfUniqueAssign, aes(x=thresh, y=freq, fill=map.sig))+
    geom_bar(stat = "identity")+
    ylim(0,1) +
    scale_fill_manual(values = sigs_palette) +
    # coord_flip() +
    labs(subtitle = "",
         x = "", y = "Relative activity") +
    theme(axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          legend.position = "none")
print(p)
dev.off()


## Removing false-positive unique events --> p-val=0.1
# Map signatures to unique segments 
dfUnique.filt=dfUnique[dfUnique$combination=="mean+2sd&mean+2sd",]
dfUnique.filt=dfUnique.filt[dfUnique.filt$t.qvalue<=0.1,]
UNIQUE=left_join(dfUnique.filt,CN.EVENTS,by=c("chromosome","start","end","segVal","sample","combination"))
SIG.EXPOSURES=NULL
SIG.DEFINITIONS=DEFS

SCORES=as.data.frame(do.call(rbind,lapply(1:nrow(UNIQUE),function(thisEvent){
    probs=UNIQUE[thisEvent,18:54]
    if(is.null(SIG.EXPOSURES)){scores=sapply(1:nrow(SIG.DEFINITIONS),FUN=function(x) sum(SIG.DEFINITIONS[x,]*probs))}
    if(!is.null(SIG.EXPOSURES)){
        scores=sapply(1:nrow(SIG.DEFINITIONS),FUN=function(x) sum(SIG.DEFINITIONS[x,]*probs)*SIG.EXPOSURES[x])}
    scores=scores/sum(scores)
    names(scores)=row.names(SIG.DEFINITIONS)
    max=names(which(scores==max(scores)))
    scores=c(scores,max)
    return(scores)
})))
colnames(SCORES)=c(row.names(SIG.DEFINITIONS),"map.sig")
SCORES[,1:16]=sapply(1:16,function(x)as.numeric(SCORES[,x]))
UNIQUE.SCORES=cbind(UNIQUE,SCORES)


# Plot similarity across definitions
UNIQUE.SCORES$map.sig=factor(UNIQUE.SCORES$map.sig,levels = row.names(DEFS))
dfUniqueAssign = UNIQUE.SCORES %>% group_by(map.sig,.drop=F) %>% summarise(n=n())
dfUniqueAssign$thresh="P<0.1"
dfUniqueAssign$freq=dfUniqueAssign$n/sum(dfUniqueAssign$n)

pdf(file=paste0(PLOTS_DIR, "/nrpbpcc/BarPlot_mapUnique_TPUnique_0.1.pdf"), width = 90/25.4, height = 45/25.4)
p = ggplot(dfUniqueAssign,aes(x=map.sig,y=freq))+
    # facet_wrap(~map.sig, ncol=4) +
    ylim(0, 0.4) +
    geom_bar(stat = "identity", position = position_dodge())+
    geom_text(aes(label=n), vjust=0, size = 1.8) +
    labs(subtitle = "Signatures assigned to unique segments",
         x = "CIN signatures", y = "Frequency of segments (%)") +
    theme(axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          legend.position = "none")
print(p)
dev.off()

pdf(file=paste0(PLOTS_DIR, "/nrpbpcc/BarPlot_mapUnique_TPUnique_0.1_color.pdf"), width = 20/25.4, height = 45/25.4)
p = ggplot(dfUniqueAssign, aes(x=thresh, y=freq, fill=map.sig))+
    geom_bar(stat = "identity")+
    ylim(0,1) +
    scale_fill_manual(values = sigs_palette) +
    # coord_flip() +
    labs(subtitle = "",
         x = "", y = "Relative activity") +
    theme(axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          legend.position = "none")
print(p)
dev.off()


## Including false-positive unique events
# Map signatures to unique segments 
dfUnique.filt=dfUnique[dfUnique$combination=="mean+2sd&mean+2sd",]
UNIQUE=left_join(dfUnique.filt,CN.EVENTS,by=c("chromosome","start","end","segVal","sample","combination"))
SIG.EXPOSURES=NULL
SIG.DEFINITIONS=DEFS

SCORES=as.data.frame(do.call(rbind,lapply(1:nrow(UNIQUE),function(thisEvent){
    probs=UNIQUE[thisEvent,18:54]
    if(is.null(SIG.EXPOSURES)){scores=sapply(1:nrow(SIG.DEFINITIONS),FUN=function(x) sum(SIG.DEFINITIONS[x,]*probs))}
    if(!is.null(SIG.EXPOSURES)){
        scores=sapply(1:nrow(SIG.DEFINITIONS),FUN=function(x) sum(SIG.DEFINITIONS[x,]*probs)*SIG.EXPOSURES[x])}
    scores=scores/sum(scores)
    names(scores)=row.names(SIG.DEFINITIONS)
    max=names(which(scores==max(scores)))
    scores=c(scores,max)
    return(scores)
})))
colnames(SCORES)=c(row.names(SIG.DEFINITIONS),"map.sig")
SCORES[,1:16]=sapply(1:16,function(x)as.numeric(SCORES[,x]))
UNIQUE.SCORES=cbind(UNIQUE,SCORES)


# Plot similarity across definitions
UNIQUE.SCORES$map.sig=factor(UNIQUE.SCORES$map.sig,levels = row.names(DEFS))
dfUniqueAssign = UNIQUE.SCORES %>% group_by(map.sig,.drop=F) %>% summarise(n=n())
dfUniqueAssign$thresh="All"
dfUniqueAssign$freq=dfUniqueAssign$n/sum(dfUniqueAssign$n)

pdf(file=paste0(PLOTS_DIR, "/nrpbpcc/BarPlot_mapUnique_FPUnique.pdf"), width = 90/25.4, height = 45/25.4)
p = ggplot(dfUniqueAssign,aes(x=map.sig,y=freq))+
    # facet_wrap(~map.sig, ncol=4) +
    ylim(0, 0.4) +
    geom_bar(stat = "identity", position = position_dodge())+
    geom_text(aes(label=n), vjust=0, size = 1.8) +
    labs(subtitle = "Signatures assigned to unique segments",
         x = "CIN signatures", y = "Frequency of segments (%)") +
    theme(axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          legend.position = "none")
print(p)
dev.off()

pdf(file=paste0(PLOTS_DIR, "/nrpbpcc/BarPlot_mapUnique_FPUnique_color.pdf"), width = 20/25.4, height = 45/25.4)
p = ggplot(dfUniqueAssign, aes(x=thresh, y=freq, fill=map.sig))+
    geom_bar(stat = "identity")+
    ylim(0,1) +
    scale_fill_manual(values = sigs_palette) +
    # coord_flip() +
    labs(subtitle = "",
         x = "", y = "Relative activity") +
    theme(axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          legend.position = "none")
print(p)
dev.off()


################################################################################
#### 5. ONGOING CIN BASED ON "REAL" UNIQUE EVENTS AT 1MB #######################
# Only for good thresholds

## Removing false-positive unique events --> p-val=0.05
# Map signatures to unique segments 
dfUnique.filt=dfUnique_mb[dfUnique_mb$combination=="mean+2sd&mean+2sd",]
dfUnique.filt=dfUnique.filt[dfUnique.filt$significant=="yes",]
UNIQUE=left_join(dfUnique.filt,CN.EVENTS,by=c("chromosome","start","end","segVal","sample","combination"))
SIG.EXPOSURES=NULL
SIG.DEFINITIONS=DEFS

SCORES=as.data.frame(do.call(rbind,lapply(1:nrow(UNIQUE),function(thisEvent){
    probs=UNIQUE[thisEvent,18:54]
    if(is.null(SIG.EXPOSURES)){scores=sapply(1:nrow(SIG.DEFINITIONS),FUN=function(x) sum(SIG.DEFINITIONS[x,]*probs))}
    if(!is.null(SIG.EXPOSURES)){
        scores=sapply(1:nrow(SIG.DEFINITIONS),FUN=function(x) sum(SIG.DEFINITIONS[x,]*probs)*SIG.EXPOSURES[x])}
    scores=scores/sum(scores)
    names(scores)=row.names(SIG.DEFINITIONS)
    max=names(which(scores==max(scores)))
    scores=c(scores,max)
    return(scores)
})))
colnames(SCORES)=c(row.names(SIG.DEFINITIONS),"map.sig")
SCORES[,1:16]=sapply(1:16,function(x)as.numeric(SCORES[,x]))
UNIQUE.SCORES=cbind(UNIQUE,SCORES)


# Plot similarity across definitions
UNIQUE.SCORES$map.sig=factor(UNIQUE.SCORES$map.sig,levels = row.names(DEFS))
dfUniqueAssign = UNIQUE.SCORES %>% group_by(map.sig,.drop=F) %>% summarise(n=n())
dfUniqueAssign$thresh="P<0.05"
dfUniqueAssign$freq=dfUniqueAssign$n/sum(dfUniqueAssign$n)

pdf(file=paste0(PLOTS_DIR, "/nrpbpcc/BarPlot_mapUnique_TPUnique_mb.pdf"), width = 90/25.4, height = 45/25.4)
p = ggplot(dfUniqueAssign,aes(x=map.sig,y=freq))+
    # facet_wrap(~map.sig, ncol=4) +
    ylim(0, 0.4) +
    geom_bar(stat = "identity", position = position_dodge())+
    geom_text(aes(label=n), vjust=0, size = 1.8) +
    labs(subtitle = "Signatures assigned to unique segments",
         x = "CIN signatures", y = "Frequency of segments (%)") +
    theme(axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          legend.position = "none")
print(p)
dev.off()

pdf(file=paste0(PLOTS_DIR, "/nrpbpcc/BarPlot_mapUnique_TPUnique_color_mb.pdf"), width = 20/25.4, height = 45/25.4)
p = ggplot(dfUniqueAssign, aes(x=thresh, y=freq, fill=map.sig))+
    geom_bar(stat = "identity")+
    ylim(0,1) +
    scale_fill_manual(values = sigs_palette) +
    # coord_flip() +
    labs(subtitle = "",
         x = "", y = "Relative activity") +
    theme(axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          legend.position = "none")
print(p)
dev.off()


## Removing false-positive unique events --> p-val=0.1
# Map signatures to unique segments 
dfUnique.filt=dfUnique_mb[dfUnique_mb$combination=="mean+2sd&mean+2sd",]
dfUnique.filt=dfUnique.filt[dfUnique.filt$t.qvalue<=0.1,]
UNIQUE=left_join(dfUnique.filt,CN.EVENTS,by=c("chromosome","start","end","segVal","sample","combination"))
SIG.EXPOSURES=NULL
SIG.DEFINITIONS=DEFS

SCORES=as.data.frame(do.call(rbind,lapply(1:nrow(UNIQUE),function(thisEvent){
    probs=UNIQUE[thisEvent,18:54]
    if(is.null(SIG.EXPOSURES)){scores=sapply(1:nrow(SIG.DEFINITIONS),FUN=function(x) sum(SIG.DEFINITIONS[x,]*probs))}
    if(!is.null(SIG.EXPOSURES)){
        scores=sapply(1:nrow(SIG.DEFINITIONS),FUN=function(x) sum(SIG.DEFINITIONS[x,]*probs)*SIG.EXPOSURES[x])}
    scores=scores/sum(scores)
    names(scores)=row.names(SIG.DEFINITIONS)
    max=names(which(scores==max(scores)))
    scores=c(scores,max)
    return(scores)
})))
colnames(SCORES)=c(row.names(SIG.DEFINITIONS),"map.sig")
SCORES[,1:16]=sapply(1:16,function(x)as.numeric(SCORES[,x]))
UNIQUE.SCORES=cbind(UNIQUE,SCORES)


# Plot similarity across definitions
UNIQUE.SCORES$map.sig=factor(UNIQUE.SCORES$map.sig,levels = row.names(DEFS))
dfUniqueAssign = UNIQUE.SCORES %>% group_by(map.sig,.drop=F) %>% summarise(n=n())
dfUniqueAssign$thresh="P<0.1"
dfUniqueAssign$freq=dfUniqueAssign$n/sum(dfUniqueAssign$n)

pdf(file=paste0(PLOTS_DIR, "/nrpbpcc/BarPlot_mapUnique_TPUnique_0.1_mb.pdf"), width = 90/25.4, height = 45/25.4)
p = ggplot(dfUniqueAssign,aes(x=map.sig,y=freq))+
    # facet_wrap(~map.sig, ncol=4) +
    ylim(0, 0.4) +
    geom_bar(stat = "identity", position = position_dodge())+
    geom_text(aes(label=n), vjust=0, size = 1.8) +
    labs(subtitle = "Signatures assigned to unique segments",
         x = "CIN signatures", y = "Frequency of segments (%)") +
    theme(axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          legend.position = "none")
print(p)
dev.off()

pdf(file=paste0(PLOTS_DIR, "/nrpbpcc/BarPlot_mapUnique_TPUnique_0.1_color_mb.pdf"), width = 20/25.4, height = 45/25.4)
p = ggplot(dfUniqueAssign, aes(x=thresh, y=freq, fill=map.sig))+
    geom_bar(stat = "identity")+
    ylim(0,1) +
    scale_fill_manual(values = sigs_palette) +
    # coord_flip() +
    labs(subtitle = "",
         x = "", y = "Relative activity") +
    theme(axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          legend.position = "none")
print(p)
dev.off()


## Including false-positive unique events
# Map signatures to unique segments 
dfUnique.filt=dfUnique_mb[dfUnique_mb$combination=="mean+2sd&mean+2sd",]
UNIQUE=left_join(dfUnique.filt,CN.EVENTS,by=c("chromosome","start","end","segVal","sample","combination"))
SIG.EXPOSURES=NULL
SIG.DEFINITIONS=DEFS

SCORES=as.data.frame(do.call(rbind,lapply(1:nrow(UNIQUE),function(thisEvent){
    probs=UNIQUE[thisEvent,18:54]
    if(is.null(SIG.EXPOSURES)){scores=sapply(1:nrow(SIG.DEFINITIONS),FUN=function(x) sum(SIG.DEFINITIONS[x,]*probs))}
    if(!is.null(SIG.EXPOSURES)){
        scores=sapply(1:nrow(SIG.DEFINITIONS),FUN=function(x) sum(SIG.DEFINITIONS[x,]*probs)*SIG.EXPOSURES[x])}
    scores=scores/sum(scores)
    names(scores)=row.names(SIG.DEFINITIONS)
    max=names(which(scores==max(scores)))
    scores=c(scores,max)
    return(scores)
})))
colnames(SCORES)=c(row.names(SIG.DEFINITIONS),"map.sig")
SCORES[,1:16]=sapply(1:16,function(x)as.numeric(SCORES[,x]))
UNIQUE.SCORES=cbind(UNIQUE,SCORES)


# Plot similarity across definitions
UNIQUE.SCORES$map.sig=factor(UNIQUE.SCORES$map.sig,levels = row.names(DEFS))
dfUniqueAssign = UNIQUE.SCORES %>% group_by(map.sig,.drop=F) %>% summarise(n=n())
dfUniqueAssign$thresh="All"
dfUniqueAssign$freq=dfUniqueAssign$n/sum(dfUniqueAssign$n)

pdf(file=paste0(PLOTS_DIR, "/nrpbpcc/BarPlot_mapUnique_FPUnique_mb.pdf"), width = 90/25.4, height = 45/25.4)
p = ggplot(dfUniqueAssign,aes(x=map.sig,y=freq))+
    # facet_wrap(~map.sig, ncol=4) +
    ylim(0, 0.4) +
    geom_bar(stat = "identity", position = position_dodge())+
    geom_text(aes(label=n), vjust=0, size = 1.8) +
    labs(subtitle = "Signatures assigned to unique segments",
         x = "CIN signatures", y = "Frequency of segments (%)") +
    theme(axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          legend.position = "none")
print(p)
dev.off()

pdf(file=paste0(PLOTS_DIR, "/nrpbpcc/BarPlot_mapUnique_FPUnique_color_mb.pdf"), width = 20/25.4, height = 45/25.4)
p = ggplot(dfUniqueAssign, aes(x=thresh, y=freq, fill=map.sig))+
    geom_bar(stat = "identity")+
    ylim(0,1) +
    scale_fill_manual(values = sigs_palette) +
    # coord_flip() +
    labs(subtitle = "",
         x = "", y = "Relative activity") +
    theme(axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          legend.position = "none")
print(p)
dev.off()