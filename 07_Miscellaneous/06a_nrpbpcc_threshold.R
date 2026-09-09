################################################################################
## This script is for evaluating different nrpbpcc for excluding unpowered cells
## Started: 16/04/2026
## Updated: 16/04/2026
##
## I donwsampled cells at 10, 15, 30, 45 and 60 nrpbpcc
## Here, I compare the number of powered cells, the noise across depths, the
## number of unique events per cell, and the ongoing signature extracted.
##
## Script for the near-diploid BRCA2-deficient model
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
depths=c("10","15","30","45","60")

### LOAD INPUT DATA
# Profiles at all nrpbpcc
sc.files=dir(paste0(INPUT_DIR,"/KOs/nrpbpcc"),full=T,pattern="nrpbpcc")
sc.files=do.call(c,lapply(sc.files, function(path){p=(file.path(path, "absoluteCN.rds"))}))


################################################################################
#### 1. GET NUMBER & NOISE OF CELLS FROM SINGLE-CELL PROFILES ##################

### 0. Extract profiles
SCELLS.CN=as.data.frame(do.call(rbind,lapply(1:length(sc.files),function(f){
    thisFile=sc.files[f]
    
    # read file
    cnobject=readRDS(thisFile)
    
    # get profiles
    if(length(cnobject)>0){
        segTabs=c()
        for(c in 1:length(cnobject)){
            cn=cnobject[[c]]
            copynumber=assayDataElement(cn, "copynumber")
            thisClon=colnames(copynumber)
            segTab=CNpare::getSegTable(cn[,thisClon])
            segTab$sample=thisClon
            segTab[,2:4]=sapply(2:4,function(x)as.numeric(segTab[,x]))
            segTabs=rbind(segTabs,segTab)
        }
    }
    
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
        segTabs$nrpbpcc=depths[f]
        return(segTabs) 
    }
})))


### 1. Get number of cells across depths
dfNCells=as.data.frame(do.call(rbind,lapply(1:length(sc.files),function(f){
    thisFile=sc.files[f]
    
    # read file
    cnobject=readRDS(thisFile)
    
    # ncells
    ncells=as.numeric(length(cnobject))
    
    # output
    out=c(depths[f],ncells)
})))
colnames(dfNCells)=c("nrpbpcc","ncells")
dfNCells$ncells=as.numeric(dfNCells$ncells)
dfNCells$freq=dfNCells$ncells/96

# Plotting number of cells
pdf(file=paste0(PLOTS_DIR, "/nrpbpcc/BarPlot_nCells_nrpbpcc.pdf"), width = 45/25.4, height = 45/25.4)
p=ggplot(dfNCells,aes(x=nrpbpcc,y=ncells))+
    geom_bar(stat = "identity")+
    labs(subtitle = "Number of powered cells",
         x = "Number of reads per bin\nper chromosome copy", y = "Number of cells") +
    theme(axis.text.x = element_text(size = 5, angle = 0, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 6),
          legend.position = "none", legend.title = element_blank())
print(p)
dev.off()


### 2. Get noise of cells across depths
dfNoise=as.data.frame(do.call(rbind,lapply(1:length(sc.files),function(f){
    print(f)
    thisFile=sc.files[f]
    
    # read file
    cnobject=readRDS(thisFile)
    
    # noise
    if(length(cnobject)>0){
        out=c()
        for(c in 1:length(cnobject)){
            cn=cnobject[[c]]
            copynumber=assayDataElement(cn, "copynumber")
            noise=apply(copynumber, 2, sdFUN, na.rm=TRUE)
            out=rbind(out,c(noise,gsub(".markdup","",names(noise)),depths[f]))
            
        }
        out=as.data.frame(out)
        colnames(out)=c("noise","cell","nrpbpcc")
        return(out)
    }
})))
dfNoise$cell[1:nrow(dfNoise)]=paste0("DLP22_",sapply(1:nrow(dfNoise),function(x)strsplit(dfNoise$cell, '_')[[x]][2]))
dfNoise$noise=as.numeric(dfNoise$noise)


# Plotting noise
png(file=paste0(PLOTS_DIR, "/nrpbpcc/BoxPlot_noise_BRCA_KOs_all.png"), width = 120/25.4, height = 80/25.4, units = "in", res = 300)
p=ggplot(dfNoise, aes(x=nrpbpcc, y=noise)) + 
    geom_boxplot(outlier.shape = NA, color="grey30", fill=NA)+
    geom_point(pch = 21, size = 1, col = "black", 
               position = position_jitter(w = 0.2, h = 0), alpha = 0.8) +
    labs(subtitle = "Comparing noise across sequencing depths") +
    xlab("Number of reads per bin per chromosome copy") +
    ylab("Noise (observed SD)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
          axis.text.y = element_text(size = 5),
          legend.position = "none")
print(p)
dev.off()

# limit to cells in 60nrpbpcc
cells=unique(SCELLS.CN$sample[SCELLS.CN$nrpbpcc=="60"])
cells[1:length(cells)]=paste0("DLP22_",sapply(1:length(cells),function(x)strsplit(cells, '_')[[x]][2]))

pdf(file=paste0(PLOTS_DIR, "/nrpbpcc/BoxPlot_noise_BRCA_KOs_filtered.pdf"), width = 80/25.4, height = 45/25.4)
p=ggplot(dfNoise[dfNoise$cell%in%cells,], aes(x=nrpbpcc, y=noise)) + 
    geom_hline(yintercept = 0.22, linetype = "dashed", colour = "red", size = 0.5) +
    geom_boxplot(outlier.shape = NA, color="grey30", fill=NA)+
    geom_point(pch = 21, size = 1, col = "black", 
               position = position_jitter(w = 0.2, h = 0), alpha = 0.8) +
    labs(subtitle = "Comparing noise across sequencing depths") +
    xlab("Number of reads per bin per chromosome copy") +
    ylab("Noise (observed SD)") +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 5),
          axis.text.y = element_text(size = 5),
          legend.position = "none")
print(p)
dev.off()

pairwise.wilcox.test(dfNoise$noise[dfNoise$cell%in%cells],dfNoise$nrpbpcc[dfNoise$cell%in%cells], p.adjust.method = "bonferroni")
#    10      15      30    45   
# 15 1.8e-05 -       -     -    
# 30 4.0e-12 2.2e-05 -     -    
# 45 7.3e-14 5.5e-07 0.260 -    
# 60 7.6e-15 1.5e-08 0.026 0.254


################################################################################
#### 1. GET UNIQUE EVENTS IN CELLS ACROSS SEQUENCING DEPTHS ####################

### PARAMETERS
SAVEINTERMEDIATE=TRUE
REMOVEQUIET=TRUE # Just to check that all samples have at least 20 CNAs. If not remove.
WINDOW=100e3 #resolution of single-cell data x2 for calling a segment
NOISE=0.21 #obtained in the "03_1a_BRCA_KO_Explore_Noise.R"
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
# Select cells for comparing --> passing all qc criteria
SEG.TABLE=readRDS(paste0(BASE, "/04_Process_scWGS/intermediate/3.1.1_BRCA_KO_Segments.rds"))
cells=unique(as.character(SEG.TABLE[["scBRCA2_DLP22"]]$sample)) 


### Process for each sequencing depth
dfUnique=c()
CN.EVENTS=c()

for(d in 1:length(depths)){
    thisDepth=depths[d]
    thisFile=sc.files[d]
    
    print(paste0("Nrpbpcc: ",thisDepth))
    
    # 1. Get profiles of selected cells
    cnprofiles=SCELLS.CN[SCELLS.CN$nrpbpcc==thisDepth & SCELLS.CN$sample%in%cells,1:5]
    cnprofiles$sample=factor(cnprofiles$sample, levels=unique(cnprofiles$sample))
    
    # 2. Get raw values per bin
    cnobject=readRDS(thisFile)
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
    cnevents=SCELLS.CN[SCELLS.CN$nrpbpcc==thisDepth & SCELLS.CN$sample%in%cells & SCELLS.CN$segVal!=2,]
    cnevents$segVal[cnevents$segVal<0]=0
    cnevents$length=cnevents$end-cnevents$start
    cnevents=cnevents[cnevents$length!=0,]
    cnevents$sample=factor(cnevents$sample, levels=unique(cnevents$sample))
    cnevents=cbind(cnevents,ExC)
    CN.EVENTS=rbind(CN.EVENTS,cnevents)
        
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
    seg.unique$nrpbpcc=thisDepth
    
    # output
    dfUnique=rbind(dfUnique,seg.unique)
}

# Stats
nUnique=dfUnique%>%group_by(nrpbpcc, significant)%>%summarise(n=n())
nUnique$freq=NA
for(depth in unique(nUnique$nrpbpcc)){nUnique$freq[nUnique$nrpbpcc==depth]=nUnique$n[nUnique$nrpbpcc==depth]/sum(nUnique$n[nUnique$nrpbpcc==depth])}

# Plotting unique events
pdf(file=paste0(PLOTS_DIR, "/nrpbpcc/BarPlot_nUnique_nrpbpcc.pdf"), width = 45/25.4, height = 45/25.4)
p=ggplot(nUnique,aes(x=nrpbpcc,y=n))+
    geom_bar(aes(fill=significant), stat = "identity") +
    scale_fill_manual(values=c("red","seagreen3")) +
    scale_y_continuous(breaks = c(seq(0,140,25))) +
    labs(subtitle = "",
         x = "Number of reads per bin\nper chromosome copy", y = "Number of unique events") +
    theme(axis.text.x = element_text(size = 5, angle = 0, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 6),
          legend.position = "none", legend.title = element_blank())
print(p)
dev.off()


################################################################################
#### 2. MAP UNIQUE EVENTS IN CELLS ACROSS SEQUENCING DEPTHS ####################

# Signature definitions
DEFS=readRDS(file=paste0(BASE,"/02_Derive_New_scCINSignatures/outputs/2_New_scCINSignatures_Definitions_Corrected.rds"))
DEFS=DEFS[row.names(DEFS)!="CX7",]

# Map signatures to unique segments 
dfUnique.filt=dfUnique[dfUnique$significant=="yes",]
UNIQUE=left_join(dfUnique.filt,CN.EVENTS,by=c("chromosome","start","end","segVal","sample","nrpbpcc"))
SIG.EXPOSURES=NULL
SIG.DEFINITIONS=DEFS

SCORES=as.data.frame(do.call(rbind,lapply(1:nrow(UNIQUE),function(thisEvent){
    probs=UNIQUE[thisEvent,16:52]
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


# Plot 
UNIQUE.SCORES$map.sig=factor(UNIQUE.SCORES$map.sig,levels = row.names(DEFS))
dfUniqueAssign = UNIQUE.SCORES %>% group_by(nrpbpcc,map.sig,.drop=F) %>% summarise(n=n())
dfUniqueAssign$freq=NA
for(depth in unique(dfUniqueAssign$nrpbpcc)){dfUniqueAssign$freq[dfUniqueAssign$nrpbpcc==depth]=dfUniqueAssign$n[dfUniqueAssign$nrpbpcc==depth]/sum(dfUniqueAssign$n[dfUniqueAssign$nrpbpcc==depth])}

pdf(file=paste0(PLOTS_DIR, "/nrpbpcc/BarPlot_mapUnique_nrpbpcc.pdf"), width = 120/25.4, height = 120/25.4)
p = ggplot(dfUniqueAssign,aes(x=nrpbpcc,y=freq))+
    facet_wrap(~map.sig, ncol=4) +
    geom_bar(stat = "identity", position = position_dodge())+
    geom_text(aes(label=n), vjust=0, size = 1.8) +
    labs(subtitle = "Signatures assigned to unique segments",
         x = "Number of reads per bin per chromosome copy", y = "Frequency of segments (%)") +
    theme(axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          legend.position = "none")
print(p)
dev.off()


sigs_palette=c("#4cbcd1", "#db9540","#b16f45", "#4d75b2","#73792f", "#d35138","#ad5fd0", "#6cb543",
               "#ca48a0", "#3a814f","#6861c4", "#b9b04b","#8b99e4", "#ce82bc","#d95070", "#a34e6a", "#5bbf8a")
names(sigs_palette)=paste0("CX", rep(1:17))

pdf(file=paste0(PLOTS_DIR, "/nrpbpcc/BarPlot_mapUnique_nrpbpcc_color.pdf"), width = 45/25.4, height = 45/25.4)
p = ggplot(dfUniqueAssign, aes(x=nrpbpcc, y=freq, fill=map.sig))+
    geom_bar(stat = "identity")+
    ylim(0,1) +
    scale_fill_manual(values = sigs_palette) +
    # coord_flip() +
    labs(subtitle = "",
         x = "Number of reads per bin\nper chromosome copy", y = "Relative activity") +
    theme(axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          legend.position = "none")
print(p)
dev.off()


################################################################################
#### COMPARE NOISE WITH ACT DATA ###############################################
Minussi.NOISE=readRDS(paste0(BASE, "/04_Process_scWGS/intermediate/1.0_scMinussi_Powered_Noise.rds"))
Minussi.NOISE=Minussi.NOISE[Minussi.NOISE$sample%in%c("TN1","TN2","TN5","MDAMB231_popp31","MDAMB231c8"),]
Minussi.NOISE$nrpbpcc="15"

data=rbind(dfNoise[dfNoise$nrpbpcc=="30",],Minussi.NOISE[,c(1:2,4)])
pairwise.wilcox.test(data$noise,data$nrpbpcc, p.adjust.method = "bonferroni") #0.059

pdf(file=paste0(PLOTS_DIR, "/nrpbpcc/BoxPlot_noise_Minussi_BRCAmodel.pdf"), width = 45/25.4, height = 45/25.4)
p=ggplot(data, aes(x=nrpbpcc, y=noise)) + 
    # geom_hline(yintercept = 0.22, linetype = "dashed", colour = "red", size = 0.5) +
    geom_boxplot(outlier.shape = NA, color="grey30", fill=NA)+
    geom_point(pch = 21, size = 1, col = "black", 
               position = position_jitter(w = 0.2, h = 0), alpha = 0.8) +
    labs(subtitle = "") +
    xlab("Number of reads per bin\nper chromosome copy") +
    ylab("Noise (observed SD)") +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 5),
          axis.text.y = element_text(size = 5),
          legend.position = "none")
print(p)
dev.off()



