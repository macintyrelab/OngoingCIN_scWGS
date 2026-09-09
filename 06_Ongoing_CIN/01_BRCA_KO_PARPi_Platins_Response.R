################################################################################
## This script is to process bulk and single-cell sequencing data from TP53+BRCA2 KO
## used in the drug screening
## Started: 01/08/24
## Updated: 02/06/26
##
## I derived copy number profiles from scWGS data generated in the lab using DLP+
## We validated that BRCA1/2 KO model acquires copy number alterations induced by
## deficiencies in the HR repair mechanism (see 05_Induced_CIN/02_BRCA_KO_IHR_Signatures.R script)
## 
## Here I compare the ability to predict response to PARPi and platins using bulk 
## vs single-cell sequencing data, and therefore ongoingCIN signatures. We treated
## BRCA1/2 clones the day before starting the treatment. These clones has also been 
## sequenced by bulk sWGS.
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
ONGO_DIR=dirname(this.path())
INPUT_DIR=file.path(SCELL_DIR,"intermediate")
# Create directories for outputs
OUTPUT_DIR=file.path(ONGO_DIR,"outputs/brca_parpi")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
PLOTS_DIR=file.path(ONGO_DIR,"plots/brca_parpi")
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)


### PLOTTING THEME
theme_set(theme_tufte(base_size = 6, base_family = "ArialMT"))
theme_update(text = element_text(size = 6),
             axis.text = element_text(size = 6),
             axis.line = element_line(size = 0.5), 
             axis.ticks = element_line(size = 0.5),
             axis.ticks.length = unit(.1, "cm"))

sigs_palette=c("#4cbcd1", "#db9540","#b16f45", "#4d75b2","#73792f", "#d35138","#ad5fd0", "#6cb543",
               "#ca48a0", "#3a814f","#6861c4", "#b9b04b","#8b99e4", "#ce82bc","#d95070", "#a34e6a", "#5bbf8a")
names(sigs_palette)=paste0("CX", rep(1:17))

### FUNCTIONS
source(file.path(BASE, "scripts/newCNSigs_functions.R"))
source(file.path(BASE, "scripts/main_functions.R"))
source(file.path(BASE, "scripts/plotSCHeatmap.R"))
source(file.path(BASE, "scripts/uniqueEvents_functions.R"))

segTab_plot <- function(dat){
    #prepare length of each chr
    chr_sizes <- CNpare:::chr_sizes[1:22,]
    chr_sizes$length <- 1e-6*chr_sizes$length #Convert sizes to Mb
    chr_sizes$offset <- cumsum(chr_sizes$length)-chr_sizes$length #Calculate the interval range of each chr
    
    #adapt coordinates for plotting
    events<-CNpare::CNconvert(dat,chr_sizes)
    
    #plot
    plot.new()
    plot.window(xlim=c(0,sum(events$length)),ylim=c(0,10),yaxs="i",xaxs="r")
    abline(v=chr_sizes$offset,lty=2,lwd=1)
    abline(v=0,lty=1,lwd=1)
    
    #title
    sample = as.character(unique(dat$sample))
    title(main="Copy Number Profile", sub = sample, xlab="Chromosome", ylab="Copy Number", cex.lab=1.1, cex.main=1.5, font.lab=2, outer=FALSE)
    
    #set the axis
    axis(1, labels=seq_len(22), tick=TRUE, at=chr_sizes$offset+chr_sizes$length/2, cex.axis=1)
    
    lablist<-as.vector(c(0:10))
    text(par("usr")[4], seq(0, 10, by=1), labels=lablist, srt=0, pos=2, xpd=TRUE, offset=1, cex=1.3)
    
    #plot segments
    for (i in seq_len(nrow(events))){
        e <- events[i,]
        segments(e$start,e$segVal,e$end,e$segVal,lwd=3.,col="red")
    }
}



### PARAMETERS
RESOLUTION="100kb"
BIN.SIZE=as.numeric(substr(RESOLUTION,1,nchar(RESOLUTION)-2))*1000
WINDOW=as.numeric(BIN.SIZE*2)
THRESH=1.00 #% for defining a shared event as clonal/subclonal
SIGN.THRESH=0.05 #for defining significant unique segments

RMNORM=TRUE
CORES=1
WIGGLE=0.1
IGNOREDELS=FALSE
SIZE=5000000
UNINFPRIOR=TRUE


### LOAD INPUT DATA
# Profiles 
SCELLS.CN=readRDS(paste0(INPUT_DIR, "/3.1.0_BRCA_KO_Segments.rds")) 
BULK.CN=readRDS(dir(paste0(BASE,"/data/KOs/bulk"),full=T,pattern="absoluteCN"))
# ECNF 
ECNF=readRDS(paste0(INPUT_DIR, "/3.1.0_BRCA_KO_ECNF.rds"))
# ExC matrix
ExC=readRDS(paste0(INPUT_DIR, "/3.1.1_BRCA_KO_ExC.rds"))

# Feature components
MIXMODELS=readRDS(paste0(BASE,"/data/SC_MIXMODELS.rds"))
# Signature definitions
DEFS=readRDS(file=paste0(NEWSIGS_DIR,"/outputs/2_New_scCINSignatures_Definitions_Corrected.rds"))
DEFS=DEFS[row.names(DEFS)!="CX7",]

# Null distribution matrix => for unique calling
lDistributions=readRDS(paste0(INPUT_DIR, "/confidence/NullDistribution_Confidence_Segsize.rds"))


################################################################################
#### 0. PREPARE DATA ###########################################################

# We perform shallow WGS of the parental clone (before treatment with PARPi)
# to derive absolute copy number profile. In addition, we sequenced single cells
# from the same parental clon to compare bulk vs ongoing signal

### 1. PROCESS BULK-DERIVED COPY NUMBER PROFILES
sample=c("ac_GM-589L_S51_L002.markdup","ac_GM-586L_S48_L002.markdup") #bulk of BRCA2 and BRCA1 KO model respectively
bulk.segTab=c()
for(s in sample){
    tab=as.data.frame(CNpare::getSegTable(BULK.CN[,s]))
    tab$sample=gsub(".markdup","",s)
    bulk.segTab=rbind(bulk.segTab,tab)
}
BULK.CN=bulk.segTab
BULK.CN[,2:4]=sapply(2:4,function(x)as.numeric(BULK.CN[,x]))
bulk.segTab=NULL

# avoiding artefacts
BULK.CN$segVal[BULK.CN$segVal<0]=0
BULK.CN=BULK.CN[(BULK.CN$end-BULK.CN$start)!=0,]

# Set everything very close to 2 to 2
BULK.CN$segVal[BULK.CN$segVal > (2-WIGGLE) & BULK.CN$segVal < (2+WIGGLE)] = 2

# Merge segments only when segments with same copy number value follow each other 
df=idSmoothingTargets(BULK.CN, WIGGLE = 0, colNameSegVal = "segVal", colNameChr = "chromosome", IGNOREDELS = IGNOREDELS)
# Split by sample name
lRaw=split(df, df$sample)
# Smooth segments by taking the weighted average of the segVal and their lengths
lSmooth=smoothSegments(lRaw, CORES = 1, SMOOTHINGFACTOR = 0, colNameMerge = "segVal", colNameChr = "chromosome",
                       colNameStart = "start", colNameEnd = "end", IGNOREDELS = IGNOREDELS, asDf = FALSE)

BULK.CN=data.frame(rbindlist(lSmooth))
BULK.CN$model[BULK.CN$sample=="ac_GM-589L_S51_L002"]="BRCA2-bulk"
BULK.CN$model[BULK.CN$sample=="ac_GM-586L_S48_L002"]="BRCA1-bulk"


### 2. PROCESS SINGLE-CELL DATA

## Prepare events table
CN.EVENTS=SCELLS.CN[SCELLS.CN$segVal!=2,]
CN.EVENTS$segVal[CN.EVENTS$segVal<0]=0
CN.EVENTS$length=CN.EVENTS$end-CN.EVENTS$start
CN.EVENTS=CN.EVENTS[CN.EVENTS$length!=0,]
CN.EVENTS=cbind(CN.EVENTS,ExC)

# select the TP53+BRCA1/2 model treated with PARPi
CN.EVENTS=CN.EVENTS[CN.EVENTS$model%in%c("scTP53_drug","scBRCA1_drug","scBRCA2_drug"),]
CN.EVENTS$sample=factor(CN.EVENTS$sample,levels=unique(as.character(CN.EVENTS$sample))) #104 cells

## Classify CNAs as clonal, subclonal & unique
SAMPS=unique(CN.EVENTS$model)
lcntables=lapply(SAMPS,function(thisSamp){
    cn=SCELLS.CN[SCELLS.CN$model==thisSamp,1:5]
    return(cn)
})
names(lcntables)=SAMPS

# Get confidence matrix per segment
confidence.files=c(dir(paste0(INPUT_DIR,"/confidence/brca_ko"),full=T,pattern=SAMPS[1]),
                   dir(paste0(INPUT_DIR,"/confidence/brca_ko"),full=T,pattern=SAMPS[2]),
                   dir(paste0(INPUT_DIR,"/confidence/brca_ko"),full=T,pattern=SAMPS[3]))
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
saveRDS(nUnique.model,file.path(OUTPUT_DIR,"stats_unique_per_model_BRCA_PARPi_Platins.rds"))


################################################################################
#### 1. DECODE CIN IN BULK DATA ################################################

# We quantify signatures using YAPSA when the sample has at least >10-20 CNAS
# However only 3 segments are observed in the bulk profile from RPE1-TP53+BRCA2 model
# & 4 segments are observed in the bulk profile from RPE1-TP53+BRCA1 model.
# Therefore I just can explore the features of the events detected at bulk...
table(BULK.CN$model[BULK.CN$segVal!=2])

## 1.1. Extract copy number features
lProfiles=split(BULK.CN, BULK.CN$sample)
ECNF=extractThreeCopynumberFeatures.sc(lProfiles, cores=CORES, prePath=paste0(BASE,"/01_Derive_New_CINSignatures/data"), rmNorm=RMNORM, size=SIZE, res=BIN.SIZE)

## Get ExC matrix
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
        if(length(id)!=0){postDatScaled[id,]=matrix(rep(c(rep(0,7),1), length(id)), ncol = 8, nrow =length(id), byrow = T)}
    }
    # Create component matrix
    colnames(postDatScaled) = paste0(thisFeature, 1:ncol(postDatScaled))
    postDatScaled=as.matrix(postDatScaled)
    
    return(postDatScaled)
})
names(lMats)=FEATS
ExC.bulk=do.call(cbind, lMats)


## 1.2. Map signatures to non-diploid segments
# To know how similar they are to signatures
sig.scores=as.data.frame(do.call(rbind,lapply(1:nrow(ExC.bulk),function(e){
    probs=ExC.bulk[e,]
    scores=sapply(1:nrow(DEFS),FUN=function(x) sum(DEFS[x,]*probs))
    scores=scores/sum(scores)
    names(scores)=row.names(DEFS)
    max=names(which(scores==max(scores)))
    scores=c(scores,max)
    return(scores)
})))
colnames(sig.scores)[17]="map.sig"
sig.scores[,1:16]=sapply(1:16,function(x)as.numeric(sig.scores[,x]))
unique(sig.scores$map.sig) 
sig.scores$model=BULK.CN$model[BULK.CN$segVal!=2]


# Barplot for manuscript
dat = sig.scores %>% group_by(model,map.sig,.drop=F) %>% summarise(n=n())
dat$freq=1

for(m in unique(sig.scores$model)){
    pdf(file=paste0(PLOTS_DIR, "/Fig3_BarPlot_",m,"_KO_PARPi_Platins.pdf"), width = 30/25.4, height = 20/25.4)
    p = ggplot(dat[dat$model==m,],aes(x=model,y=freq,fill=map.sig))+
        geom_bar(stat = "identity")+
        ylim(0,1) +
        scale_fill_manual(values = sigs_palette) +
        coord_flip() +
        labs(subtitle = "",
             x = "", y = "") +
        theme(axis.text.x = element_text(size = 5),
              axis.text.y = element_blank(),
              legend.position = "none")
    print(p)
    dev.off()
}


# Conclusion! We detect low number of events at bulk level (only 4 and 3) & all events map to CX1 (prob>0.2 (random is 1/17))
# These model represent a premalignant stage of tumor evolution, where we will not
# detect the ongoing accumulation of CNAs due to impaired homologous recombination
# We know that this model respond to PARPi and platins, so it opens a door for earlier interventions


################################################################################
#### 2. DISSECT ONGOING CIN USING SINGLE-CELL DATA #############################

# We already confirmed that the patterns of induced copy number alterations in the
# RPE1-TP53+BRCA1/2 KO model where similar to CX2/CX3 definitions (IHR signatures)

# Here we are going to map unique events (induced) to signatures. Our signature
# mapping method includes a correction by the signature activities found in the
# sample. As we cannot quantify activities at sample level (not enough alterations
# at bulk level), we mapped only to signatures mapped to all cells and events.
# We thus get probabilities for mapping correction by taking into account signature scores 
# in all events of all cells


## 2.1. Detect which signatures are present => mapping all segments in all cells
sig.scores=as.data.frame(do.call(rbind,lapply(1:nrow(CN.EVENTS),function(e){
    probs=CN.EVENTS[e,9:45]
    scores=sapply(1:nrow(DEFS),FUN=function(x) sum(DEFS[x,]*probs))
    scores=scores/sum(scores)
    names(scores)=row.names(DEFS)
    max=names(which(scores==max(scores)))
    scores=c(scores,max)
    return(scores)
})))
colnames(sig.scores)=c(row.names(DEFS),"map.sig")
sig.scores[,1:16]=sapply(1:16,function(x)as.numeric(sig.scores[,x]))

sigs=unique(sig.scores$map.sig) #7 signatures detected
acts=colMeans(sig.scores[1:16])
acts[!names(acts)%in%sigs]=0 #set to 0 values those are not detected
acts[names(acts)%in%sigs]=1 #set to 1 values those are detected


## 2.2. Map signatures to unique segments 
UNIQUE=left_join(dfUnique.filt,CN.EVENTS,by=c("chromosome","start","end","segVal","sample","model")) #100 unique segments
SIG.EXPOSURES=acts #same result after correcting or not by acts
SIG.DEFINITIONS=DEFS

SCORES=as.data.frame(do.call(rbind,lapply(1:nrow(UNIQUE),function(thisEvent){
    probs=UNIQUE[thisEvent,17:53]
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

# Construct the output file including genomic coordinates & matched signature & save
UNIQUE.SCORES=cbind(UNIQUE,SCORES)


# Plot 
UNIQUE.SCORES$map.sig=factor(UNIQUE.SCORES$map.sig,levels = row.names(DEFS))
dfUniqueAssign = UNIQUE.SCORES %>% group_by(model,map.sig,.drop=F) %>% summarise(n=n())
dfUniqueAssign$freq=NA
for(mod in unique(dfUniqueAssign$model)){dfUniqueAssign$freq[dfUniqueAssign$model==mod]=dfUniqueAssign$n[dfUniqueAssign$model==mod]/sum(dfUniqueAssign$n[dfUniqueAssign$model==mod])}

pdf(file=paste0(PLOTS_DIR, "/BarPlot_BRCA_KO_PARPi_Platins_unique.pdf"), width = 100/25.4, height = 100/25.4)
p = ggplot(dfUniqueAssign,aes(x=map.sig,y=freq))+
    facet_wrap(~model, ncol=1) +
    geom_bar(stat = "identity", position = position_dodge())+
    geom_text(aes(label=n), vjust=0, size = 1.8) +
    labs(title = "Signatures assigned to induced segments",
         x = "CIN signatures", y = "Frequency of segments (%)") +
    theme(axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          legend.position = "right")
print(p)
dev.off()

# Barplot for manuscript
# Remove putative random ongoing CIN detection => >0.05% & >2 events
dat=dfUniqueAssign
idx=which(dat$freq<0.05 | dat$n<3)
dat$freq[idx]=0
for(mod in unique(dat$model)){dat$freq[dat$model==mod]=dat$freq[dat$model==mod]/sum(dat$freq[dat$model==mod])}

for (mod in SAMPS){
    pdf(file=paste0(PLOTS_DIR, "/Fig3_BarPlot_",mod,"_KO_PARPi_Platins_unique.pdf"), width = 30/25.4, height = 20/25.4)
    p = ggplot(dat[dat$model==mod,],aes(x=model,y=freq,fill=map.sig))+
        geom_bar(stat = "identity")+
        ylim(0,1) +
        scale_fill_manual(values = sigs_palette) +
        coord_flip() +
        labs(subtitle = "",
             x = "", y = "") +
        theme(axis.text.x = element_text(size = 5),
              axis.text.y = element_blank(),
              legend.position = "none")
    print(p)
    dev.off()
}


################################################################################
#### 3. PLOTTING BULK & ONGOING FEATURES #######################################

#BRCA1
mFeatProbs=rbind(colMeans(ExC.bulk[1:4,]),colMeans(UNIQUE.SCORES[UNIQUE.SCORES$model=="scBRCA1_drug",17:53]))
dtFeatProbs = data.table(melt(mFeatProbs))

# Prepare plot and produce the figure
dtFeatProbs$Col = gsub('[[:digit:]]+', '', dtFeatProbs$Var2 )
dtFeatProbs$Col = factor(dtFeatProbs$Col, levels = c("segsize", "changepoint", "bpwindow"),
                         labels = c("Segment size", "Change point", "Breakpoints per 5MB side"))
dtFeatProbs$Var1 = factor(as.character(dtFeatProbs$Var1), levels = c(1,2), labels = c("bulk","ongoing"))

for(samp in unique(dtFeatProbs$Var1)){
    p = ggplot(dtFeatProbs[dtFeatProbs$Var1==samp,], aes(x = Var2, y = value, fill = Col)) + geom_col() + facet_wrap(. ~ Var1, ncol = 1) +
        scale_fill_manual(values = c("#ffbe0b", "#fb5607", "#ff006e")) +
        scale_x_discrete(guide = guide_axis(angle = 90)) +
        labs(x = "Feature components", y = "Weights", fill = "Feature") +
        theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "bottom") +
        guides(fill = guide_legend(nrow = 1, byrow = TRUE))
    ggsave(paste0(PLOTS_DIR, "/FeatureProbs_Barplot_BRCA1_",samp,".svg"), p, width = 150, height = 60, units = "mm", bg = "white")
    ggsave(paste0(PLOTS_DIR, "/FeatureProbs_Barplot_BRCA1_",samp,".png"), p, width = 150, height = 60, units = "mm", bg = "white")
}

# for correcting BRCA1+Rev experiment!
brca1.feats=colMeans(UNIQUE.SCORES[UNIQUE.SCORES$model=="scBRCA1_drug",17:53])
saveRDS(brca1.feats, file.path(OUTPUT_DIR,"BRCA1_features_unique.rds"))


#BRCA2
mFeatProbs=rbind(colMeans(ExC.bulk[5:7,]),colMeans(UNIQUE.SCORES[UNIQUE.SCORES$model=="scBRCA2_drug",17:53]))
dtFeatProbs = data.table(melt(mFeatProbs))

# Prepare plot and produce the figure
dtFeatProbs$Col = gsub('[[:digit:]]+', '', dtFeatProbs$Var2 )
dtFeatProbs$Col = factor(dtFeatProbs$Col, levels = c("segsize", "changepoint", "bpwindow"),
                         labels = c("Segment size", "Change point", "Breakpoints per 5MB side"))
dtFeatProbs$Var1 = factor(as.character(dtFeatProbs$Var1), levels = c(1,2), labels = c("bulk","ongoing"))

for(samp in unique(dtFeatProbs$Var1)){
    p = ggplot(dtFeatProbs[dtFeatProbs$Var1==samp,], aes(x = Var2, y = value, fill = Col)) + geom_col() + facet_wrap(. ~ Var1, ncol = 1) +
        scale_fill_manual(values = c("#ffbe0b", "#fb5607", "#ff006e")) +
        scale_x_discrete(guide = guide_axis(angle = 90)) + 
        labs(x = "Feature components", y = "Weights", fill = "Feature") +
        theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "bottom") + 
        guides(fill = guide_legend(nrow = 1, byrow = TRUE)) 
    ggsave(paste0(PLOTS_DIR, "/FeatureProbs_Barplot_BRCA2_",samp,".svg"), p, width = 150, height = 60, units = "mm", bg = "white")
    ggsave(paste0(PLOTS_DIR, "/FeatureProbs_Barplot_BRCA2_",samp,".png"), p, width = 150, height = 60, units = "mm", bg = "white")
}



################################################################################
#### 4. PLOTTING BULK & SINGLE-CELL PROFILES ###################################

## Bulk
for(m in unique(BULK.CN$model)){
    pdf(file.path(PLOTS_DIR, paste0("CNProfile_",m,"_KO_bulk.pdf")), width = 180/25.4, height = 120/25.4)
    print(segTab_plot(BULK.CN[BULK.CN$model==m,]))
    dev.off() 
}

## Heatmap with single-cell copy number profiles
for(m in SAMPS){
    plotSCHeatmap(cntable = SCELLS.CN[SCELLS.CN$model==m,], bin.size = BIN.SIZE, path = file.path(BASE, "scripts"),
                  cluster = TRUE, cell_names = FALSE, format = "pdf", outdir = PLOTS_DIR, outname = paste0(m,"_KO_scWGS")) 
}


###### END OF THE SCRIPT! ######################################################