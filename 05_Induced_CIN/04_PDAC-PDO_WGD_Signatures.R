################################################################################
## This script is to dissect ongoing CIN in a pancreatic cancer organoid that
## has two different clones (ongoing WGD)
## Started: 07/08/25
## Updated: 21/10/25
##
## I derived copy number profiles from single-cell sequencing data generated in 
## the lab using the DLP protocol. I selected those cells that have >30nrpbpcc, 
## and I excluded cells with noisy and/or overfragmented profiles
##
## Here I aim to characterize the WGD signature as the process causing the 
## diversification between both sunclones in the sample
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
INDCIN_DIR=dirname(this.path())
INPUT_DIR=file.path(SCELL_DIR,"intermediate")
# Create directories for outputs
OUTPUT_DIR=file.path(INDCIN_DIR,"outputs")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
PLOTS_DIR=file.path(INDCIN_DIR,"plots/pdac_wgd")
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)


### PLOTTING THEME
theme_set(theme_tufte(base_size = 6, base_family = "ArialMT"))
theme_update(text = element_text(size = 6),
             axis.text = element_text(size = 6),
             axis.line = element_line(size = 0.5), 
             axis.ticks = element_line(size = 0.5),
             axis.ticks.length = unit(.1, "cm"))

set.seed(433)
pal = circlize::colorRamp2(c(0, 0.7, 0.85, 1), c("white", "#FFE5E0", "#FFBEB2", "#E44244"))(seq(0,1,by=0.01))

sigs_palette=c("#4cbcd1", "#db9540","#b16f45", "#4d75b2","#73792f", "#d35138","#ad5fd0", "#6cb543",
               "#ca48a0", "#3a814f","#6861c4", "#b9b04b","#8b99e4", "#ce82bc","#d95070", "#a34e6a", "#5bbf8a")
names(sigs_palette)=paste0("CX", rep(1:17))


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
SIGN.THRESH=0.05 #for defining significant unique segments

RMNORM=TRUE
CORES=1
WIGGLE=0.1
IGNOREDELS=FALSE
SIZE=5000000
UNINFPRIOR=TRUE


### LOAD INPUT DATA
# Profiles 
SCELLS.CN=readRDS(paste0(INPUT_DIR, "/4.1.0_PDAC-PDOs_Segments.rds")) 
BULK.CN=readRDS(dir(paste0(BASE,"/data/pdos_pdac"),full=T,pattern="pseudobulks"))

# ECNF 
ECNF=readRDS(paste0(INPUT_DIR, "/4.1.0_PDAC-PDOs_ECNF.rds"))
# ExC matrix
ExC=readRDS(paste0(INPUT_DIR, "/4.1.1_PDAC-PDOs_ExC.rds"))
# Bins table
SCELLS.BIN=readRDS(paste0(INPUT_DIR, "/4.1.1_PDAC-PDOs_Bins.rds"))[4:5]

# Feature components
MIXMODELS=readRDS(paste0(BASE,"/data/SC_MIXMODELS.rds"))
# Signature definitions
DEFS=readRDS(file=paste0(NEWSIGS_DIR,"/outputs/2_New_scCINSignatures_Definitions_Corrected.rds"))
DEFS=DEFS[row.names(DEFS)!="CX7",] #no subclonal for mapping

# Null distribution matrix => for unique calling
lDistributions=readRDS(paste0(INPUT_DIR, "/confidence/NullDistribution_Confidence_Segsize.rds"))


################################################################################
#### 0. PREPARE DATA ###########################################################

### 1. PROCESS BULK SEQUENCING DATA 
# Already processed & extracted segTables
BULK.CN$model="MUR79"
BULK.CN$model[BULK.CN$sample=="MUR79-pseudo-WGD"]="MUR79-tetr"


### 2. PROCESS SINGLE-CELL DATA

## Prepare events table
CN.EVENTS=SCELLS.CN[SCELLS.CN$segVal!=2,]
CN.EVENTS$segVal[CN.EVENTS$segVal<0]=0
CN.EVENTS$length=CN.EVENTS$end-CN.EVENTS$start
CN.EVENTS=CN.EVENTS[CN.EVENTS$length!=0,]
CN.EVENTS=cbind(CN.EVENTS,ExC)

# Select MUR79
CN.EVENTS=CN.EVENTS[CN.EVENTS$model%in%c("MUR79"),]
CN.EVENTS$sample=factor(CN.EVENTS$sample,levels=unique(as.character(CN.EVENTS$sample)))

# Split copy number profiles per clone => one is diploid & another is tetraploid (with 2 subclones)
tetr.A=c("DLP73-D3","DLP73-F10","DLP73-H5","DLP73-E7","DLP73-B9","DLP73-G8","DLP73-G5","DLP73-E4","DLP73-D1","DLP73-H7")
tetr.B=c("DLP73-C7","DLP73-B7","DLP73-B1","DLP73-E9","DLP73-B3","DLP73-G9","DLP73-B5")
tetr=c(tetr.A,tetr.B)
CN.EVENTS$model[CN.EVENTS$sample%in%tetr]="MUR79-tetr"
SCELLS.CN$model[SCELLS.CN$sample%in%tetr]="MUR79-tetr"


## 2.1. Classify CNAs as clonal, subclonal & unique
# In this case, we compare profiles of the cells in one subclone with the
# pseudobulk of the other subclone to capture the process that it's causing
# the diversification between brances
# So, let's do this

lUnique=lapply(unique(CN.EVENTS$model), function(thisModel){
  print(paste0("Model: ",thisModel))
  refprofile=BULK.CN[BULK.CN$model!=thisModel,]
  refsamp=unique(refprofile$sample)
  scprofiles=SCELLS.CN[SCELLS.CN$model==thisModel,]
  dfUnique=as.data.frame(do.call(rbind,lapply(unique(scprofiles$sample),function(thisCell){
    scprofile=scprofiles[scprofiles$sample==thisCell,]
    cntable=as.data.frame(rbind(refprofile,scprofile))
    
    ## Unify boundaries
    print("Unifying matrix to classify events...")
    # unify segments
    unify=unifySegments.multiple(profiles=cntable[,c(1:5)],window=WINDOW,bin.size=BIN.SIZE)
    # round copy number values
    unify[,4:ncol(unify)]=sapply(4:ncol(unify),function(x)round(unify[,x]))
    # get unique & shared events
    events=rbindlist(classEvents.sc(unify,threshold=THRESH),idcol = T)
    events=as.data.frame(cbind(events[,c(2:5,1,6)]))
    colnames(events)[5]="sample"
    # remove normal segments
    events=events[events$segVal!=2,]
    events=events[events$sample!=refsamp,]
    
    # Get unique segments with high confidence
    unique=events[events$class=="unique",]
    unique.gr <- GRanges(seqnames = unique$chromosome, IRanges(start = unique$start, end = unique$end))
    #genomic locations of segments
    segments.gr = GRanges(seqnames = scprofile$chromosome, IRanges(start = scprofile$start, end = scprofile$end))
    #get overlaps
    overlaps <- findOverlaps(unique.gr, segments.gr)
    # Get segments that map with the unique event in all samples
    seg.unique <- do.call(rbind,lapply(1:nrow(unique), function(thisUnq){
      # get the index of the events of all clones in the same location of unique event 
      idxs=overlaps@to[overlaps@from == thisUnq]
      segs=scprofile[idxs,]
      samp=unique[thisUnq,5]
      usegs=segs[segs$sample==samp,]
      return(usegs)
    }))
    #remove duplicates => unique events can map in the same segment
    seg.unique=seg.unique[!duplicated(seg.unique),]
    seg.unique=seg.unique[order(seg.unique$chromosome,seg.unique$start,seg.unique$end),]
    
    ## Get confidence matrix
    print("Testing unique events...")
    if(nrow(seg.unique)>0){
      # get confidence scores
      confis = as.data.frame(do.call(rbind,lapply(1:nrow(seg.unique), function(thisSeg){
        # get genomic location of segment
        segs=seg.unique[thisSeg,]
        segments.gr = GRanges(seqnames = segs$chromosome, IRanges(start = segs$start, end = segs$end))
        
        # get raw value in specific cell
        bintable=SCELLS.BIN[[thisModel]]
        bin=bintable[bintable$sample%in%segs$sample & bintable$chromosome%in%segs$chromosome,]
        bins.gr=GRanges(seqnames = bin$chromosome, IRanges(start = bin$start, end = bin$end))
        overlaps <- findOverlaps(segments.gr, bins.gr)
        rawVal=mean(bin[overlaps@to,4])
        
        # get values in the refprofile
        ref.gr=GRanges(seqnames = refprofile$chromosome, IRanges(start = refprofile$start, end = refprofile$end))
        overlaps <- findOverlaps(segments.gr, ref.gr)
        refVal=mean(as.data.frame(refprofile[overlaps@to,4])$segVal)
        
        # confidence score => mode of difference with value in other cells
        confi=abs(rawVal-round(refVal))
        min=min(confi)
        confi=as.numeric(names(sort(-table(confi))[1]))
        
        #output
        segs$confidence=confi
        segs$minchange=min
        return(segs)
      })))
      colnames(confis)=c("chromosome","start","end","segVal","sample","model","confidence","minchange")
      
      # Check that unique segments are really unique
      seg.unique=left_join(seg.unique,confis,by=c("chromosome","start","end","segVal","sample","model"))
      
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
  })))
  return(dfUnique)
})
names(lUnique)=unique(CN.EVENTS$model)
dfUnique=rbindlist(lUnique,idcol = F)

# Significant unique segments
dfUnique$t.qvalue=p.adjust(as.numeric(dfUnique$t.pvalue), method = "BH")
dfUnique$significant=ifelse(dfUnique$t.qvalue<SIGN.THRESH,"yes","no")
dfUnique$significant[dfUnique$minchange<0.7]="no" #minimal difference with other cell may be >0.7 to be considered as real
# Retain only high confident unique
dfUnique.filt=dfUnique[dfUnique$significant=="yes",]

# Remove diploid segments & duplicates
dfUnique.filt=dfUnique.filt[dfUnique.filt$segVal!=2,]
dfUnique.filt=dfUnique.filt[!duplicated(dfUnique.filt[,1:4]),]
# Nunique per branch
dfUnique.model=dfUnique.filt%>%group_by(model)%>%summarise(n=n())


################################################################################
#### 1. DECODE CIN IN PSEUDOBULK DATA ##########################################
set.seed(7777)

# We quantify signatures using YAPSA when the sample has at least >10-20 CNAS
table(BULK.CN$model[BULK.CN$segVal!=2])

lProfiles=split(BULK.CN, BULK.CN$sample)
ECNF=extractThreeCopynumberFeatures.sc(lProfiles, cores=CORES, prePath=paste0(BASE,"/01_Derive_New_CINSignatures/data"), rmNorm=RMNORM, size=SIZE)

## Get SxC matrix
FEATS=names(MIXMODELS)
lMats=lapply(FEATS, function(thisFeature) {
  print(thisFeature)
  thisEcnf=ECNF[[thisFeature]]
  thisModel=MIXMODELS[[thisFeature]]
  
  dat=as.numeric(thisEcnf[,2])
  # We want a posterior, hence likelihood (density) times prior (weight)
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
  # Create component matrix per sample
  postDatScaled$Sample = thisEcnf[,1]
  matSxC = aggregate(. ~ Sample, postDatScaled, sum)
  rownames(matSxC) = matSxC$Sample
  matSxC$Sample = NULL
  matSxC = as.matrix(matSxC)
  colnames(matSxC) = paste0(thisFeature, 1:ncol(matSxC))
  
  return(matSxC)
})
names(lMats)=FEATS
SxC.psbulk=do.call(cbind, lMats)


## Quantify sigs
EXPS.BULK=YAPSA::LCD(t(SxC.psbulk),YAPSA:::normalize_df_per_dim(t(DEFS),2))
EXPS.BULK=as.matrix(t(apply(EXPS.BULK,2,function(x){x/sum(x)}))) #normalised per sample
row.names(EXPS.BULK)=SAMPS

threshs=fread(file.path(NEWSIGS_DIR, "/MCsimulations/MC3_Boxplots_per_sig_scCINSigs_fullTCGA_1000sims_10pGaussian_10pSamplePoisson.txt"))
allSigs=colnames(EXPS.BULK)
EXPS.BULK = do.call(cbind,lapply(1:length(allSigs), function(thisSig){
  sig=allSigs[thisSig]
  acts=EXPS.BULK[,thisSig]
  thrs=threshs$Thresh_ZeroGMM_0.05[threshs$Sig==sig]
  acts[acts<thrs]=0
  return(acts)
}))
colnames(EXPS.BULK)=allSigs
EXPS.BULK=t(apply(EXPS.BULK,1,function(x) x/sum(x)))

## Get ExC matrix
FEATS=names(MIXMODELS)
lMats=lapply(FEATS, function(thisFeature) {
  print(thisFeature)
  thisEcnf=ECNF[[thisFeature]]
  thisModel=MIXMODELS[[thisFeature]]
  
  dat=as.numeric(thisEcnf[,2])
  # We want a posterior, hence likelihood (density) times prior (weight)
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
ExC.psbulk=do.call(cbind, lMats)


################################################################################
#### 2. DISSECT ONGOING CIN USING IN EACH BRANCHE ##############################

# Here we are going to map unique events to signatures. 
# Our signature mapping method includes a correction by the signature 
# probabilities found in the sample. 
# We allow to map events to those signatures that are mapped to all cells and events,
# and also to those detected at bulk level.
# The positive point (respect of using only bulk-quantified sigs) is that we are going to allow 
# to map unique events to signature that may be very recent & not detected at bulk level.
# We thus get probabilities for mapping correction by taking into account signature scores 
# in all events of all cells
SAMPS=c("MUR79","MUR79-tetr")
row.names(EXPS.BULK)=SAMPS


## 1.2. Map signatures to unique segments 
UNIQUE=left_join(dfUnique.filt,CN.EVENTS,by=c("chromosome","start","end","segVal","sample","model")) #456 unique events
SIG.EXPOSURES=acts
SIG.DEFINITIONS=DEFS

SCORES=as.data.frame(do.call(rbind,lapply(1:nrow(UNIQUE),function(thisEvent){
  probs=UNIQUE[thisEvent,16:52]
  acts=EXPS.BULK[row.names(EXPS.BULK)==UNIQUE$model[thisEvent],]
  acts[acts>0.1]=1
  scores=sapply(1:nrow(SIG.DEFINITIONS),FUN=function(x) sum(SIG.DEFINITIONS[x,]*probs)*acts[x])
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

# Barplot
ymax=max(dfUniqueAssign$freq)
for (mod in SAMPS){
  pdf(file=paste0(PLOTS_DIR, "/BarPlot_PDAC-PDO_",mod,"_unique.pdf"), width = 100/25.4, height = 45/25.4)
  p = ggplot(dfUniqueAssign[dfUniqueAssign$model==mod,],aes(x=map.sig,y=freq))+
    facet_wrap(~model, ncol=1) +
    geom_bar(stat = "identity", position = position_dodge())+
    ylim(0,ymax+0.01) +
    geom_text(aes(label=n), vjust=0, size = 1.8) +
    labs(title = "Signatures assigned to unique segments",
         x = "CIN signatures", y = "Frequency of segments (%)") +
    theme(axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          legend.position = "right")
  print(p)
  dev.off()
}


# Barplot for manuscript
for (mod in SAMPS){
  # Remove putative random ongoing CIN detection => >0.05% & >2 events
  dat=dfUniqueAssign[dfUniqueAssign$model==mod,]
  idx=which(dat$freq<0.05 | dat$n<3)
  dat$freq[idx]=0
  dat$freq=dat$freq/sum(dat$freq)
  
  pdf(file=paste0(PLOTS_DIR, "/Fig2_BarPlot_PDAC-PDO_",mod,"_unique.pdf"), width = 30/25.4, height = 20/25.4)
  p = ggplot(dat,aes(x=model,y=freq,fill=map.sig))+
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
#### 3. PLOTTING ONGOING FEATURES IN EACH BRANCH ################################

### Diploid clone
# We compute the sum-of-posterior probabilities for each feature component
SAMPS="MUR79"
  
# Get unique events from the model
events=UNIQUE[UNIQUE$model==SAMPS,]

# Get SxC per clone
SxC=do.call(rbind,lapply(SAMPS,function(thisSamp){
  sxc=colMeans(events[events$model==thisSamp,16:52])
  return(sxc)
}))
row.names(SxC)=SAMPS

# Ploting 
dtFeatProbs = data.table(melt(SxC))
dtFeatProbs$Col = gsub('[[:digit:]]+', '', dtFeatProbs$Var2 )
dtFeatProbs$Col = factor(dtFeatProbs$Col, levels = c("segsize", "changepoint", "bpwindow"),
                         labels = c("Segment size", "Change point", "Breakpoints per 5MB side"))

p = ggplot(dtFeatProbs, aes(x = Var2, y = value, fill = Col)) + geom_col() +
  scale_fill_manual(values = c("#ffbe0b", "#fb5607", "#ff006e")) +
  scale_x_discrete(guide = guide_axis(angle = 90)) + 
  labs(x = "Feature components", y = "Weights", fill = "Feature") +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "bottom") + 
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) 
ggsave(paste0(PLOTS_DIR, "/4.3.1_BarPlot_MUR79dip_signature_definition_unique.svg"), p, width = 150, height = 60, units = "mm", bg = "white")
ggsave(paste0(PLOTS_DIR, "/4.3.1_BarPlot_MUR79dip_signature_definition_unique.png"), p, width = 150, height = 60, units = "mm", bg = "white")



### Tetraplod clone
# We compute the sum-of-posterior probabilities for each feature component
SAMPS="MUR79-tetr"

# Get unique events from the model
events=UNIQUE[UNIQUE$model==SAMPS,]

# Get SxC per clone
SxC=do.call(rbind,lapply(SAMPS,function(thisSamp){
  sxc=colMeans(events[events$model==thisSamp,16:52])
  return(sxc)
}))
row.names(SxC)=SAMPS

# Ploting 
dtFeatProbs = data.table(melt(SxC))
dtFeatProbs$Col = gsub('[[:digit:]]+', '', dtFeatProbs$Var2 )
dtFeatProbs$Col = factor(dtFeatProbs$Col, levels = c("segsize", "changepoint", "bpwindow"),
                         labels = c("Segment size", "Change point", "Breakpoints per 5MB side"))

p = ggplot(dtFeatProbs, aes(x = Var2, y = value, fill = Col)) + geom_col() +
  scale_fill_manual(values = c("#ffbe0b", "#fb5607", "#ff006e")) +
  scale_x_discrete(guide = guide_axis(angle = 90)) + 
  labs(x = "Feature components", y = "Weights", fill = "Feature") +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "bottom") + 
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) 
ggsave(paste0(PLOTS_DIR, "/4.3.1_BarPlot_MUR79tetr_signature_definition_unique.svg"), p, width = 150, height = 60, units = "mm", bg = "white")
ggsave(paste0(PLOTS_DIR, "/4.3.1_BarPlot_MUR79tetr_signature_definition_unique.png"), p, width = 150, height = 60, units = "mm", bg = "white")


###### END OF THE SCRIPT! ######################################################

