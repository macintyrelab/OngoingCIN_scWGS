################################################################################
## This script is to dissect ongoing CIN in a TNBC tumours from Minussi et al 2022
## Started: 18/04/25
## Updated: 22/10/25
##
## Single-cell DNA sequencing data from Navin's lab (https://www.nature.com/articles/s41586-021-03357-x)
## They applied ACT to sequence 9,765 cells from 8 TNBC tumours, and used this data
## to reconstruct phylogenetic trees and evolution of these tumours.
##
## Here, I derived copy number profiles of the subclones found in each sample by
## concatenating fasta files of all cells included in each subclone. Then, I call
## unique events found in each subclone and map signatures to these unique events.
## This allows me to know the mutational signature driving the acquisition of new events.
##
## Note: we explore signatures at subclone level because the number of reads per cell
## it's not enough for having a significant number of cells per tumour. 
## We have only the desired resolution of TN1 and TN2, which are used for deriving
## the null distribution for calling unique events
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
OUTPUT_DIR=file.path(ONGO_DIR,"outputs/minussi")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
PLOTS_DIR=file.path(ONGO_DIR,"plots/minussi")
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
SCELLS.CN=readRDS(paste0(INPUT_DIR, "/6.1.0_subMinussi_Powered_Segments.rds")) 
BULK.CN=readRDS(paste0(BASE,"/data/minussi_psbk/minussi_pseudobulk_segment_tables.rds"))
# ECNF 
ECNF=readRDS(paste0(INPUT_DIR, "/6.1.0_subMinussi_ECNF.rds"))
# ExC matrix
ExC=readRDS(paste0(INPUT_DIR, "/6.1.1_subMinussi_ExC.rds"))

# Feature components
MIXMODELS=readRDS(paste0(BASE,"/data/SC_MIXMODELS.rds"))
MODELS=readRDS(paste0(BASE,"/01_Derive_New_CINSignatures/data/NEW_MIXMODELS.rds"))
# Signature definitions
DEFS=readRDS(file=paste0(NEWSIGS_DIR,"/outputs/2_New_scCINSignatures_Definitions_Corrected.rds"))
DEFS=DEFS[row.names(DEFS)!="CX7",]

# Null distribution matrix => for unique calling
lDistributions=readRDS(paste0(INPUT_DIR, "/confidence/NullDistribution_Confidence_Segsize.rds"))


################################################################################
#### 0. PREPARE DATA ###########################################################

### 1. PROCESS BULK SEQUENCING DATA 
# Already processed & extracted segTables
BULK.CN$model=gsub("_pseudobulk_R1.markdup","",BULK.CN$sample)


### 2. PROCESS SINGLE-CELL DATA

## Prepare events table
CN.EVENTS=SCELLS.CN[SCELLS.CN$segVal!=2,]
CN.EVENTS$segVal[CN.EVENTS$segVal<0]=0
CN.EVENTS$length=CN.EVENTS$end-CN.EVENTS$start
CN.EVENTS=CN.EVENTS[CN.EVENTS$length!=0,]
CN.EVENTS=cbind(CN.EVENTS,ExC)

## Classify CNAs as clonal, subclonal & unique
SAMPS=unique(CN.EVENTS$model)
lcntables=lapply(SAMPS,function(thisSamp){
  cn=SCELLS.CN[SCELLS.CN$model==thisSamp,1:5]
  return(cn)
})
names(lcntables)=SAMPS

# Get confidence matrix per segment
confidence.files=dir(paste0(INPUT_DIR,"/confidence/minussi_sub"),full=T,pattern="Confidence")
confidence.files=confidence.files[c(3,6:12)] #remove branches for TN1 & TN2
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
# Stats of unique
quantile(table(as.character(dfUnique.filt$sample[dfUnique.filt$model==SAMPS[1]]))) #52 unique events (median=3) in 17 clones
quantile(table(as.character(dfUnique.filt$sample[dfUnique.filt$model==SAMPS[2]]))) #36 unique events (median=2) in 14 clones
quantile(table(as.character(dfUnique.filt$sample[dfUnique.filt$model==SAMPS[3]]))) #27 unique events (median=3.5) in 9 clones
quantile(table(as.character(dfUnique.filt$sample[dfUnique.filt$model==SAMPS[4]]))) #128 unique events (median=3.5) in 23 clones
quantile(table(as.character(dfUnique.filt$sample[dfUnique.filt$model==SAMPS[5]]))) #35 unique events (median=7) in 5 clones
quantile(table(as.character(dfUnique.filt$sample[dfUnique.filt$model==SAMPS[6]]))) #13 unique events (median=1) in 17 clones
quantile(table(as.character(dfUnique.filt$sample[dfUnique.filt$model==SAMPS[7]]))) #32 unique events (median=1.5) in 22 clones
quantile(table(as.character(dfUnique.filt$sample[dfUnique.filt$model==SAMPS[8]]))) #20 unique events (median=2) in 15 clones

# Stats
nUnique.model=dfUnique.filt[dfUnique.filt$significant=="yes",]%>%group_by(model)%>%summarise(n=n())
saveRDS(nUnique.model,file.path(OUTPUT_DIR,"stats_unique_per_model_Minussi.rds"))


################################################################################
#### 1. DECODE CIN IN PSEUDOBULK DATA ##########################################
set.seed(7777)

# We quantify signatures using YAPSA when the sample has at least >10-20 CNAS
table(BULK.CN$model[BULK.CN$segVal!=2])

lProfiles=split(BULK.CN, BULK.CN$sample)
ECNF=extractThreeCopynumberFeatures.sc(lProfiles, cores=CORES, prePath=paste0(BASE,"/01_Derive_New_CINSignatures/data"), rmNorm=RMNORM, size=SIZE, res=BIN.SIZE)

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


# Quantify sigs
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
#### 2. DISSECT ONGOING CIN USING SINGLE-CELL DATA #############################

# Here we are going to map unique events (recently acquired due to an ongoing process) 
# to signatures. Our signature mapping method includes a correction by the signature 
# probabilities found in the sample. 
# We allow to map events to those signatures that are mapped to all cells and events,
# and also to those detected at bulk level.
# The positive point (respect of using only bulk-quantified sigs) is that we are going to allow 
# to map unique events to signature that may be very recent & not detected at bulk level.
# We thus get probabilities for mapping correction by taking into account signature scores 
# in all events of all cells


## 2.1. Detect which signatures are present => mapping all segments in all cells
sig.scores=as.data.frame(do.call(rbind,lapply(1:nrow(CN.EVENTS),function(e){
  probs=CN.EVENTS[e,8:44]
  model=CN.EVENTS$model[e]
  scores=sapply(1:nrow(DEFS),FUN=function(x) sum(DEFS[x,]*probs))
  scores=scores/sum(scores)
  names(scores)=row.names(DEFS)
  max=names(which(scores==max(scores)))
  scores=c(scores,max,model)
  return(scores)
})))
colnames(sig.scores)=c(row.names(DEFS),"map.sig","model")
sig.scores[,1:16]=sapply(1:16,function(x)as.numeric(sig.scores[,x]))

acts=as.matrix(do.call(rbind,lapply(SAMPS,function(thisMod){
  acts=colMeans(sig.scores[sig.scores$model==thisMod,1:16])
  map.sigs=unique(sig.scores$map.sig[sig.scores$model==thisMod]) #mapped using all events => highest probability in at least one event
  act.sigs=EXPS.BULK[row.names(EXPS.BULK)==thisMod,] #detected at bulk level
  act.sigs=names(act.sigs[act.sigs>0])
  acts[!names(acts)%in%c(map.sigs,act.sigs)]=0 #set to 0 values those are not detected/mapped 
  acts[names(acts)%in%c(map.sigs,act.sigs)]=1 #set to 1 values those are detected/mapped 
  return(acts)
})))
row.names(acts)=SAMPS


## 1.2. Map signatures to unique segments 
UNIQUE=left_join(dfUnique.filt,CN.EVENTS,by=c("chromosome","start","end","segVal","sample","model")) #153 unique segments
SIG.EXPOSURES=acts
SIG.DEFINITIONS=DEFS

SCORES=as.data.frame(do.call(rbind,lapply(1:nrow(UNIQUE),function(thisEvent){
  probs=UNIQUE[thisEvent,16:52]
  if(is.null(SIG.EXPOSURES)){scores=sapply(1:nrow(SIG.DEFINITIONS),FUN=function(x) sum(SIG.DEFINITIONS[x,]*probs))}
  if(!is.null(SIG.EXPOSURES)){
    acts=SIG.EXPOSURES[row.names(SIG.EXPOSURES)==UNIQUE$model[thisEvent],]
    scores=sapply(1:nrow(SIG.DEFINITIONS),FUN=function(x) sum(SIG.DEFINITIONS[x,]*probs)*acts[x])}
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
  pdf(file=paste0(PLOTS_DIR, "/BarPlot_subMinussi_",mod,"_unique_permodel.pdf"), width = 100/25.4, height = 45/25.4)
  p = ggplot(dfUniqueAssign[dfUniqueAssign$model==mod,],aes(x=map.sig,y=freq))+
    facet_wrap(~model, ncol=1) +
    geom_bar(stat = "identity", position = position_dodge())+
    ylim(0,ymax+0.01) +
    geom_text(aes(label=n), vjust=0, size = 1.8) +
    labs(title = "Signatures assigned to unique segments",
         x = "CIN signatures", y = "Frequency of segments (%)") +
    theme(axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          legend.position = "none")
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
  
  pdf(file=paste0(PLOTS_DIR, "/Fig4_BarPlot_subMinussi_",mod,"_unique_permodel.pdf"), width = 30/25.4, height = 20/25.4)
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


# Heatmap
set.seed(433)
pal = circlize::colorRamp2(c(0, 0.1, 0.25, 0.3, 0.5, 0.6, 0.7), c("white", "#f5babb", "#eb7577", "#e65355","#e23133","#ce1d20","#ab181a"))(seq(0,0.7,by=0.01))
pal = circlize::colorRamp2(breaks=seq(0,0.7,by=0.01),colors=pal)

for(samp in unique(UNIQUE.SCORES$model)){
  # Remove putative random ongoing CIN detection => >0.05% & >2 events
  dat=dfUniqueAssign[dfUniqueAssign$model==samp,]
  idx=which(dat$freq<0.05 | dat$n<3)
  dat$freq[idx]=0
  acts=dat$freq/sum(dat$freq)
  
  names(acts)=row.names(DEFS)
  pdf(file=paste0(PLOTS_DIR, "/CINSig_Frequency_Heatmap_",samp,"_ongoing.pdf"), width = 250/25.4, height = 80/25.4)
  p=Heatmap(t(acts), col = pal, column_title = samp, cluster_rows = FALSE, cluster_columns = FALSE,
            heatmap_legend_param = list(
              title = "SigFrq"
            ))
  print(p)
  dev.off()
}

# Heatmap with single-clone copy number profiles per sample
for(samp in unique(UNIQUE.SCORES$model)){
  cn=SCELLS.CN[SCELLS.CN$model==samp,1:5]
  plotSCHeatmap(cntable = cn, bin.size = BIN.SIZE, path = file.path(BASE, "scripts"),
                cluster = TRUE, cell_names = FALSE, format = "pdf", outdir = PLOTS_DIR, outname = paste0("subMinussi_",samp))
}

################################################################################
#### 3. PLOTTING BULK & ONGOING FEATURES #######################################

## Pseudobulk
dtFeatProbs = data.table(melt(SxC.psbulk))

# Prepare plot and produce the figure
dtFeatProbs$Col = gsub('[[:digit:]]+', '', dtFeatProbs$Var2 )
dtFeatProbs$Col = factor(dtFeatProbs$Col, levels = c("segsize", "changepoint", "bpwindow"),
                         labels = c("Segment size", "Change point", "Breakpoints per 5MB side"))
dtFeatProbs$Var1 = gsub("_pseudobulk_R1.markdup","",dtFeatProbs$Var1)
dtFeatProbs$Var1 = factor(as.character(dtFeatProbs$Var1), levels = SAMPS)

for(samp in unique(dtFeatProbs$Var1)){
  p = ggplot(dtFeatProbs[dtFeatProbs$Var1==samp,], aes(x = Var2, y = value, fill = Col)) + geom_col() + facet_wrap(. ~ Var1, ncol = 1) +
    scale_fill_manual(values = c("#ffbe0b", "#fb5607", "#ff006e")) +
    scale_x_discrete(guide = guide_axis(angle = 90)) + 
    labs(x = "Feature components", y = "Weights", fill = "Feature") +
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "bottom") +  
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) 
  ggsave(paste0(PLOTS_DIR, "/FeatureProbs_Barplot_",samp,"_psbk.svg"), p, width = 150, height = 60, units = "mm", bg = "white")
  ggsave(paste0(PLOTS_DIR, "/FeatureProbs_Barplot_",samp,"_psbk.png"), p, width = 150, height = 60, units = "mm", bg = "white")
}


## Ongoing
mFeatProbs=rbind(colMeans(UNIQUE.SCORES[UNIQUE.SCORES$model==SAMPS[1],16:52]),
                 colMeans(UNIQUE.SCORES[UNIQUE.SCORES$model==SAMPS[2],16:52]),
                 colMeans(UNIQUE.SCORES[UNIQUE.SCORES$model==SAMPS[3],16:52]),
                 colMeans(UNIQUE.SCORES[UNIQUE.SCORES$model==SAMPS[4],16:52]),
                 colMeans(UNIQUE.SCORES[UNIQUE.SCORES$model==SAMPS[5],16:52]),
                 colMeans(UNIQUE.SCORES[UNIQUE.SCORES$model==SAMPS[6],16:52]),
                 colMeans(UNIQUE.SCORES[UNIQUE.SCORES$model==SAMPS[7],16:52]),
                 colMeans(UNIQUE.SCORES[UNIQUE.SCORES$model==SAMPS[8],16:52]))
dtFeatProbs = data.table(melt(mFeatProbs))

# Prepare plot and produce the figure
dtFeatProbs$Col = gsub('[[:digit:]]+', '', dtFeatProbs$Var2 )
dtFeatProbs$Col = factor(dtFeatProbs$Col, levels = c("segsize", "changepoint", "bpwindow"),
                         labels = c("Segment size", "Change point", "Breakpoints per 5MB side"))

dtFeatProbs$Var1 = factor(as.character(dtFeatProbs$Var1), levels = 1:8, labels = SAMPS)

ymax=max(dtFeatProbs$value)
for(samp in unique(dtFeatProbs$Var1)){
  p = ggplot(dtFeatProbs[dtFeatProbs$Var1==samp,], aes(x = Var2, y = value, fill = Col)) + geom_col() + facet_wrap(. ~ Var1, ncol = 1) +
    scale_fill_manual(values = c("#ffbe0b", "#fb5607", "#ff006e")) +
    scale_x_discrete(guide = guide_axis(angle = 90)) + 
    labs(x = "Feature components", y = "Weights", fill = "Feature") +
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "bottom") +  
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) 
  ggsave(paste0(PLOTS_DIR, "/FeatureProbs_Barplot_",samp,".svg"), p, width = 150, height = 60, units = "mm", bg = "white")
  ggsave(paste0(PLOTS_DIR, "/FeatureProbs_Barplot_",samp,".png"), p, width = 150, height = 60, units = "mm", bg = "white")
}



################################################################################
#### 4. PLOTTING BULK & ONGOING SIGNATURES #####################################


# Map pseudobulk events
BULK.EVENTS=BULK.CN[BULK.CN$segVal!=2,]
BULK.EVENTS=cbind(BULK.EVENTS,ExC.psbulk)
SCORES=as.data.frame(do.call(rbind,lapply(1:nrow(BULK.EVENTS),function(thisEvent){
  probs=BULK.EVENTS[thisEvent,7:43]
  acts=EXPS.BULK[row.names(EXPS.BULK)==BULK.EVENTS$model[thisEvent],]
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
BULK.SCORES=cbind(BULK.EVENTS,SCORES)

# Plot 
BULK.SCORES$map.sig=factor(BULK.SCORES$map.sig,levels = row.names(DEFS))
dfBulkAssign = BULK.SCORES %>% group_by(model,map.sig,.drop=F) %>% summarise(n=n())
dfBulkAssign$freq=NA
for(mod in unique(dfBulkAssign$model)){dfBulkAssign$freq[dfBulkAssign$model==mod]=dfBulkAssign$n[dfBulkAssign$model==mod]/sum(dfBulkAssign$n[dfBulkAssign$model==mod])}


# Barplot for manuscript
for (samp in unique(BULK.SCORES$model)){
  dat=dfBulkAssign[dfBulkAssign$model==samp,]
  
  pdf(file=paste0(PLOTS_DIR, "/Fig4_BarPlot_subMinussi_",samp,"_pseudobulk_permodel.pdf"), width = 30/25.4, height = 20/25.4)
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

# Heatmap
set.seed(433)
pal = circlize::colorRamp2(c(0, 0.1, 0.25, 0.3, 0.5, 0.6, 0.7), c("white", "#f5babb", "#eb7577", "#e65355","#e23133","#ce1d20","#ab181a"))(seq(0,0.7,by=0.01))
pal = circlize::colorRamp2(breaks=seq(0,0.7,by=0.01),colors=pal)

for(samp in unique(BULK.SCORES$model)){
  acts=dfBulkAssign$freq[dfBulkAssign$model==samp]
  names(acts)=row.names(DEFS)
  pdf(file=paste0(PLOTS_DIR, "/CINSig_Frequency_Heatmap_",samp,"_pseudobulk.pdf"), width = 250/25.4, height = 80/25.4)
  p=Heatmap(t(acts), col = pal, column_title = samp, cluster_rows = FALSE, cluster_columns = FALSE,
            heatmap_legend_param = list(
              title = "SigFrq"
            ))
  print(p)
  dev.off()
}


####### END OF THE SCRIPT! #####################################################

