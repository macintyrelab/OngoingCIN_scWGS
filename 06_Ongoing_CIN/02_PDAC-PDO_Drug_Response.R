################################################################################
## This script is to dissect ongoing CIN in a pancreatic cancer organoids for 
## explaining oxiplatin response
## Started: 06/04/25
## Updated: 19/10/25
##
## I derived copy number profiles from single-cell sequencing data generated in 
## the lab using the DLP protocol. I selected those cells that have >30nrpbpcc, 
## and I excluded cells with noisy and/or overfragmented profiles
##
## Here I explore the ability to predict response to oxiplatin using single-cell
## sequencing data, and therefore ongoing CIN signatures. As control, we treated
## these organoids with MRTX113 (KRAS-G12D inhibitor) because one organoid is G12V/WT
## and the other is G12D. So, one should be resistant and the other sensitive.
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
OUTPUT_DIR=file.path(ONGO_DIR,"outputs/pdac_drug")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
PLOTS_DIR=file.path(ONGO_DIR,"plots/pdac_drug")
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
SCELLS.CN=readRDS(paste0(INPUT_DIR, "/4.1.0_PDAC-PDOs_Segments.rds")) 
BULK.CN=readRDS(dir(paste0(BASE,"/data/pdos_pdac"),full=T,pattern="bulk_absoluteCN"))
# ECNF 
ECNF=readRDS(paste0(INPUT_DIR, "/4.1.0_PDAC-PDOs_ECNF.rds"))
# ExC matrix
ExC=readRDS(paste0(INPUT_DIR, "/4.1.1_PDAC-PDOs_ExC.rds"))

# Feature components
MIXMODELS=readRDS(paste0(BASE,"/data/SC_MIXMODELS.rds"))
MODELS=readRDS(paste0(BASE,"/01_Derive_New_CINSignatures/data/NEW_MIXMODELS.rds"))
# Signature definitions
DEFS=readRDS(file=paste0(NEWSIGS_DIR,"/outputs/2_New_scCINSignatures_Definitions_Corrected.rds"))
DEFS=DEFS[row.names(DEFS)!="CX7",]

# Null distribution matrix => for unique calling
lDistributions=readRDS(paste0(INPUT_DIR, "/confidence/NullDistribution_Confidence_Segsize.rds"))

# Drug response
DRUG.RESP=fread(paste0(BASE,"/data/pdos_pdac/drug_response_pdos_pdac.txt"))


################################################################################
#### 0. PREPARE DATA ###########################################################

### 1. PROCESS BULK SEQUENCING DATA 
# Already processed & extracted segTables
BULK.CN$model="MUR54"
BULK.CN$model[BULK.CN$sample=="GM512.markdup"]="MUR67"
BULK.CN$model[BULK.CN$sample=="GM510.markdup"]="MUR59"
BULK.CN$model[BULK.CN$sample=="GM514.markdup"]="MUR72"


### 2. PROCESS SINGLE-CELL DATA

## Prepare events table
CN.EVENTS=SCELLS.CN[SCELLS.CN$segVal!=2,]
CN.EVENTS$segVal[CN.EVENTS$segVal<0]=0
CN.EVENTS$length=CN.EVENTS$end-CN.EVENTS$start
CN.EVENTS=CN.EVENTS[CN.EVENTS$length!=0,]
CN.EVENTS=cbind(CN.EVENTS,ExC)

# select the 2 organoids included in the drug screening
CN.EVENTS=CN.EVENTS[CN.EVENTS$model%in%c("MUR54","MUR67","MUR59", "MUR72"),]
CN.EVENTS$sample=factor(CN.EVENTS$sample,levels=unique(as.character(CN.EVENTS$sample)))

## Classify CNAs as clonal, subclonal & unique
SAMPS=unique(CN.EVENTS$model)
lcntables=lapply(SAMPS,function(thisSamp){
  cn=SCELLS.CN[SCELLS.CN$model==thisSamp,1:5]
  return(cn)
})
names(lcntables)=SAMPS

# Get confidence matrix per segment
confidence.files=dir(paste0(INPUT_DIR,"/confidence/pdac_pdos"),full=T,pattern="Confidence")
confidence.files=confidence.files[c(1:4)]
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
quantile(table(as.character(dfUnique.filt$sample[dfUnique.filt$model==SAMPS[1]]))) #111 unique events (median=2) in 19 cells
quantile(table(as.character(dfUnique.filt$sample[dfUnique.filt$model==SAMPS[2]]))) #24 unique events (median=3) in 17 cells
quantile(table(as.character(dfUnique.filt$sample[dfUnique.filt$model==SAMPS[3]]))) #24 unique events (median=2) in 17 cells
quantile(table(as.character(dfUnique.filt$sample[dfUnique.filt$model==SAMPS[4]]))) #42 unique events (median=2) in 29 cells

# Stats
nUnique.model=dfUnique.filt[dfUnique.filt$significant=="yes",]%>%group_by(model)%>%summarise(n=n())
saveRDS(nUnique.model,file.path(OUTPUT_DIR,"stats_unique_per_model_PDAC-PDO_drug.rds"))



################################################################################
#### 1. PREDICT PLATINE RESISTANCE IN BULK DATA ################################
set.seed(7777)

# We quantify signatures using YAPSA when the sample has at least >15-20 CNAS
# MUR67 has not CIN and therefore may be resistant
table(BULK.CN$model[BULK.CN$segVal!=2])

## 1.1. Extract copy number features
lProfiles=split(BULK.CN, BULK.CN$sample)
ECNF=extractThreeCopynumberFeatures.new(lProfiles, cores=CORES, prePath=paste0(BASE,"/01_Derive_New_CINSignatures/data"), rmNorm=RMNORM, size=SIZE)

## Get SxC matrix
FEATS=names(MODELS)
lMats=lapply(FEATS, function(thisFeature) {
  print(thisFeature)
  thisEcnf=ECNF[[thisFeature]]
  thisModel=MODELS[[thisFeature]]
  
  dat=as.numeric(thisEcnf[,2])
  # We want a posterior, hence likelihood (density) times prior (weight)
  print("Calculating probabilities...")
  if(thisFeature=="bpwindow"){
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
  # Make sure that focal amplifications with huge number of copies have weight in the chp=10 
  # Make sure that windows with a lot of segments have weight in the bpw=9
  id = which(rowSums(is.na(postDatScaled)) == ncol(postDatScaled))
  if (thisFeature == "changepoint"){
    if(length(id)!=0){postDatScaled[id,]=matrix(rep(c(rep(0,9),1), length(id)), ncol = 10, nrow =length(id), byrow = T)}
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
DEFS.BULK=readRDS(file=paste0(BASE, "/01_Derive_New_CINSignatures/outputs/1_NewSignatures_Definitions.rds"))
EXPS.BULK=YAPSA::LCD(t(SxC.psbulk),YAPSA:::normalize_df_per_dim(t(DEFS.BULK),2))
EXPS.BULK=as.matrix(t(apply(EXPS.BULK,2,function(x){x/sum(x)}))) #normalised per sample
row.names(EXPS.BULK)=SAMPS

threshs=fread(file.path(BASE, "/01_Derive_New_CINSignatures/MCsimulations/MC3_Boxplots_per_sig_newCNSigs_fullTCGA_1000sims_10pGaussian_10pSamplePoisson.txt"))
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


# Apply biomarker (Thompson et al. Nat Genetics 2025)
denovo_dir="C:/Users/bhernando/Desktop/CNIO/Projects/scCINSignatures/XX_DeNovo_CINSignatures/outputs"
tcga_exposures <- readRDS(file.path(denovo_dir, 'Activities_NewCINSigs_TCGA_yapsa_k17_THRESH95.rds'))

brca_status <- readRDS(file.path(denovo_dir, 'TCGA_BRCA_MutStatus.rds'))
brca1 = unique(as.character(brca_status$Sample[ grepl("[germline|somatic] BRCA1", 
                                                      brca_status$Status) ]))
brca2 = unique(as.character(brca_status$Sample[ grepl("[germline|somatic] BRCA2", 
                                                      brca_status$Status) ]))
brcaMuts = unique(c(brca1,brca2))
mBRCA = tcga_exposures[row.names(tcga_exposures) %in% brcaMuts, ]
mBRCA = as.matrix(mBRCA[ , c("CX3", "CX2")])
plat_scaling_model = list(mean = attributes(scale(mBRCA))$`scaled:center`, 
                          scale = attributes(scale(mBRCA))$`scaled:scale`)


bulk_plat_prediction <- as.data.frame(EXPS.BULK) %>%
  mutate(sCX2 = (CX2 - plat_scaling_model$mean[['CX2']]) / plat_scaling_model$scale[['CX2']],
         sCX3 = (CX3 - plat_scaling_model$mean[['CX3']]) / plat_scaling_model$scale[['CX3']]) %>%
  mutate(prediction = ifelse(sCX3>sCX2, 'Sensitive', 'Resistant'))


# Conclusion! According to our signature-based biomarker for platinum resistance,
# MUR67 is sensitive (CX3>CX2) while the other 3 organoids are resistant to oxiplatin. 
# This is not what we observe in our drug response data.

## Relative contribution of CX2/CX3
bulk_plat_prediction$model=row.names(bulk_plat_prediction)
dfRELIHR=as.data.frame(melt(bulk_plat_prediction[,c(21,2:3)]))
for(mod in unique(dfRELIHR$model)){dfRELIHR$freq[dfRELIHR$model==mod]=dfRELIHR$value[dfRELIHR$model==mod]/sum(dfRELIHR$value[dfRELIHR$model==mod])}

# Plot CX3/CX2
dfRELIHR$model=factor(dfRELIHR$model,levels = c("MUR72","MUR54","MUR67","MUR59"))
pdf(file.path(PLOTS_DIR, "BarPlot_PDAC-PDOs_CX2_vs_CX3_bulk.pdf"), height = 40/25.4, width = 40/25.4)
p_sigs=ggplot(dfRELIHR, aes(x=model, y=freq, fill=variable)) + 
  geom_bar(stat="identity") +
  ylim(0,1)+
  xlab("") +
  ylab("Relative activity") +
  scale_fill_manual(values=c("#f3e967","#4cbb9d","#0047E6")) + #yellow #green #blue="#0047E6"
  theme(legend.position = "none", legend.title = element_blank()) +
  coord_flip()
print(p_sigs)
dev.off()




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


## 1.1. Detect which signatures are present => mapping all segments in all cells
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
  map.sigs=unique(sig.scores$map.sig[sig.scores$model==thisMod])
  act.sigs=EXPS.BULK[row.names(EXPS.BULK)==thisMod,]
  act.sigs=names(act.sigs[act.sigs>0])
  acts[!names(acts)%in%c(map.sigs,act.sigs)]=0 #set to 0 values those are not detected
  acts[names(acts)%in%c(map.sigs,act.sigs)]=1 #set to 1 values those are detected
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

for (mod in SAMPS){
  pdf(file=paste0(PLOTS_DIR, "/BarPlot_PDAC-PDO_Drug_",mod,"_unique.pdf"), width = 100/25.4, height = 45/25.4)
  p = ggplot(dfUniqueAssign[dfUniqueAssign$model==mod,],aes(x=map.sig,y=freq))+
    facet_wrap(~model, ncol=1) +
    geom_bar(stat = "identity", position = position_dodge())+
    geom_text(aes(label=n), vjust=0, size = 1.8) +
    labs(title = "Signatures assigned to unique segments",
         x = "CIN signatures", y = "Frequency of segments (%)") +
    theme(axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          legend.position = "right")
  print(p)
  dev.off()
}


## 1.3. Relative contribution of ongoing CX2/CX3
dfRELIHR=dfUniqueAssign[dfUniqueAssign$map.sig%in%c("CX2","CX3"),]
for(mod in unique(dfRELIHR$model)){dfRELIHR$freq[dfRELIHR$model==mod]=dfRELIHR$n[dfRELIHR$model==mod]/sum(dfRELIHR$n[dfRELIHR$model==mod])}

# Plot CX3/CX2
dfRELIHR$model=factor(dfRELIHR$model,levels = c("MUR72","MUR54","MUR67","MUR59"))
pdf(file.path(PLOTS_DIR, "BarPlot_PDAC-PDOs_CX2_vs_CX3.pdf"), height = 40/25.4, width = 40/25.4)
p_sigs=ggplot(dfRELIHR, aes(x=model, y=freq, fill=map.sig)) + 
  geom_bar(stat="identity") +
  ylim(0,1)+
  xlab("") +
  ylab("Relative activity") +
  scale_fill_manual(values=c("#f3e967","#4cbb9d","#0047E6")) + #yellow #green #blue="#0047E6"
  theme(legend.position = "none", legend.title = element_blank()) +
  coord_flip()
print(p_sigs)
dev.off()


## 1.4. Drug response differences
# Plot ED50
DRUG.RESP$model=factor(DRUG.RESP$model,levels = c("MUR72","MUR54","MUR67","MUR59"))
p_plat <- ggplot(DRUG.RESP[DRUG.RESP$drug=="Oxaliplatin",], aes(x=ED50, y=model)) +
  geom_point() +
  xlab("ED50 (uM)") +
  theme_bw() +
  scale_x_continuous(trans='log10') +
  theme(legend.position = "none",
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(size = 0.5), 
        axis.ticks = element_line(size = 0.5),
        axis.ticks.length = unit(.1, "cm"))


# Plot ED50
p_KRASi <- ggplot(DRUG.RESP[DRUG.RESP$drug!="Oxaliplatin",], aes(x=ED50, y=model)) +
  geom_point() +
  xlab("ED50 (uM)") +
  theme_bw() +
  scale_x_continuous(trans='log10') +
  theme(legend.position = "none",
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(size = 0.5), 
        axis.ticks = element_line(size = 0.5),
        axis.ticks.length = unit(.1, "cm"))

p=ggarrange(p_plat, p_KRASi, nrow=1, ncol=2, widths=c(0.50, 0.50))

pdf(file.path(PLOTS_DIR, "DrugResponse_ongoing_ed50_PDAC-PDOs.pdf"), height = 40/25.4, width = 60/25.4)
print(p)
dev.off()




################################################################################
#### 3. PLOTTING ONGOING FEATURES ##############################################

mFeatProbs=rbind(colMeans(UNIQUE.SCORES[UNIQUE.SCORES$model==SAMPS[1],16:52]),
                 colMeans(UNIQUE.SCORES[UNIQUE.SCORES$model==SAMPS[2],16:52]),
                 colMeans(UNIQUE.SCORES[UNIQUE.SCORES$model==SAMPS[3],16:52]),
                 colMeans(UNIQUE.SCORES[UNIQUE.SCORES$model==SAMPS[4],16:52]))
dtFeatProbs = data.table(melt(mFeatProbs))

# Prepare plot and produce the figure
dtFeatProbs$Col = gsub('[[:digit:]]+', '', dtFeatProbs$Var2 )
dtFeatProbs$Col = factor(dtFeatProbs$Col, levels = c("segsize", "changepoint", "bpwindow"),
                         labels = c("Segment size", "Change point", "Breakpoints per 5MB side"))
dtFeatProbs$Var1 = factor(as.character(dtFeatProbs$Var1), levels = c(1,2,3,4), labels = SAMPS)

ymax=max(dtFeatProbs$value)
for(samp in unique(dtFeatProbs$Var1)){
  p = ggplot(dtFeatProbs[dtFeatProbs$Var1==samp,], aes(x = Var2, y = value, fill = Col)) + geom_col() + facet_wrap(. ~ Var1, ncol = 1) +
    scale_fill_manual(values = c("#ffbe0b", "#fb5607", "#ff006e")) +
    scale_x_discrete(guide = guide_axis(angle = 90)) + 
    ylim(0,ymax+0.01)+
    labs(x = "Feature components", y = "Weights", fill = "Feature") +
    theme(legend.position = "bottom") + 
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) 
  ggsave(paste0(PLOTS_DIR, "/FeatureProbs_Barplot_",samp,".svg"), p, width = 150, height = 60, units = "mm", bg = "white")
  ggsave(paste0(PLOTS_DIR, "/FeatureProbs_Barplot_",samp,".png"), p, width = 150, height = 60, units = "mm", bg = "white")
}




################################################################################
#### 4. PLOTTING BULK & SINGLE-CELL PROFILES ###################################

## Bulk
pdf(file.path(PLOTS_DIR, "CNProfile_PDAC-PDOs_drug_bulk.pdf"), width = 180/25.4, height = 120/25.4)
for(samp in unique(BULK.CN$model)){
    dat = BULK.CN[BULK.CN$model==samp,]
    print(segTab_plot(dat))
}
dev.off()


## Heatmap with single-cell copy number profiles
for(samp in unique(SCELLS.CN$model)){
  print(samp)
  cn=SCELLS.CN[SCELLS.CN$model==samp,1:5]
  plotSCHeatmap(cntable = cn, bin.size = BIN.SIZE, path = file.path(BASE, "scripts"),
                cluster = TRUE, cell_names = FALSE, format = "pdf", outdir = PLOTS_DIR, outname = paste0("PDAC-PDOs_",samp))
}


###### END OF THE SCRIPT! ######################################################