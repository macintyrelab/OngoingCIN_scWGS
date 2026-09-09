################################################################################
## This script is to process single-cell sequencing data treated with doxorubicin
## Started: 09/06/25
## Updated: 22/07/25
## 
## We treated RPE1 cells with doxorubicin, which is known to induce micronuclei.
## Cells were cultured for some days to allow the amplification of micronuclei.
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
PLOTS_DIR=file.path(INDCIN_DIR,"plots/doxo")
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
source(file.path(BASE, "scripts/uniqueEvents_functions.R"))

set.seed(433)
pal = circlize::colorRamp2(c(0, 0.7, 0.85, 1), c("white", "#FFE5E0", "#FFBEB2", "#E44244"))(seq(0,1,by=0.01))


### PARAMETERS
RESOLUTION="100kb"
BIN.SIZE=as.numeric(substr(RESOLUTION,1,nchar(RESOLUTION)-2))*1000
WINDOW=as.numeric(BIN.SIZE*2)
THRESH=1.00 #% for defining a shared event as clonal/subclonal
SIGN.THRESH=0.01 #for defining significant unique segments


### LOAD INPUT DATA
# Profiles 
SCELLS.CN=readRDS(paste0(INPUT_DIR, "/5.1.0_Drug_Segments.rds")) #53 cells
# ECNF 
ECNF=readRDS(paste0(INPUT_DIR, "/5.1.0_Drug_ECNF.rds"))
# ExC matrix
ExC=readRDS(paste0(INPUT_DIR, "/5.1.1_Drug_ExC.rds"))

# Feature components
MIXMODELS=readRDS(paste0(BASE,"/data/SC_MIXMODELS.rds"))
# Signature definitions
DEFS=readRDS(file=paste0(NEWSIGS_DIR,"/outputs/2_New_scCINSignatures_Definitions_Corrected.rds"))

# Null distribution matrix => for unique calling
lDistributions=readRDS(paste0(INPUT_DIR, "/confidence/NullDistribution_Confidence_Segsize.rds"))


################################################################################
###### 0. PREPARE DATA #########################################################

# Prepare events table
CN.EVENTS=SCELLS.CN[SCELLS.CN$segVal!=2,]
CN.EVENTS$segVal[CN.EVENTS$segVal<0]=0
CN.EVENTS$length=CN.EVENTS$end-CN.EVENTS$start
CN.EVENTS=CN.EVENTS[CN.EVENTS$length!=0,]
CN.EVENTS=cbind(CN.EVENTS,ExC)

# Compute number of CNAs per cell per model => For correct substraction of background components distributions
num_cnas=as.data.frame(do.call(rbind,lapply(unique(CN.EVENTS$model),function(thisMod){
    dat=CN.EVENTS[CN.EVENTS$model==thisMod,]
    ncnas=nrow(dat)
    ncells=length(as.character(unique(dat$sample)))
    out=c(thisMod,ncnas/ncells)
    return(out)
})))
colnames(num_cnas)=c("model","ncnas_cell")
num_cnas$ncnas_cell=as.numeric(num_cnas$ncnas_cell)

ctr.ncnas=mean(num_cnas$ncnas_cell[num_cnas$model=="DMSO"])
dox.ncnas=mean(num_cnas$ncnas_cell[num_cnas$model=="Doxo"])
hu.ncnas=mean(num_cnas$ncnas_cell[num_cnas$model%in%c("HU_100","HU_1.5")])


### Classify CNAs as clonal, subclonal & unique
SAMPS=unique(SCELLS.CN$model)
lcntables=lapply(SAMPS,function(thisSamp){
    cn=SCELLS.CN[SCELLS.CN$model==thisSamp,1:5]
    return(cn)
})
names(lcntables)=SAMPS

# Get confidence matrix per segment
confidence.files=dir(paste0(INPUT_DIR,"/confidence/drug_rpe1"),full=T,pattern="Confidence")
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
# There is already an amplification in the control in the chrm 12
dfUnique.filt=dfUnique.filt[dfUnique.filt$chromosome!=12,]

# Stats
nUnique.model=dfUnique.filt[dfUnique.filt$significant=="yes",]%>%group_by(model)%>%summarise(n=n())
saveRDS(nUnique.model,file.path(OUTPUT_DIR,"stats_unique_per_model_Doxo.rds"))


################################################################################
###### 1. CHARACTERIZING SIGNATURES IN UNTREATED BACKGROUND ####################
# We compute the sum-of-posterior probabilities for each feature component
# This will be compared the background 
SAMPS=c("DMSO")

### 1.1. Sum posterior probabilities of induced events in the TP53 model

# Get events from the DMSO model
events=dfUnique.filt[dfUnique.filt$model%in%SAMPS,]
events=left_join(events,CN.EVENTS,by=c("chromosome","start","end","segVal","sample","model"))
events$sample=as.character(events$sample)
table(events$model)
quantile(table(events$sample))
dmso.ncnas=mean(table(events$sample))
nrow(events) #14 recent events! in 11 cells => these cells are TP53-/- treated with DMSO

# Plot ExC
exc=events[,16:52]
png(file=paste0(PLOTS_DIR, "/3.1.1_Heatmap_ExC_DMSO.png"), width = 250/25.4, height = 80/25.4, units = "in", res = 300)
p=Heatmap(exc, column_title = "posterior feature components matrix -- RPE1 DMSO unique events", cluster_rows = TRUE, cluster_columns = FALSE)
print(p)
dev.off()

# Save DMSO background signature!
dmso.sig=colSums(exc)
dmso.sig[1:20]=dmso.sig[1:20]/sum(dmso.sig[1:20])
dmso.sig[21:28]=dmso.sig[21:28]/sum(dmso.sig[21:28])
dmso.sig[29:37]=dmso.sig[29:37]/sum(dmso.sig[29:37])
names(dmso.sig)=colnames(exc)
saveRDS(dmso.sig,file.path(OUTPUT_DIR,"3.1.1_DMSO_signature_definition.rds"))

png(file=paste0(PLOTS_DIR, "/3.1.1_Heatmap_SxC_DMSO.png"), width = 250/25.4, height = 80/25.4, units = "in", res = 300)
p=Heatmap(t(as.matrix(dmso.sig)), column_title = "sum-of-posterior feature components matrix -- RPE1 DMSO unique events", cluster_rows = FALSE, cluster_columns = FALSE)
print(p)
dev.off()


# Ploting definition
dtFeatProbs = data.table(melt(dmso.sig))
dtFeatProbs$Var2 = names(dmso.sig)
dtFeatProbs$Col = gsub('[[:digit:]]+', '', dtFeatProbs$Var2 )
dtFeatProbs$Var2 = factor(dtFeatProbs$Var2,levels = names(dmso.sig))
dtFeatProbs$Col = factor(dtFeatProbs$Col, levels = c("segsize", "changepoint", "bpwindow"),
                         labels = c("Segment size", "Change point", "Breakpoints per 5MB side"))

p = ggplot(dtFeatProbs, aes(x = Var2, y = value, fill = Col)) + geom_col() +
    scale_fill_manual(values = c("#ffbe0b", "#fb5607", "#ff006e")) +
    scale_x_discrete(guide = guide_axis(angle = 90)) + 
    labs(x = "Feature components", y = "Weights", fill = "Feature") +
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "bottom") + 
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) 
ggsave(paste0(PLOTS_DIR, "/3.1.1_BarPlot_DMSO_signature_definition_unique.svg"), p, width = 150, height = 60, units = "mm", bg = "white")
ggsave(paste0(PLOTS_DIR, "/3.1.1_BarPlot_DMSO_signature_definition_unique.png"), p, width = 150, height = 60, units = "mm", bg = "white")



################################################################################
###### 2. CHARACTERIZING SIGNATURES LINKED TO APHIDICOLIN TREATMENT ############
# We compute the sum-of-posterior probabilities for each feature component
SAMPS=c("Doxo")

### 2.1. Sum posterior probabilities of induced events by aphidicolin

# Get events from the DMSO model
events=dfUnique.filt[dfUnique.filt$model%in%SAMPS,]
events=left_join(events,CN.EVENTS,by=c("chromosome","start","end","segVal","sample","model"))
events$sample=as.character(events$sample)
table(events$model)
quantile(table(events$sample))
doxo.ncnas=mean(table(events$sample[events$segVal>4]))
nrow(events) #17 recent events! in 2 cells but most of them are due to TP53-/- background => allows for accumulation of more CNAs


# Plot ExC
exc=events[,16:52]
png(file=paste0(PLOTS_DIR, "/3.2.1_Heatmap_ExC_Doxo.png"), width = 250/25.4, height = 200/25.4, units = "in", res = 300)
p=Heatmap(exc, column_title = "Posterior feature components matrix -- RPE1 Doxorubicin unique events", cluster_rows = TRUE, cluster_columns = FALSE)
print(p)
dev.off()


# Get SxC per clone
SxC=do.call(rbind,lapply(SAMPS,function(thisSamp){
    sxc=colMeans(events[events$model==thisSamp,16:52])
    return(sxc)
}))
row.names(SxC)=SAMPS

# Plot SxC
png(file=paste0(PLOTS_DIR, "/3.2.1_Heatmap_SxC_DMSO_Doxo.png"), width = 250/25.4, height = 150/25.4, units = "in", res = 300)
p=Heatmap(SxC, column_title = "sum-of-posterior feature components matrix -- RPE1 Doxorubicin unique events", cluster_rows = FALSE, cluster_columns = FALSE)
print(p)
dev.off()

# Plotting 
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
ggsave(paste0(PLOTS_DIR, "/3.2.1_BarPlot_DMSO+Doxo_signature_definition_unique.svg"), p, width = 150, height = 60, units = "mm", bg = "white")
ggsave(paste0(PLOTS_DIR, "/3.1.1_BarPlot_DMSO+Doxo_signature_definition_unique.png"), p, width = 150, height = 60, units = "mm", bg = "white")


# Subtract DMSO background 
doxo.ncnas=doxo.ncnas/(doxo.ncnas+dmso.ncnas)
doxo.sig=(as.numeric(SxC)*doxo.ncnas)-(dmso.sig*(1-doxo.ncnas))
doxo.sig[doxo.sig<0]=0
doxo.sig[1:20]=doxo.sig[1:20]/sum(doxo.sig[1:20])
doxo.sig[21:28]=doxo.sig[21:28]/sum(doxo.sig[21:28])
doxo.sig[29:37]=doxo.sig[29:37]/sum(doxo.sig[29:37])
# Save Doxo signature!
saveRDS(doxo.sig,file.path(OUTPUT_DIR,"3.2.1_Doxo_signature_definition_unique.rds"))

png(file=paste0(PLOTS_DIR, "/3.2.1_Heatmap_SxC_Doxo.png"), width = 250/25.4, height = 80/25.4, units = "in", res = 300)
p=Heatmap(t(as.matrix(doxo.sig)), column_title = "sum-of-posterior feature components matrix -- RPE1 Doxorubicin (-DMSO) unique events", cluster_rows = FALSE, cluster_columns = FALSE)
print(p)
dev.off()


# Ploting definition
dtFeatProbs = data.table(melt(doxo.sig))
dtFeatProbs$Var2 = names(doxo.sig)
dtFeatProbs$Col = gsub('[[:digit:]]+', '', dtFeatProbs$Var2 )
dtFeatProbs$Var2 = factor(dtFeatProbs$Var2,levels = names(doxo.sig))
dtFeatProbs$Col = factor(dtFeatProbs$Col, levels = c("segsize", "changepoint", "bpwindow"),
                         labels = c("Segment size", "Change point", "Breakpoints per 5MB side"))

p = ggplot(dtFeatProbs, aes(x = Var2, y = value, fill = Col)) + geom_col() +
    scale_fill_manual(values = c("#ffbe0b", "#fb5607", "#ff006e")) +
    scale_x_discrete(guide = guide_axis(angle = 90)) + 
    labs(x = "Feature components", y = "Weights", fill = "Feature") +
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "bottom") + 
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) 
ggsave(paste0(PLOTS_DIR, "/3.2.1_BarPlot_Doxo_signature_definition_unique.svg"), p, width = 150, height = 60, units = "mm", bg = "white")
ggsave(paste0(PLOTS_DIR, "/3.1.1_BarPlot_Doxo_signature_definition_unique.png"), p, width = 150, height = 60, units = "mm", bg = "white")


# Similarity between Doxo & DMSO background
lsa::cosine(doxo.sig,dmso.sig) #0.4367666


### 2.2. Compare with signature definitions
# Pure signatures are defined in a noCIN scenario, whereas the CX signatures are extracted from CIN samples
# This is particularly important for the bpwindow feature
# this feature is specifically designed to explore the genomic context of CNAs in CIN samples
# this is not the case in this models, which have isolated alterations (no matter the mechanism driving it's accumulation)
# For comparisons, we used 2 strategies => raw comparison & exclude bpw

## No bpw comparison
# Cosine similarities
cosine<-do.call(rbind,lapply(1:nrow(DEFS), function(thisSig){
    cos=round(lsa::cosine(DEFS[thisSig,1:28], doxo.sig[1:28]),2)}))
row.names(cosine)<-row.names(DEFS)

# Plot similarities
pdf(file=paste0(PLOTS_DIR, "/3.2.2_Heatmap_CXSigs_vs_Doxo_noBPW.pdf"), width = 250/25.4, height = 80/25.4)
p=Heatmap(t(cosine), col = pal, column_title = "RPE1-TP53 Doxorubicin (-DMSO)", cluster_rows = FALSE, cluster_columns = FALSE,
          heatmap_legend_param = list(
              title = "Cosine",
              at = c(0, 0.5, 0.75, 1),  # fixed legend ticks
              labels = c("0", "0.5", "0.75", "1")  # labels for legend
          ))
print(p)
dev.off()


# The unique events generated are mapped to CX8 - an amplification-generating signature !!

############ END SCRIPT ! ######################################################