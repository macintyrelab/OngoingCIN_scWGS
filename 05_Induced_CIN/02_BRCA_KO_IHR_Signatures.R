################################################################################
## This script is to process single-cell sequencing data from BRCA1/2 KOs
## Started: 01/08/24
## Updated: 08/05/26
## 
## We have generated BRCA1/2 KO models to induce IHR-related alterations. These
## models were generated using a normal diploid RPE1 cell line, under a TP53 null 
## background. We aim to characterize a pure signature linked to IHR. 
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
PLOTS_DIR=file.path(INDCIN_DIR,"plots/brca_ko")
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


### LOAD INPUT DATA
# Profiles 
SCELLS.CN=readRDS(paste0(INPUT_DIR, "/3.1.0_BRCA_KO_Segments.rds")) #141 cells
# ECNF 
ECNF=readRDS(paste0(INPUT_DIR, "/3.1.0_BRCA_KO_ECNF.rds"))
# ExC matrix
ExC=readRDS(paste0(INPUT_DIR, "/3.1.1_BRCA_KO_ExC.rds"))

# Feature components
MIXMODELS=readRDS(paste0(BASE,"/data/SC_MIXMODELS.rds"))
# Signature definitions
DEFS=readRDS(file=paste0(NEWSIGS_DIR,"/outputs/2_New_scCINSignatures_Definitions_Corrected.rds"))
# DEFS=DEFS[row.names(DEFS)!="CX7",]

# Null distribution matrix => for unique calling
lDistributions=readRDS(paste0(INPUT_DIR, "/confidence/NullDistribution_Confidence_Segsize.rds"))


################################################################################
###### 0. PREPARE DATA #########################################################

### Prepare events table
CN.EVENTS=SCELLS.CN[SCELLS.CN$segVal!=2,]
CN.EVENTS$segVal[CN.EVENTS$segVal<0]=0
CN.EVENTS$length=CN.EVENTS$end-CN.EVENTS$start
CN.EVENTS=CN.EVENTS[CN.EVENTS$length!=0,]
CN.EVENTS=cbind(CN.EVENTS,ExC)


### Compute number of CNAs per cell per model => For correct substraction of background components distributions
num_cnas=as.data.frame(do.call(rbind,lapply(unique(CN.EVENTS$model),function(thisMod){
    dat=CN.EVENTS[CN.EVENTS$model==thisMod,]
    ncnas=nrow(dat)
    ncells=length(as.character(unique(dat$sample)))
    out=c(thisMod,ncnas/ncells)
    return(out)
})))
colnames(num_cnas)=c("model","ncnas_cell")
num_cnas$ncnas_cell=as.numeric(num_cnas$ncnas_cell)

wt.ncnas=num_cnas$ncnas_cell[num_cnas$model=="hTERT_RPE1_WT_48c"]
tp53.ncnas=mean(num_cnas$ncnas_cell[num_cnas$model%in%c("hTERT_RPE1_TP53_B9","hTERT_RPE1_TP53_B3",
                                                        "hTERT_RPE1_TP53_79c","hTERT_RPE1_TP53_72c")])
brca.ncnas=mean(num_cnas$ncnas_cell[num_cnas$model%in%c("hTERT_RPE1_TP53_BRCA1_C8","hTERT_RPE1_BRCA2_D1",
                                                        "hTERT_RPE1_BRCA2_C11","scBRCA2_DLP22")]) #no include hTERT_RPE1_BRCA1 because low efficiency of the CRISPR screening & not induction

### Classify CNAs as clonal, subclonal & unique
SAMPS=unique(SCELLS.CN$model)
lcntables=lapply(SAMPS,function(thisSamp){
    cn=SCELLS.CN[SCELLS.CN$model==thisSamp,1:5]
    return(cn)
})
names(lcntables)=SAMPS

# Get confidence matrix per segment
confidence.files=dir(paste0(INPUT_DIR,"/confidence/brca_ko"),full=T,pattern="Confidence")
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
saveRDS(nUnique.model,file.path(OUTPUT_DIR,"stats_unique_per_model_BRCA_KO.rds"))


################################################################################
###### 1. CHARACTERIZING SIGNATURES IN WT BACKGROUND ###########################
# We compute the sum-of-posterior probabilities for each feature component
# This will be compared the background 
SAMPS=c("hTERT_RPE1_WT_48c")


### 1.1. Sum posterior probabilities of induced events in the WT model
# Only one unique event => no induction of CNAs in culture => no background alterations due to culture => no need to substract
events=dfUnique.filt[dfUnique.filt$model%in%SAMPS,]
events=left_join(events,CN.EVENTS,by=c("chromosome","start","end","segVal","sample","model"))
nrow(events) #1 ! => not accumulation of CNAs in the WT model

# Plot ExC
exc=events[,17:53]
png(file=paste0(PLOTS_DIR, "/2.1.1_Heatmap_ExC_WT.png"), width = 250/25.4, height = 80/25.4, units = "in", res = 300)
p=Heatmap(exc, column_title = "posterior feature components matrix -- RPE1-WT unique events", cluster_rows = TRUE, cluster_columns = FALSE)
print(p)
dev.off()

# Save WT background signature!
wt.sig=as.numeric(exc)
names(wt.sig)=colnames(exc)


################################################################################
###### 2. CHARACTERIZING SIGNATURES LINKED TO TP53-/- BACKGROUND ###############
# We compute the sum-of-posterior probabilities for each feature component
# This will be compared the background 
SAMPS=c("hTERT_RPE1_TP53_79c","hTERT_RPE1_TP53_72c","hTERT_RPE1_TP53_B3","hTERT_RPE1_TP53_B9") #"hTERT_RPE1_TP53_B9_35c" => it is not used for the background because it's only 35 cell cycles & models have ~80-100 cycles => no alterations induced


### 2.1. Sum posterior probabilities of induced events in the TP53 model

# Get events from the TP53-/- model
events=dfUnique.filt[dfUnique.filt$model%in%SAMPS,]
events=left_join(events,CN.EVENTS,by=c("chromosome","start","end","segVal","sample","model"))
events$sample=as.character(events$sample)
table(events$model)
quantile(table(events$sample))
tp53.ncnas=mean(table(events$sample))
nrow(events) #31 recent events! in 20 cells from 3 models

# Plot ExC
exc=events[,17:53]
png(file=paste0(PLOTS_DIR, "/2.2.1_Heatmap_ExC_TP53.png"), width = 250/25.4, height = 200/25.4, units = "in", res = 300)
p=Heatmap(exc, column_title = "posterior feature components matrix -- RPE1-TP53 unique events", cluster_rows = TRUE, cluster_columns = FALSE)
print(p)
dev.off()

# Get SxC per clone
SxC=do.call(rbind,lapply(SAMPS,function(thisSamp){
    sxc=colMeans(events[events$model==thisSamp,17:53])
    return(sxc)
}))
row.names(SxC)=gsub("hTERT_RPE1_","",SAMPS)

# Plot SxC
png(file=paste0(PLOTS_DIR, "/2.2.1_Heatmap_SxC_TP53_clones.png"), width = 250/25.4, height = 150/25.4, units = "in", res = 300)
p=Heatmap(SxC, column_title = "sum-of-posterior feature components matrix -- RPE1-TP53 unique events", cluster_rows = FALSE, cluster_columns = FALSE)
print(p)
dev.off()

# Save TP53 background signature!
tp53.sig=colSums(SxC)
tp53.sig[1:20]=tp53.sig[1:20]/sum(tp53.sig[1:20])
tp53.sig[21:28]=tp53.sig[21:28]/sum(tp53.sig[21:28])
tp53.sig[29:37]=tp53.sig[29:37]/sum(tp53.sig[29:37])
saveRDS(tp53.sig,file.path(OUTPUT_DIR,"2.2.1_TP53_signature_definition_unique.rds"))

# Ploting definition
dtFeatProbs = data.table(melt(tp53.sig))
dtFeatProbs$Var2 = names(tp53.sig)
dtFeatProbs$Col = gsub('[[:digit:]]+', '', dtFeatProbs$Var2 )
dtFeatProbs$Var2 = factor(dtFeatProbs$Var2,levels = names(tp53.sig))
dtFeatProbs$Col = factor(dtFeatProbs$Col, levels = c("segsize", "changepoint", "bpwindow"),
                         labels = c("Segment size", "Change point", "Breakpoints per 5MB side"))

p = ggplot(dtFeatProbs, aes(x = Var2, y = value, fill = Col)) + geom_col() +
    scale_fill_manual(values = c("#ffbe0b", "#fb5607", "#ff006e")) +
    scale_x_discrete(guide = guide_axis(angle = 90)) + 
    labs(x = "Feature components", y = "Weights", fill = "Feature") +
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "bottom") + 
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) 
ggsave(paste0(PLOTS_DIR, "/2.2.1_BarPlot_TP53_signature_definition_unique.svg"), p, width = 150, height = 60, units = "mm", bg = "white")
ggsave(paste0(PLOTS_DIR, "/2.2.1_BarPlot_TP53_signature_definition_unique.png"), p, width = 150, height = 60, units = "mm", bg = "white")



################################################################################
###### 3. CHARACTERIZING SIGNATURES INDUCED BY BRCA1-/- ########################
# We compute the sum-of-posterior probabilities for each feature component
SAMPS=c("hTERT_RPE1_TP53_BRCA1_C8") #"hTERT_RPE1_BRCA1") #=> from CRISPR screening & there is not induced alterations... efficiency of editing?

### 3.1. Sum posterior probabilities of induced events in the TP53+BRCA1 model

# Get unique events from the model
events=dfUnique.filt[dfUnique.filt$model%in%SAMPS,]
events=left_join(events,CN.EVENTS,by=c("chromosome","start","end","segVal","sample","model"))
events$sample=as.character(events$sample)
table(events$model)
quantile(table(events$sample))
brca.ncnas=mean(table(events$sample))
nrow(events) #16 events! in 3 cells of 1 model

# Plot ExC
exc=events[,17:53]
png(file=paste0(PLOTS_DIR, "/2.3.1_Heatmap_ExC_TP53_BRCA1.png"), width = 250/25.4, height = 200/25.4, units = "in", res = 300)
p=Heatmap(exc, column_title = "Posterior feature components matrix -- RPE1-BRCA1+TP53 unique events", cluster_rows = TRUE, cluster_columns = FALSE)
print(p)
dev.off()


# Get SxC per clone
SxC=do.call(rbind,lapply(SAMPS,function(thisSamp){
    sxc=colMeans(events[events$model==thisSamp,17:53])
    return(sxc)
}))
row.names(SxC)=gsub("TP53_","",gsub("hTERT_RPE1_","",SAMPS))

# Plot SxC
png(file=paste0(PLOTS_DIR, "/2.3.1_Heatmap_SxC_TP53_BRCA1_clones.png"), width = 250/25.4, height = 150/25.4, units = "in", res = 300)
p=Heatmap(SxC, column_title = "sum-of-posterior feature components matrix -- RPE1-BRCA1+TP53 unique events", cluster_rows = FALSE, cluster_columns = FALSE)
print(p)
dev.off()

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
ggsave(paste0(PLOTS_DIR, "/2.3.1_BarPlot_TP53+BRCA1_signature_definition_unique.svg"), p, width = 150, height = 60, units = "mm", bg = "white")
ggsave(paste0(PLOTS_DIR, "/2.3.1_BarPlot_TP53+BRCA1_signature_definition_unique.png"), p, width = 150, height = 60, units = "mm", bg = "white")


# Subtract TP53 background 
brca.ncnas=brca.ncnas/(brca.ncnas+tp53.ncnas)
brca1.sig=(as.numeric(SxC)*brca.ncnas)-(tp53.sig*(1-brca.ncnas))
brca1.sig[brca1.sig<0]=0
brca1.sig[1:20]=brca1.sig[1:20]/sum(brca1.sig[1:20])
brca1.sig[21:28]=brca1.sig[21:28]/sum(brca1.sig[21:28])
brca1.sig[29:37]=brca1.sig[29:37]/sum(brca1.sig[29:37])
# Save BRCA1 signature!
saveRDS(brca1.sig,file.path(OUTPUT_DIR,"2.3.1_BRCA1_signature_definition_unique.rds"))

# Ploting definition
dtFeatProbs = data.table(melt(brca1.sig))
dtFeatProbs$Var2 = names(brca1.sig)
dtFeatProbs$Col = gsub('[[:digit:]]+', '', dtFeatProbs$Var2 )
dtFeatProbs$Var2 = factor(dtFeatProbs$Var2,levels = names(tp53.sig))
dtFeatProbs$Col = factor(dtFeatProbs$Col, levels = c("segsize", "changepoint", "bpwindow"),
                         labels = c("Segment size", "Change point", "Breakpoints per 5MB side"))

p = ggplot(dtFeatProbs, aes(x = Var2, y = value, fill = Col)) + geom_col() +
    scale_fill_manual(values = c("#ffbe0b", "#fb5607", "#ff006e")) +
    scale_x_discrete(guide = guide_axis(angle = 90)) + 
    labs(x = "Feature components", y = "Weights", fill = "Feature") +
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "bottom") + 
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) 
ggsave(paste0(PLOTS_DIR, "/2.3.1_BarPlot_BRCA1_signature_definition_unique.svg"), p, width = 150, height = 60, units = "mm", bg = "white")
ggsave(paste0(PLOTS_DIR, "/2.3.1_BarPlot_BRCA1_signature_definition_unique.png"), p, width = 150, height = 60, units = "mm", bg = "white")


# Similarity between BRCA1 & TP53 background
lsa::cosine(brca1.sig,tp53.sig) #0.8626458


### 3.2. Compare with signature definitions
# Pure signatures are defined in a noCIN scenario, whereas the CX signatures are extracted from CIN samples
# This is particularly important for the bpwindow feature
# this feature is specifically designed to explore the genomic context of CNAs in CIN samples
# this is not the case in this models, which have isolated alterations (no matter the mechanism driving it's accumulation)
# For comparisons, we used 2 strategies => raw comparison & exclude bpw

## No bpw comparison
# Cosine similarities
cosine<-do.call(rbind,lapply(1:nrow(DEFS), function(thisSig){
    cos=lsa::cosine(DEFS[thisSig,1:28], brca1.sig[1:28])}))
row.names(cosine)<-row.names(DEFS)

pdf(file=paste0(PLOTS_DIR, "/2.3.2_Heatmap_CXSigs_vs_BRCA1_noBPW.pdf"), width = 250/25.4, height = 80/25.4)
p=Heatmap(t(cosine), col = pal, column_title = "RPE1-TP53+BRCA1 (-TP53)", cluster_rows = FALSE, cluster_columns = FALSE,
          heatmap_legend_param = list(
              title = "Cosine"
          ))
print(p)
dev.off()


################################################################################
###### 4. CHARACTERIZING SIGNATURES INDUCED BY BRCA2-/- ########################
# We compute the sum-of-posterior probabilities for each feature component
SAMPS=c("hTERT_RPE1_BRCA2_D1","scBRCA2_DLP22") #"hTERT_RPE1_BRCA2_C11" => low alterations & with divergent profiles => efficiency of ko? 

### 4.1. Sum posterior probabilities of induced events in the TP53+BRCA2 model

# Get unique events from the model
events=dfUnique.filt[dfUnique.filt$model%in%SAMPS,]
events=left_join(events,CN.EVENTS,by=c("chromosome","start","end","segVal","sample","model"))
# events=events[events$sample!="brca2-B7",]
events$sample=as.character(events$sample)
table(events$model)
quantile(table(events$sample))
brca.ncnas=mean(table(events$sample))
nrow(events) #137 recent events! in 32 cells from 2 models

# Plot ExC
exc=events[,17:53]
png(file=paste0(PLOTS_DIR, "/2.4.1_Heatmap_ExC_TP53_BRCA2.png"), width = 250/25.4, height = 200/25.4, units = "in", res = 300)
p=Heatmap(exc, column_title = "Posterior feature components matrix -- RPE1-BRCA2+TP53 unique events", cluster_rows = TRUE, cluster_columns = FALSE)
print(p)
dev.off()


# Get SxC per clone
SxC=do.call(rbind,lapply(SAMPS,function(thisSamp){
    sxc=colMeans(events[events$model==thisSamp,17:53])
    return(sxc)
}))
row.names(SxC)=gsub("TP53_","",gsub("hTERT_RPE1_","",SAMPS))

# Ploting 
dtFeatProbs = data.table(melt(SxC))
dtFeatProbs$Col = gsub('[[:digit:]]+', '', dtFeatProbs$Var2 )
dtFeatProbs$Col = factor(dtFeatProbs$Col, levels = c("segsize", "changepoint", "bpwindow"),

                                                  labels = c("Segment size", "Change point", "Breakpoints per 5MB side"))
for(samp in unique(dtFeatProbs$Var1)){
    p = ggplot(dtFeatProbs[dtFeatProbs$Var1==samp,], aes(x = Var2, y = value, fill = Col)) + geom_col() +
        scale_fill_manual(values = c("#ffbe0b", "#fb5607", "#ff006e")) +
        scale_x_discrete(guide = guide_axis(angle = 90)) + 
        labs(x = "Feature components", y = "Weights", fill = "Feature") +
        theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "bottom") + 
        guides(fill = guide_legend(nrow = 1, byrow = TRUE)) 
    ggsave(paste0(PLOTS_DIR, "/2.4.1_BarPlot_TP53+BRCA2_signature_definition_unique_", samp,".svg"), p, width = 150, height = 60, units = "mm", bg = "white")
    ggsave(paste0(PLOTS_DIR, "/2.4.1_BarPlot_TP53+BRCA2_signature_definition_unique_", samp,".png"), p, width = 150, height = 60, units = "mm", bg = "white")
}


# Subtract TP53 background 
brca.ncnas=brca.ncnas/(brca.ncnas+tp53.ncnas)
brca2.sig=t(sapply(1:nrow(SxC),function(x) (as.numeric(SxC[x,])*brca.ncnas)-(tp53.sig*(1-brca.ncnas))))
brca2.sig[brca2.sig<0]=0
brca2.sig[,1:20]=t(apply(brca2.sig[,1:20],1,function(x){x/sum(x)}))
brca2.sig[,21:28]=t(apply(brca2.sig[,21:28],1,function(x){x/sum(x)}))
brca2.sig[,29:37]=t(apply(brca2.sig[,29:37],1,function(x){x/sum(x)}))
# Save BRCA2 signature!
saveRDS(brca2.sig,file.path(OUTPUT_DIR,"2.4.1_BRCA2_signature_definition_unique.rds"))

# Ploting definition
dtFeatProbs = data.table(melt(brca2.sig))
dtFeatProbs$Col = gsub('[[:digit:]]+', '', dtFeatProbs$Var2 )
dtFeatProbs$Col = factor(dtFeatProbs$Col, levels = c("segsize", "changepoint", "bpwindow"),
                         labels = c("Segment size", "Change point", "Breakpoints per 5MB side"))

for(samp in unique(dtFeatProbs$Var1)){
    p = ggplot(dtFeatProbs[dtFeatProbs$Var1==samp,], aes(x = Var2, y = value, fill = Col)) + geom_col() +
        scale_fill_manual(values = c("#ffbe0b", "#fb5607", "#ff006e")) +
        scale_x_discrete(guide = guide_axis(angle = 90)) + 
        labs(x = "Feature components", y = "Weights", fill = "Feature") +
        theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "bottom") + 
        guides(fill = guide_legend(nrow = 1, byrow = TRUE)) 
    ggsave(paste0(PLOTS_DIR, "/2.4.1_BarPlot_BRCA2_signature_definition_unique_", samp,".svg"), p, width = 150, height = 60, units = "mm", bg = "white")
    ggsave(paste0(PLOTS_DIR, "/2.4.1_BarPlot_BRCA2_signature_definition_unique_", samp,".png"), p, width = 150, height = 60, units = "mm", bg = "white")
}

# Similarity between BRCA2 & TP53 background
cosine<-do.call(c,lapply(1:nrow(brca2.sig), function(thisClon){
    cos=lsa::cosine(as.numeric(brca2.sig[thisClon,]),tp53.sig)}))
names(cosine)=row.names(SxC)


### 4.2. Compare with signature definitions
# Pure signatures are defined in a noCIN scenario, whereas the CX signatures are extracted from CIN samples
# This is particularly important for the bpwindow feature
# this feature is specifically designed to explore the genomic context of CNAs in CIN samples
# this is not the case in this models, which have isolated alterations (no matter the mechanism driving it's accumulation)
# For comparisons, we used 2 strategies => raw comparison & exclude bpw


## No bpw comparison
# Cosine similarities
cosine<-do.call(rbind,lapply(1:nrow(DEFS), function(thisSig){
    cos=lsa::cosine(DEFS[thisSig,1:28], t(brca2.sig[,1:28]))}))
row.names(cosine)<-row.names(DEFS)

# Plot similarities
png(file=paste0(PLOTS_DIR, "/2.4.2_Heatmap_CXSigs_vs_BRCA2_noBPW.png"), width = 120/25.4, height = 120/25.4, units = "in", res = 300)
p=Heatmap(cosine, column_title = "RPE1-BRCA2 (-TP53) \n cosine similarities -- noBPW", cluster_rows = FALSE, cluster_columns = FALSE,
          cell_fun = function(j, i, x, y, width, height, fill) {
              if(cosine[i, j] > 0.7)
                  grid.text(sprintf("%.2f", cosine[i, j]), x, y, gp = gpar(fontsize = 6))
          })
print(p)
dev.off()


for(samp in 1:ncol(cosine)){
    pdf(file=paste0(PLOTS_DIR, "/2.3.2_Heatmap_CXSigs_vs_BRCA2_noBPW_clone",samp,".pdf"), width = 250/25.4, height = 80/25.4)
    p=Heatmap(t(cosine[,samp]), col = pal, column_title = paste0("RPE1-TP53+BRCA2 (-TP53) -- clone #",samp), cluster_rows = FALSE, cluster_columns = FALSE,
              heatmap_legend_param = list(
                  title = "Cosine"
              ))
    print(p)
    dev.off() 
}


# The unique events generated are mapped to IHR-related signatures!
# > 0.90 in all cases

############ END SCRIPT ! ######################################################