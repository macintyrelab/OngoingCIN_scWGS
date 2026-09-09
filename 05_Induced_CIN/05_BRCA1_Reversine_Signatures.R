################################################################################
## This script is to dissect single-cell sequencing data from BRCA1 KOs treated
## with reversine
## Started: 05/06/26
## Updated: 08/06/26
## 
## We have treated our BRCA1 KO model with reversine to induce mitotic errors.
## The model was generated using a normal diploid RPE1 cell line, under a TP53 null 
## background. We aim to characterize a pure signature linked to mitotic errors
## as the active process causing new alterations. 
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
PLOTS_DIR=file.path(INDCIN_DIR,"plots/brca_rev")
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


### LOAD INPUT DATA
# Profiles 
SCELLS.CN=readRDS(paste0(INPUT_DIR, "/8.1.0_BRCA1-Rev_Segments.rds"))
# ECNF 
ECNF=readRDS(paste0(INPUT_DIR, "/8.1.0_BRCA1-Rev_ECNF.rds"))
# ExC matrix
ExC=readRDS(paste0(INPUT_DIR, "/8.1.1_BRCA1-Rev_ExC.rds"))

# Feature components
MIXMODELS=readRDS(paste0(NEWSIGS_DIR,"/data/SC_MIXMODELS.rds"))
# Signature definitions
DEFS=readRDS(file=paste0(NEWSIGS_DIR,"/outputs/2_New_scCINSignatures_Definitions_Corrected.rds"))

# Null distribution matrix => for unique calling
lDistributions=readRDS(paste0(INPUT_DIR, "/confidence/NullDistribution_Confidence_Segsize.rds"))

# Background signal of TP53+BRCA1
BACK.SIG=readRDS(paste0(BASE,"/06_Ongoing_CIN/outputs/brca_parpi/BRCA1_features_unique.rds"))


################################################################################
###### 1. OBSERVATION OF SINGLE-CELL COPY NUMBER PROFILES ######################

# We observe that reversine treatment induce WGD in some cells. This means that it has
# been an induction of missegregation errors in our TP53+BRCA1 ko
# We also observed gains/losses of whole/arm chromosomes in the diploid cells

## Heatmap with single-cell copy number profiles
plotSCHeatmap(cntable = SCELLS.CN, bin.size = BIN.SIZE, path = file.path(BASE, "scripts"),
              cluster = TRUE, cell_names = TRUE, format = "pdf", outdir = PLOTS_DIR, outname = "BRCA1-Rev_scWGS") 


################################################################################
###### 1. CHARACTERIZING SIGNATURES INDUCED BY REVERSINE #######################

# Reversine signature --> We perform this analysis in the diploid subclone, 
# where we can substract the TP53+BRCA1 background signal 
# Ongoing CIN --> We perform the analysis in each subclone

##### 1.1. PREPARE DATA ########################################################

### Prepare events table
CN.EVENTS=SCELLS.CN[SCELLS.CN$segVal!=2,]
CN.EVENTS$segVal[CN.EVENTS$segVal<0]=0
CN.EVENTS$length=CN.EVENTS$end-CN.EVENTS$start
CN.EVENTS=CN.EVENTS[CN.EVENTS$length!=0,]
CN.EVENTS=cbind(CN.EVENTS,ExC)


### Classify CNAs as clonal, subclonal & unique --> only in diploid cells
SAMPS="BRCA1-Rev"
lcntables=lapply(SAMPS,function(thisSamp){
    cn=SCELLS.CN[SCELLS.CN$model==thisSamp,1:5]
    return(cn)
})
names(lcntables)=SAMPS

# Get confidence matrix per segment
confidence.files=dir(paste0(INPUT_DIR,"/confidence/brca1_rev"),full=T,pattern="Confidence")
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
dfUnique.filt=dfUnique[dfUnique$significant=="yes",] #37 unique events in 18 diploid cells + 56 unique in 3 tetraploid cells


# Add info of tetraploid cells
rm=c("dlp426_C12","dlp426_B15","dlp426_I10")
dfUnique.filt$model[dfUnique.filt$sample%in%rm]="BRCA1-Rev-tetr"
CN.EVENTS$model[CN.EVENTS$sample%in%rm]="BRCA1-Rev-tetr"

##### 1.2. ONGOING CIN #########################################################
defs=DEFS[row.names(DEFS)!="CX7",]

## 1.2.1. Detect which signatures are present => mapping all segments in all cells
sig.scores=as.data.frame(do.call(rbind,lapply(1:nrow(CN.EVENTS),function(e){
    probs=CN.EVENTS[e,8:44]
    scores=sapply(1:nrow(defs),FUN=function(x) sum(defs[x,]*probs))
    scores=scores/sum(scores)
    names(scores)=row.names(defs)
    max=names(which(scores==max(scores)))
    scores=c(scores,max)
    return(scores)
})))
colnames(sig.scores)=c(row.names(defs),"map.sig")
sig.scores[,1:16]=sapply(1:16,function(x)as.numeric(sig.scores[,x]))

sigs=unique(sig.scores$map.sig) #9 signatures detected
acts=colMeans(sig.scores[1:16])
acts[!names(acts)%in%sigs]=0 #set to 0 values those are not detected
acts[names(acts)%in%sigs]=1 #set to 1 values those are detected


## 1.2.2. Map signatures to unique segments 
UNIQUE=left_join(dfUnique.filt,CN.EVENTS,by=c("chromosome","start","end","segVal","sample","model")) #93 unique segments
SIG.EXPOSURES=acts #same result after correcting or not by acts
SIG.DEFINITIONS=defs

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

# Construct the output file including genomic coordinates & matched signature & save
UNIQUE.SCORES=cbind(UNIQUE,SCORES)


# Plot 
UNIQUE.SCORES$map.sig=factor(UNIQUE.SCORES$map.sig,levels = row.names(DEFS))
dfUniqueAssign = UNIQUE.SCORES %>% group_by(model,map.sig,.drop=F) %>% summarise(n=n())
dfUniqueAssign$freq=NA
for(mod in unique(dfUniqueAssign$model)){dfUniqueAssign$freq[dfUniqueAssign$model==mod]=dfUniqueAssign$n[dfUniqueAssign$model==mod]/sum(dfUniqueAssign$n[dfUniqueAssign$model==mod])}

pdf(file=paste0(PLOTS_DIR, "/BarPlot_BRCA1-Rev_unique.pdf"), width = 100/25.4, height = 70/25.4)
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

for (mod in unique(dfUniqueAssign$model)){
    pdf(file=paste0(PLOTS_DIR, "/SuppFigX_BarPlot_",mod,"_unique.pdf"), width = 30/25.4, height = 20/25.4)
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


##### 1.3. REVERSINE SIGNATURE #################################################

### 1.3.1. Sum posterior probabilities of induced events in the TP53+BRCA1 model
# Get unique events from the model
events=left_join(dfUnique.filt[dfUnique.filt$model=="BRCA1-Rev",],CN.EVENTS[CN.EVENTS$model=="BRCA1-Rev",],by=c("chromosome","start","end","segVal","sample","model"))
events$sample=as.character(events$sample)
table(events$model) #37
quantile(table(events$sample))

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
ggsave(paste0(PLOTS_DIR, "/5.1.3_BarPlot_TP53+BRCA1+Rev_signature_definition_unique.svg"), p, width = 150, height = 60, units = "mm", bg = "white")
ggsave(paste0(PLOTS_DIR, "/5.1.3_BarPlot_TP53+BRCA1+Rev_signature_definition_unique.png"), p, width = 150, height = 60, units = "mm", bg = "white")


# Subtract background 
rev.ncnas=0.5 # same indunction of alterations
rev.sig=(as.numeric(SxC)*rev.ncnas)-(BACK.SIG*(1-rev.ncnas))
rev.sig[rev.sig<0]=0
rev.sig[1:20]=rev.sig[1:20]/sum(rev.sig[1:20])
rev.sig[21:28]=rev.sig[21:28]/sum(rev.sig[21:28])
rev.sig[29:37]=rev.sig[29:37]/sum(rev.sig[29:37])


# Ploting definition
dtFeatProbs = data.table(melt(rev.sig))
dtFeatProbs$Var2 = names(rev.sig)
dtFeatProbs$Col = gsub('[[:digit:]]+', '', dtFeatProbs$Var2 )
dtFeatProbs$Var2 = factor(dtFeatProbs$Var2,levels = names(rev.sig))
dtFeatProbs$Col = factor(dtFeatProbs$Col, levels = c("segsize", "changepoint", "bpwindow"),
                         labels = c("Segment size", "Change point", "Breakpoints per 5MB side"))

p = ggplot(dtFeatProbs, aes(x = Var2, y = value, fill = Col)) + geom_col() +
    scale_fill_manual(values = c("#ffbe0b", "#fb5607", "#ff006e")) +
    scale_x_discrete(guide = guide_axis(angle = 90)) + 
    labs(x = "Feature components", y = "Weights", fill = "Feature") +
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "bottom") + 
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) 
ggsave(paste0(PLOTS_DIR, "/5.1.3_BarPlot_Rev_signature_definition_unique.svg"), p, width = 150, height = 60, units = "mm", bg = "white")
ggsave(paste0(PLOTS_DIR, "/5.1.3_BarPlot_Rev_signature_definition_unique.png"), p, width = 150, height = 60, units = "mm", bg = "white")


### 1.3.2. Compare with signature definitions

# Cosine similarities
cosine<-do.call(rbind,lapply(1:nrow(DEFS), function(thisSig){
    cos=lsa::cosine(DEFS[thisSig,1:37], rev.sig[1:37])}))
row.names(cosine)<-row.names(DEFS)

pdf(file=paste0(PLOTS_DIR, "/5.1.3_Heatmap_CXSigs_vs_BRCA1-Rev.pdf"), width = 250/25.4, height = 80/25.4)
p=Heatmap(t(cosine), col = pal, column_title = "RPE1-TP53+BRCA1+Rev (-TP53&BRCA1)", cluster_rows = FALSE, cluster_columns = FALSE,
          heatmap_legend_param = list(
              title = "Cosine"
          ))
print(p)
dev.off()

############ END SCRIPT ! ######################################################