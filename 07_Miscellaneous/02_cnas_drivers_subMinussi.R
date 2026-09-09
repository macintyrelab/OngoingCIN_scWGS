################################################################################
## This script is to explore amps/dels in subclones from Minussi data
## Started: 24/10/25
## Updated: 24/10/25
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
MISC_DIR=dirname(this.path())
SCELL_DIR=file.path(BASE,"04_Process_scWGS")
INPUT_DIR=file.path(BASE,"data")
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
source(file.path(BASE, "scripts/main_functions.R"))


### LOAD INPUT DATA
# Profiles 
CN.TABLE=readRDS(paste0(SCELL_DIR, "/intermediate/6.1.0_subMinussi_Powered_Segments.rds")) 

# Genes of interest
CANCER.GENES=readRDS(paste0(INPUT_DIR,"/Reactome_and_DriverGenes_ConsensusAnnotations.rds"))

################################################################################

## Extract copy number values of drivers per subclone
samps.gr = GRanges(seqnames = CN.TABLE$chromosome, IRanges(start = CN.TABLE$start, end = CN.TABLE$end))

# Gene region
genes.gr <- GRanges(seqnames = CANCER.GENES$chromosome, IRanges(start = CANCER.GENES$start))
genes.gr$gene_name <- CANCER.GENES$name

# Find overlaps
overlaps <- findOverlaps(genes.gr, samps.gr)

# Get copy number value of gene per sample
lCNs <- lapply(1:nrow(CANCER.GENES), function(thisGene){
    # get the index of the samples with alterations in the different genes (interesting for later)
    idxs=overlaps@to[overlaps@from == thisGene]
    cn=as.data.frame(CN.TABLE[idxs,4:5])
    return(cn)
})
names(lCNs) <- genes.gr$gene_name
CNAS.GOI = rbindlist(lCNs, idcol=T)
colnames(CNAS.GOI) = c("gene", "cn","sample")


## Get ploidies of each subclone
dfPloidy=as.data.frame(do.call(rbind,lapply(unique(CN.TABLE$sample), function(s){
    segTab<-CN.TABLE[CN.TABLE$sample==s,]
    segLen<-(as.numeric(segTab$end)-as.numeric(segTab$start))
    ploidy<-sum((segLen/sum(segLen))*as.numeric(segTab$segVal))
    return(c(unique(segTab$model),s,ploidy))
})))
colnames(dfPloidy)=c("model","sample","ploidy")
dfPloidy$ploidy=as.numeric(dfPloidy$ploidy)


## Classify as amps or dels based on consensus value
CNAS.GOI=left_join(CNAS.GOI,dfPloidy,by="sample")

dfDiffs=c()
for(m in unique(CNAS.GOI$model)){
    df=CNAS.GOI[CNAS.GOI$model==m,]
    diffs=as.data.frame(do.call(rbind,lapply(unique(df$gene), function(thisGene){
        vals=df$cn[df$gene==thisGene]
        vals=vals-(median(vals))
        samps=df$sample[df$gene==thisGene]
        return(cbind(thisGene,vals,samps))
    })))
    colnames(diffs)=c("gene","diff","sample")
    diffs$model=m
    dfDiffs=rbind(dfDiffs,diffs)
}
dfDiffs$diff=as.numeric(dfDiffs$diff)

# Remove no changes
dfDiffs = dfDiffs[dfDiffs$diff!=0,]
dfDiffs$consequence="DEL"
dfDiffs$consequence[dfDiffs$diff>0]="AMP"

# Gene type
dfDiffs$gene_status="TSG"
dfDiffs$gene_status[dfDiffs$gene%in%CANCER.GENES$name[CANCER.GENES$Status == "Oncogene" ]]="OG"
dfDiffs$gene_status[dfDiffs$gene== "CCNE1"]="OG"


## Cases per sample
stats=dfDiffs%>%group_by(model,gene,gene_status,consequence)%>%summarise(n=n())

# Frequency per model
stats$freq=NA
for(m in unique(CN.TABLE$model)){
    subs=length(unique(CN.TABLE$sample[CN.TABLE$model==m]))
    stats$freq[stats$model==m]=stats$n[stats$model==m]/subs
}


# Filter by GoI
stats_filt=stats[stats$gene%in%CANCER.GENES$name[CANCER.GENES$Status!="Unknown"],]
stats_filt=stats_filt[stats_filt$gene%in%CANCER.GENES$name[CANCER.GENES$Status!="Ambivalent"],] #538 genes
write.table(stats_filt,file.path(OUTPUT_DIR,"stats_driver_alterations_Minussi.txt"),quote = F, sep="\t", col.names = T, row.names = F)

stats_filt_unique=stats_filt[stats_filt$freq<0.3,]



