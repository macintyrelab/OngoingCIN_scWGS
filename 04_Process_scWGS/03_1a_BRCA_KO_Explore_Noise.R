################################################################################
## This script is to explore noise in single-cell sequencing data from BRCA1/2 KOs
## Started: 01/08/24
## Updated: 02/06/26
##
################################################################################

### LIBRARIES
rm(list=ls(all=TRUE))

library(this.path)
library(data.table)
library(QDNAseq)
library(GenomicRanges)
library(ComplexHeatmap)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(lemon)
library(RColorBrewer)
library(lsa)
library(dplyr)

### PATHS
BASE=dirname(dirname(this.path()))
SCELL_DIR=dirname(this.path())
INPUT_DIR=file.path(BASE,"data")
INTER_DIR=file.path(SCELL_DIR,"intermediate")
dir.create(INTER_DIR, showWarnings = FALSE, recursive = TRUE)
# Create directories for outputs
OUTPUT_DIR=file.path(SCELL_DIR,"outputs")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
PLOTS_DIR=file.path(SCELL_DIR,"plots/brca_ko")
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

# For evaluation of observed noise
sdDiffTrim <- function(x, ..., trim=0.001, scale=TRUE) {
    if (scale)
        x <- x / mean(x, na.rm=TRUE)
    matrixStats::sdDiff(x, ..., trim=trim)
}
sdFUN <- match.fun(sdDiffTrim)


### LOAD INPUT DATA --> raw data will be available at publication
# Profiles all cells
sc.files.picoplex=dir(paste0(INPUT_DIR,"/KOs/picoplex"),full=T,pattern="absoluteCN")
sc.files.dlp=dir(paste0(INPUT_DIR,"/KOs/dlp"),full=T,pattern="absoluteCN")
# Profiles drug screening
sc.files.dlp367=dir(paste0(INPUT_DIR,"/KOs/dlp/DLP367"),full=T,pattern="absoluteCN") #TP53
sc.files.dlp369=dir(paste0(INPUT_DIR,"/KOs/dlp/DLP369"),full=T,pattern="absoluteCN") #BRCA1
sc.files.dlp371=dir(paste0(INPUT_DIR,"/KOs/dlp/DLP371"),full=T,pattern="absoluteCN") #BRCA2


################################################################################
#### 0. GET NOISE FROM SINGLE-CELL PROFILES ####################################

### 1. Powered cells check noise (are downsampled but I excluded those with <15nrpbpcc)
## Picoplex
SCELLS.PPLEX=as.data.frame(do.call(rbind,lapply(sc.files.picoplex,function(thisFile){
    print(basename(thisFile))
    
    # read file
    cnobject=readRDS(thisFile)
    
    # ncells
    ncells=as.numeric(length(cnobject))
    
    # noise
    if(length(cnobject)>0){
        out=c()
        for(c in 1:length(cnobject)){
            cn=cnobject[[c]]
            copynumber=assayDataElement(cn, "copynumber")
            noise=apply(copynumber, 2, sdFUN, na.rm=TRUE)
            out=rbind(out,c(noise,gsub(".markdup","",names(noise)),
                            gsub("_absoluteCN_100kb.rds","",basename(thisFile))))
            
        }
        out=as.data.frame(out)
        colnames(out)=c("noise","cell","sample")
        out$sample=gsub("sc_","",out$sample)
        out$tech="picoplex"
        return(out)
    }
})))
SCELLS.PPLEX[,1]=as.numeric(SCELLS.PPLEX[,1])
SCELLS.PPLEX$cell[1:nrow(SCELLS.PPLEX)]=sapply(1:nrow(SCELLS.PPLEX),function(x)strsplit(SCELLS.PPLEX$cell, '_')[[x]][1])

# Add info
#add model info in cell lines
info=rbind(fread(paste0(INPUT_DIR,"/KOs/sc_brca2_picoplex_15082024.txt")),
           fread(paste0(INPUT_DIR,"/KOs/sc_brca1_screen_picoplex_15082024.txt")),
           fread(paste0(INPUT_DIR,"/KOs/sc_brca1_picoplex_12062024.txt")))
colnames(info)[2]="cell"
SCELLS.PPLEX=left_join(SCELLS.PPLEX,info,by="cell")

## DLP
sc.files.dlp = c(sc.files.dlp, sc.files.dlp367, sc.files.dlp369, sc.files.dlp371)
SCELLS.DLP=as.data.frame(do.call(rbind,lapply(sc.files.dlp,function(thisFile){
    print(basename(thisFile))
    
    # read file
    cnobject=readRDS(thisFile)
    
    # ncells
    ncells=as.numeric(length(cnobject))
    
    # noise
    if(length(cnobject)>0){
        out=c()
        for(c in 1:length(cnobject)){
            cn=cnobject[[c]]
            copynumber=assayDataElement(cn, "copynumber")
            noise=apply(copynumber, 2, sdFUN, na.rm=TRUE)
            out=rbind(out,c(noise,gsub(".markdup","",names(noise)),
                            gsub("_absoluteCN_100kb.rds","",basename(thisFile))))
            
        }
        out=as.data.frame(out)
        colnames(out)=c("noise","cell","sample")
        out$sample=gsub("sc_","",out$sample)
        out$tech="dlp"
        return(out)
    }
})))
SCELLS.DLP[,1]=as.numeric(SCELLS.DLP[,1])
DLP22=SCELLS.DLP[SCELLS.DLP$sample=="scBRCA2_DLP22",]
DLP22$cell[1:nrow(DLP22)]=sapply(1:nrow(DLP22),function(x)strsplit(DLP22$cell, '_')[[x]][2])
SCELLS.DLP=rbind(DLP22,SCELLS.DLP[SCELLS.DLP$sample!="scBRCA2_DLP22",])

# Add info
info=rbind(fread(paste0(INPUT_DIR,"/KOs/sc_info_dlp75-76_20250519.txt")),
           fread(paste0(INPUT_DIR,"/KOs/sc_info_dlp367-372_drug_screenings.txt")))
colnames(info)[2]="cell"
SCELLS.DLP$model=info$model[match(SCELLS.DLP$cell,info$cell)]
SCELLS.DLP$model[SCELLS.DLP$sample=="scBRCA2_DLP22"]="BRCA2"

#### 1. COMPARE NOISE FROM DIFFERENT FITTING APPROACHES ########################
SCELLS=rbind(SCELLS.PPLEX,SCELLS.DLP)
noise.stats.approach=SCELLS%>%group_by(tech)%>%summarise(n=n(),mean=mean(noise),sd=sd(noise),median=median(noise),q60=quantile(noise,0.6),
                                                         q70=quantile(noise,0.7),q75=quantile(noise)[4])

noise.stats=SCELLS%>%group_by(model,tech)%>%summarise(n=n(),mean=mean(noise),sd=sd(noise),median=median(noise),q60=quantile(noise,0.6),
                                                       q70=quantile(noise,0.7),q75=quantile(noise)[4])

# Plotting noise
thresh=0.19 #Minussi data
thresh1=0.26 #Picoplex
thresh2=0.21 #DLP 30 nrpbpcc

png(file=paste0(PLOTS_DIR, "/BoxPlot_noise_BRCA_KOs.png"), width = 120/25.4, height = 80/25.4, units = "in", res = 300)
p=ggplot(SCELLS, aes(x=model, y=noise)) + 
    ylim(0.1,0.6) +
    geom_hline(yintercept = thresh, linetype = "dashed", colour = "red", size = 0.5) +
    geom_hline(yintercept = thresh1, linetype = "dashed", colour = "green", size = 0.5) +
    geom_hline(yintercept = thresh2, linetype = "dashed", colour = "blue", size = 0.5) +
    geom_boxplot(outlier.shape = NA, color="grey30", fill=NA)+
    geom_point(aes(fill = tech), pch = 21, size = 1, col = "black", 
               position = position_jitter(w = 0.2, h = 0), alpha = 0.8) +
    labs(subtitle = "Comparing noise across data processing samples - picoplex & dlp") +
    xlab("") +
    ylab("Noise (observed SD)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
          axis.text.y = element_text(size = 5),
          legend.position = "none")
print(p)
dev.off()


#### 2. GET THRESHOLD FOR EXCLUDING DIVING CELLS ###############################
## Picoplex
SCELLS.DIV.PPLEX=as.data.frame(do.call(rbind,lapply(sc.files.picoplex,function(thisFile){
    print(basename(thisFile))
    
    ## Step 1: Process data
    # read file
    cnobject=readRDS(thisFile)
    
    # filter out noisy profiles
    if(length(cnobject)>0){
        segTabs=c()
        for(c in 1:length(cnobject)){
            cn=cnobject[[c]]
            copynumber=assayDataElement(cn, "copynumber")
            noise=apply(copynumber, 2, sdFUN, na.rm=TRUE)
            if(noise<=0.27){ #the highest threshold that I will use later on + 0.01
                thisClon=colnames(copynumber)
                segTab=CNpare::getSegTable(cn[,thisClon])
                segTab$sample=thisClon
                segTab[,2:4]=sapply(2:4,function(x)as.numeric(segTab[,x]))
                segTabs=rbind(segTabs,segTab)
            }
        }
        nsegs=segTabs%>%group_by(sample)%>%summarise(n=n())
        thresh=mean(nsegs$n)+(2*sd(nsegs$n))
        out=c(thresh,mean(nsegs$n),sd(nsegs$n),gsub("_absoluteCN_100kb.rds","",basename(thisFile)))
        return(out)
    }
})))
colnames(SCELLS.DIV.PPLEX)=c("thresh","mean","sd","sample")
SCELLS.DIV.PPLEX[,1:3]=sapply(1:3,function(x)as.numeric(SCELLS.DIV.PPLEX[,x]))
SCELLS.DIV.PPLEX$tech="picoplex"

## DLP
SCELLS.DIV.DLP=as.data.frame(do.call(rbind,lapply(sc.files.dlp,function(thisFile){
    print(basename(thisFile))
    
    ## Step 1: Process data
    # read file
    cnobject=readRDS(thisFile)
    
    # filter out noisy profiles
    if(length(cnobject)>0){
        segTabs=c()
        for(c in 1:length(cnobject)){
            cn=cnobject[[c]]
            copynumber=assayDataElement(cn, "copynumber")
            noise=apply(copynumber, 2, sdFUN, na.rm=TRUE)
            if(noise<=0.21){ #the highest threshold that I will use later on + 0.01
                thisClon=colnames(copynumber)
                segTab=CNpare::getSegTable(cn[,thisClon])
                segTab$sample=thisClon
                segTab[,2:4]=sapply(2:4,function(x)as.numeric(segTab[,x]))
                segTabs=rbind(segTabs,segTab)
            }
        }
        nsegs=segTabs%>%group_by(sample)%>%summarise(n=n())
        thresh=mean(nsegs$n)+(2*sd(nsegs$n))
        out=c(thresh,mean(nsegs$n),sd(nsegs$n),gsub("_absoluteCN_100kb.rds","",basename(thisFile)))
        return(out)
    }
})))
colnames(SCELLS.DIV.DLP)=c("thresh","mean","sd","sample")
SCELLS.DIV.DLP[,1:3]=sapply(1:3,function(x)as.numeric(SCELLS.DIV.DLP[,x]))
SCELLS.DIV.DLP$tech="dlp"

#Save
SCELLS.DIV=rbind(SCELLS.DIV.PPLEX,SCELLS.DIV.DLP)
saveRDS(SCELLS.DIV, paste0(INTER_DIR, "/3.0_BRCA_KO_Sphase_THRESH.rds"))
