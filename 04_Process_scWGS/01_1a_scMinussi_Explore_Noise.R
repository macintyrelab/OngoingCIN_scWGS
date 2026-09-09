################################################################################
## This script is to explore noise in single-cell sequencing data from
## Minussi et al. 2022
## Started: 06/06/24
## Updated: 23/01/25
##
## Here I compare the noise from different fitting approaches
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
PLOTS_DIR=file.path(SCELL_DIR,"plots/minussi")
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


### LOAD INPUT DATA
# Raw data can be downloaded from https://www.ncbi.nlm.nih.gov/sra?term=PRJNA629885
# Raw data is processed as indicated in the Methods section
# Profiles downsampled & powered cells
sc.files=dir(paste0(INPUT_DIR,"/minussi_sc"),full=T,pattern="absoluteCN")


################################################################################
#### 0. GET NOISE FROM SINGLE-CELL PROFILES FROM MINUSSI SAMPLES ###############

## Powered cells check noise (are downsampled but I excluded those with <15nrpbpcc)
SCELLS.NOISE=as.data.frame(do.call(rbind,lapply(sc.files,function(thisFile){
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
            out=rbind(out,c(noise,gsub("_R1.markdup","",names(noise)),
                            gsub("_absoluteCN_100kb.rds","",basename(thisFile))))
            
        }
        out=as.data.frame(out)
        colnames(out)=c("noise","cell","sample")
        out$sample=gsub("sc_","",out$sample)
        return(out)
    }
})))
SCELLS.NOISE[,1]=as.numeric(SCELLS.NOISE[,1])

#### 1. COMPARE NOISE FROM DIFFERENT FITTING APPROACHES ########################

noise.stats.approach=SCELLS.NOISE%>%summarise(n=n(),mean=mean(noise),sd=sd(noise),
                                              median=median(noise),q60=quantile(noise,0.6),
                                              q70=quantile(noise,0.7),q75=quantile(noise)[4])
noise.stats.samps=SCELLS.NOISE%>%group_by(sample)%>%summarise(n=n(),mean=mean(noise),sd=sd(noise),
                                                              median=median(noise),q60=quantile(noise,0.6),
                                                              q70=quantile(noise,0.7),q75=quantile(noise)[4])

thresh=0.18 

png(file=paste0(PLOTS_DIR, "/BoxPlot_noise_samples.png"), width = 100/25.4, height = 60/25.4, units = "in", res = 300)
p=ggplot(SCELLS.NOISE, aes(x=sample, y=noise)) + 
    geom_hline(yintercept = thresh, linetype = "dashed", colour = "red", size = 0.5) +
    geom_boxplot(outlier.shape = NA, color="grey30", fill=NA)+
    geom_point(aes(fill=sample),pch = 21, size = 1, col = "black", 
               position = position_jitter(w = 0.2, h = 0), alpha = 0.8) +
    labs(subtitle = "Noise across cells/samples") +
    xlab("") +
    ylab("Noise (observed SD)") +
    theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 5, angle=90),
          axis.text.y = element_text(size = 5),
          legend.position = "none")
print(p)
dev.off()


#### 2. GET THRESHOLD FOR EXCLUDING DIVING CELLS ###############################

SCELLS.DIV=as.data.frame(do.call(rbind,lapply(sc.files,function(thisFile){
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
            if(noise<=0.20){ #the threshold used + 0.01
                thisClon=colnames(copynumber)
                segTab=CNpare::getSegTable(cn[,thisClon])
                segTab$sample=thisClon
                segTab[,2:4]=sapply(2:4,function(x)as.numeric(segTab[,x]))
                segTabs=rbind(segTabs,segTab)
            }
        }
        if(!is.null(segTabs)){
            nsegs=segTabs%>%group_by(sample)%>%summarise(n=n())
            thresh=mean(nsegs$n)+2*sd(nsegs$n) 
            out=c(thresh,mean(nsegs$n),sd(nsegs$n),gsub("_absoluteCN_100kb.rds","",basename(thisFile)))
            return(out)
        }
    }
})))
colnames(SCELLS.DIV)=c("thresh","mean","sd","sample")
SCELLS.DIV[,1:3]=sapply(1:3,function(x)as.numeric(SCELLS.DIV[,x]))
saveRDS(SCELLS.DIV, paste0(INTER_DIR, "/1.0_scMinussi_Powered_Sphase_THRESH.rds"))
