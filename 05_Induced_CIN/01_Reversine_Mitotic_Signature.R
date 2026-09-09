################################################################################
## This script is to process single clone sequencing data from Watson 2024
## Started: 08/08/24
## Updated: 06/08/25
##
## In this work, normal mammary cell lines are treated with reversine to induce
## missegregation errors, and then explore evolution in vitro. 
## This data is to characterize a signature linked to mitotic errors,
## which is the mutational process that has been induced. 
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
PLOTS_DIR=file.path(INDCIN_DIR,"plots/reversine")
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


### LOAD INPUT DATA
# Profiles 
SCELLS.CN=readRDS(paste0(INPUT_DIR, "/2.1.0_Watson_Segments.rds")) # 41 clones => all 4N
# ECNF 
ECNF=readRDS(paste0(INPUT_DIR, "/2.1.0_Watson_ECNF.rds"))
# ExC matrix
ExC=readRDS(paste0(INPUT_DIR, "/2.1.1_Watson_Reversine_ExC.rds"))

# Feature components
MIXMODELS=readRDS(paste0(BASE,"/data/SC_MIXMODELS.rds"))
# Signature definitions
DEFS=readRDS(file=paste0(NEWSIGS_DIR,"/outputs/2_New_scCINSignatures_Definitions_Corrected.rds"))


################################################################################
### CHARACTERIZING SIGNATURES LINKED TO MITOTIC ERRORS #########################

# We compute the sum-of-posterior probabilities for each feature component
# This will be compared to the definitions of the signatures to see the best match

### 1. Extract feature components distributions per sample (SxC)
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
SxC=do.call(cbind, lMats)
row.names(SxC)=unique(SCELLS.CN$clone)

# Plot SxC
pdf(file=paste0(PLOTS_DIR, "/1.1_Heatmap_SxC_Reversine.pdf"), width = 250/25.4, height = 200/25.4)
p=Heatmap(SxC, column_title = "Sum-of-posterior feature components matrix", cluster_rows = FALSE, cluster_columns = FALSE)
print(p)
dev.off()


### 2. Compare with signature definitions

# Cosine similarities
cosine<-do.call(rbind,lapply(1:nrow(DEFS), function(thisSig){
    cos=lsa::cosine(DEFS[thisSig,], t(SxC[,]))}))
row.names(cosine)<-row.names(DEFS)

# Plot similarities
pdf(file=paste0(PLOTS_DIR, "/1.2_Heatmap_scCINSignatures_vs_Reversine.pdf"), width = 250/25.4, height = 200/25.4)
p=Heatmap(cosine, col = pal, column_title = "Reversine treatment \n cosine similarities -- noBPW", cluster_rows = FALSE, cluster_columns = FALSE,
          cell_fun = function(j, i, x, y, width, height, fill) {
              if(cosine[i, j] > 0.7)
                  grid.text(sprintf("%.2f", cosine[i, j]), x, y, gp = gpar(fontsize = 6))
              }
          )
print(p)
dev.off()


# Plot means
cos=rowMeans(cosine)
pdf(file=paste0(PLOTS_DIR, "/1.2_Heatmap_scCINSignatures_vs_Reversine_mean_clones.pdf"), width = 250/25.4, height = 80/25.4)
p=Heatmap(t(cos), col = pal, column_title = "Reversine treatment", cluster_rows = FALSE, cluster_columns = FALSE,
          heatmap_legend_param = list(
              title = "Cosine"
          ))
print(p)
dev.off()


### 3. Ploting definition
# reversine signature
rev.sig=colMeans(SxC)
rev.sig[1:20]=rev.sig[1:20]/sum(rev.sig[1:20])
rev.sig[21:28]=rev.sig[21:28]/sum(rev.sig[21:28])
rev.sig[29:37]=rev.sig[29:37]/sum(rev.sig[29:37])


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
ggsave(paste0(PLOTS_DIR, "/1.3_BarPlot_reversine_signature_definition.svg"), p, width = 150, height = 60, units = "mm", bg = "white")
ggsave(paste0(PLOTS_DIR, "/1.3_BarPlot_reversine_signature_definition.png"), p, width = 150, height = 60, units = "mm", bg = "white")


# The most similar signatures are related to mitotic errors!

############ END SCRIPT ! ######################################################






