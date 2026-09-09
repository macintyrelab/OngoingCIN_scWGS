################################################################################
## This script quantify our compendium of copy number signatures but after applying the new
## feature encoding for the changepoint feature
## 
## Procedure:
## I quantify exposures in samples using YAPSA and the original definitions. 
## Then, I will compare sample exposures with the original exposures to see if there are
## significant differences. In case exposures are similar, this will be my new gold-standard
## exposures for future comparisons
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
library(mclust)

### PATHS
BASE=dirname(dirname(this.path()))
NEWSIGS_DIR=dirname(this.path())
INPUT_DIR=file.path(BASE,"data")
INTER_DIR=file.path(NEWSIGS_DIR,"intermediate")

### LOAD INPUT DATA
# Gold-standard definitions
GOLD.CNSIGS.DEFS=readRDS(file.path(INPUT_DIR, "tcga/Signature_Compendium_v5_Cosine-0.74_Signatures_NAMESAPRIL21.rds"))
# Gold-standard exposures
GOLD.CNSIGS.EXPS=readRDS(file.path(INPUT_DIR, "tcga/Signature_Compendium_v5_Cosine-0.74_Activities_THRESH95_NAMESAPRIL21.rds"))


################################################################################
############ Quantify original gold-standard signatures ########################
### Quantify original pan-cancer signatures after applying the new changepoint encoding
set.seed(7777)

## YAPSA needs:
# Left matrix W        sigCatalogue        signature (rows) by components (cols)  <= HAVE
# Right matrix H       expCatalogue        samples (rows) by signatures (cols)    <= WANT
# Full matrix V        mutCatalogue        samples (rows) by components (cols)    <= HAVE

V = t(readRDS(paste0(INTER_DIR, "/0.2_TCGA_SxC_newCHP.rds")))
W = GOLD.CNSIGS.DEFS

# Sanity check signature matrix
if(nrow(W) < ncol(W)) W = t(W)

# Check order of components and fix if necessary
if(!identical(rownames(W), rownames(V))) {
    W = W[match(rownames(V), rownames(W)), ]
}

### YAPSA needs:
## Full matrix V        mutCatalogue        components (rows) by samples (cols)    <= HAVE
## Left matrix W        sigCatalogue        components (rows) by signature (cols)   <= HAVE
## Right matrix H       expCatalogue        signature (rows) by samples (cols)     <= WANT

# in_mutation_catalogue_df => NxM => N - Features, M - Samples => Component by Sample matrix
# in_signatures_df => NxL => N - Features, L - Signatures => Component by Signature matrix
Hraw = as.matrix(LCD(in_mutation_catalogue_df = V, in_signatures_df = W, in_per_sample_cutoff = 0))
H = t(apply(Hraw, 2, function(x) x/sum(x)))

### Correlation with TCGA activities
H=H[order(row.names(GOLD.CNSIGS.EXPS)),]
correlation<-do.call(rbind,lapply(1:ncol(H), function(thisSig){
    corr=cor(H[,thisSig], GOLD.CNSIGS.EXPS, method="spearman")}))
row.names(correlation)<-colnames(H)

# plot correlations
# Supplementary figure 1
png(file=paste0(NEWSIGS_DIR, "/plots/0_SuppNoteFig3_TCGA_SpearmanCorrelations_Exposures_OriginalCINSigs_gold_vs_newCHP.png"), width = 8, height = 6, units = "in", res = 300)
p=Heatmap(correlation, column_title = "Spearman's rho between CINSigs exposures using original vs new encoding", cluster_rows = FALSE, cluster_columns = FALSE,
          cell_fun = function(j, i, x, y, width, height, fill){
              grid.text(sprintf("%.3f", correlation[i, j]), x, y, gp = gpar(fontsize = 6))
              })
print(p)
dev.off()

pdf(file.path(NEWSIGS_DIR, "/plots/0_SuppNoteFig3_TCGA_Pearson_Exposures_OriginalCINSigs_gold_vs_newCHP.pdf"), height = 150/25.4, width = 200/25.4)
print(p)
dev.off()


######################### END OF SCRIPT ########################################

