################################################################################
## This WORKFLOW documents the derivation CIN signatures which can be genomically mapped
## to a specific copy number change
## Started: 26/07/22
## Updated: 27/07/24
##
## Here I'm going to derive signatures using two pan-cancer features (SEGSIZE and CHANGEPOINT) 
## and one new feature (BPWINDOW) using TCGA SNP6 data. This new feature recapitulates
## information of clustering/oscillating events.
## The modification of signature encoding is needed for matching each event to a signature, since
## the old approach includes features getting information from the genomic context. In this case, 
## each event is isolated from the genomic context.
##
## Then, I'll correlate new and original CINSigs activities to see if our new procedure 
## collect the different CIN types as our compendium of CIN signatures (Drews et al. 2022)
##
## As stated before, we also modified the changepoint feature encoding to enable a 
## better event-to-signature mapping. Therefore, the gold-standard exposures are those obtained by
## requantifying the compedium of signatures in the SxC input matrix by applying the
## original 5 features but with the new changepoint encoding (see 00_WORKFLOW)
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
library(reshape2)

### PLOTTING THEME
theme_set(theme_tufte(base_size = 6, base_family = "ArialMT"))
theme_update(text = element_text(size = 6),
             axis.text = element_text(size = 6))


### PATHS
BASE=dirname(dirname(this.path()))
NEWSIGS_DIR=dirname(this.path())
INPUT_DIR=file.path(BASE,"data")
INTER_DIR=file.path(NEWSIGS_DIR,"intermediate")
OUTPUT_DIR=file.path(NEWSIGS_DIR,"outputs")
PLOTS_DIR=file.path(NEWSIGS_DIR,"plots")

### FUNCTIONS
source(file.path(BASE, "scripts/newCNSigs_functions.R"))

### LOAD INPUT DATA
# All data needed for running this script can be found downloaded from https://github.com/markowetzlab/Drews2022_CIN_Compendium
# Please, download all files and include them in a /data folder within the BASE path 

# Segment tables of samples with detectable CIN
TCGA.TABLE=readRDS(paste0(INPUT_DIR, "/tcga/0_TCGA_Segments_dCIN.rds"))
TCGA.TABLE$sample=factor(TCGA.TABLE$sample)
# Mixture models
MIXMODELS=readRDS(paste0(INPUT_DIR, "/Mixmodels_merged_components.rds"))
# Gold-standard exposures => obtained from running "00_WORKFLOW_Derive_GoldSignatures_newCHP.R" 
GOLD.CNSIGS.EXPS=readRDS(paste0(OUTPUT_DIR, "/0_Signature_17CINSigs_Activities_THRESH95_TCGA_newCHP.rds"))
# Pan-cancer signatures definitions
GOLD.CNSIGS.DEFS=readRDS(paste0(INPUT_DIR, "/tcga/Signature_Compendium_v5_Cosine-0.74_Signatures_NAMESAPRIL21.rds"))


################################################################################
################# Chapter 1: Fit mixture models ################################
### Step 1.1: Extract copy-number features
# In this case, I extract only two pan-cancer features --> SEGSIZE and CHANGEPOINT
# Additionally, I also extract an additional feature that informs about clustered/oscillating events
lTCGA.TABLE=split(TCGA.TABLE, TCGA.TABLE$sample) #6335 samples
RMNORM=TRUE #remove normals
CORES=1
SIZE=5000000

TCGA.ECNF=extractThreeCopynumberFeatures.new(lTCGA.TABLE, cores=CORES, prePath=INPUT_DIR, rmNorm=RMNORM, size=SIZE)
saveRDS(TCGA.ECNF, file=paste0(INTER_DIR,"/1.1_TCGA_ThreeECNF.rds"))

### Step 1.2: Fit and choose optimal mixture models
# First, I need to explore the data distribution of this new feature. To do so, I run the "01a_Explore_ECNF.Rmd" script
# This analysis indicates that 10% segments have not breakpoints within a 5MB window (0), with another 30% encoding 1's and 25% encoding 2's. 
# Only 10% of segments has more than 5 breakpoints within a 5MB window, being extreme values represented only in 1%.

# Then, I run the "01b_Fit_MixtureModel_BPwindow.R" script for performing the mixture modeling.
# The new feature is modeled by poisson mixtures. To have a correct representation of all data,
# I define 7 components representing from 0 to 6 counts. Then, I applied mixture modeling to the 
# remaining 4% of the data and derived 2 additional components, resulting in a total of 9 components.

# Finally, I combine the mixture modeling solution in one R object
FEATS=c("segsize","changepoint","bpwindow")
MODELS=MIXMODELS[names(MIXMODELS)%in%FEATS]
MODELS[["bpwindow"]]=readRDS(paste0(NEWSIGS_DIR, "/intermediate/MIXMODEL_BPWINDOW.rds"))
saveRDS(MODELS, file=paste0(INPUT_DIR,"/NEW_MIXMODELS.rds"))


################################################################################
############ Chapter 2: Generate input matrix for NMF ##########################

### Step 2.1: Calculate probability of each segment to belong to each category & create count matrix
TCGA.ECNF=readRDS(paste0(INTER_DIR,"/1.1_TCGA_ThreeECNF.rds"))
MODELS=readRDS(paste0(INPUT_DIR,"/NEW_MIXMODELS.rds"))
FEATS=names(MODELS)
UNINFPRIOR=TRUE

lMats=lapply(FEATS, function(thisFeature) {
    print(thisFeature)
    thisEcnf=TCGA.ECNF[[thisFeature]]
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

# Create input sum-of-posterior matrix
SxC=do.call(cbind, lMats)
saveRDS(SxC, file=paste0(INTER_DIR, "/1.2_TCGA_SxC_newCNSigs.rds"))


################################################################################
############ Chapter 3: Derive signatures using the new approach ###############
# Since we have the exposures of each TCGA sample to the 17 copy-number signatures, 
# here we used the YAPSA package to get signature definitions for the new feature

## YAPSA needs:
# Left matrix W        sigCatalogue        signature (rows) by components (cols)  <= WANT
# Right matrix H       expCatalogue        samples (rows) by signatures (cols)    <= HAVE
# Full matrix V        mutCatalogue        samples (rows) by components (cols)    <= HAVE

set.seed(7777)
SxC=readRDS(file=paste0(INTER_DIR, "/1.2_TCGA_SxC_newCNSigs.rds"))

DEFS=YAPSA::LCD(SxC,YAPSA:::normalize_df_per_dim(GOLD.CNSIGS.EXPS,1))
DEFS=as.matrix(t(apply(DEFS,1,function(x){x/sum(x)})))
saveRDS(DEFS, file=paste0(OUTPUT_DIR, "/1_NewSignatures_Definitions.rds"))

# plot definitions
pdf(file=paste0(PLOTS_DIR, "/1.3_TCGA_NewSignatures_Definitions_YAPSA.pdf"), width = 8, height = 4)
Heatmap(DEFS, column_title = "NewCNSigs Definitions - 3 features with 41 components (22+10+9)", cluster_rows = FALSE, cluster_columns = FALSE)
dev.off()


# We compute the cosine similarity across new signature definition to test if
# each signature represent a different pattern of copy number alterations
sim=cosineSimForMatrices(t(DEFS))
colnames(sim)=rownames(DEFS)
rownames(sim)=rownames(DEFS)
png(file=paste0(PLOTS_DIR, "/1.3_TCGA_NewSignatures_Definitions_Similaritiess_YAPSA.png"), width = 8, height = 6, units = "in", res = 300)
Heatmap(sim, column_title = "Cosine Similarities across NewCNSigs Definitions", cluster_rows = FALSE, cluster_columns = FALSE,
        cell_fun = function(j, i, x, y, width, height, fill) {
            if(sim[i, j] > 0.74)
                grid.text(sprintf("%.3f", sim[i, j]), x, y, gp = gpar(fontsize = 6))
        })
dev.off()

################################################################################
############ Chapter 4: Re-quantify exposures using new definitions ############
# After deriving the new signature definitions, we quantify their activity in each sample.

## YAPSA needs:
# Left matrix W        sigCatalogue        signature (rows) by components (cols)  <= HAVE
# Right matrix H       expCatalogue        samples (rows) by signatures (cols)    <= WANT
# Full matrix V        mutCatalogue        samples (rows) by components (cols)    <= HAVE

## Quantify all 17 signatures
EXPS=YAPSA::LCD(t(SxC),YAPSA:::normalize_df_per_dim(t(DEFS),2))
EXPS=as.matrix(t(apply(EXPS,2,function(x){x/sum(x)}))) #normalised per sample
saveRDS(EXPS, file=paste0(OUTPUT_DIR, "/1_NewSignatures_Activities.rds"))


################################################################################
############ Chapter 5: Calculate the signature-specific threshold #############
# To trust in small signature activities, we defined signature-specific thresholds 
# based on the 1,000 simulations. As in Drews et al. (2022), I fitted a single
# Gaussian distribution to the non-zero values and took the quantile 0.95 as a threshold

### Step 5.1: Simulate feature distributions with noise
### Step 5.2: Derive SxC matrices
### Step 5.3: Calculate signature activities
### Step 5.4: Visualize results & save activity thresholds
# These four steps are done in the cluster by running the `01c_Threshold_SigActivities.R`
# This script has been slightly modified from the one used in the pan-cancer paper


################################################################################
############ Chapter 6: Apply signature-specific thresholds ####################

# We set to zero those activity values lower than the signature-specific threshold
EXPS=readRDS(file=paste0(OUTPUT_DIR, "/1_NewSignatures_Activities.rds"))
threshs=fread(file.path(NEWSIGS_DIR, "MCsimulations/MC3_Boxplots_per_sig_newCNSigs_fullTCGA_1000sims_10pGaussian_10pSamplePoisson.txt"))

allSigs=colnames(EXPS)
EXPS.CORRECTED = do.call(cbind,lapply(1:length(allSigs), function(thisSig){
    sig=allSigs[thisSig]
    acts=EXPS[,thisSig]
    thrs=threshs$Thresh_ZeroGMM_0.05[threshs$Sig==sig]
    acts[acts<thrs]=0
    return(acts)
}))
colnames(EXPS.CORRECTED)=allSigs
EXPS.CORRECTED=t(apply(EXPS.CORRECTED,1,function(x) x/sum(x)))
saveRDS(EXPS.CORRECTED,file=paste0(OUTPUT_DIR, "/1_NewSignatures_Activities_THRESH95.rds"))

# Now, we perform a correlation between the pan-cancer exposures (gold-standard)
# and the activity levels of the new signatures. This will allow us to confirm that
# the new signature encoding collects the different genomic footprints leave by each CIN type

GOLD.CNSIGS.EXPS=GOLD.CNSIGS.EXPS[row.names(GOLD.CNSIGS.EXPS)%in%row.names(EXPS.CORRECTED),]
EXPS.CORRECTED=EXPS.CORRECTED[order(row.names(GOLD.CNSIGS.EXPS)),]
correlation<-do.call(rbind,lapply(1:ncol(EXPS.CORRECTED), function(thisSig){
    corr=cor(EXPS.CORRECTED[,thisSig], GOLD.CNSIGS.EXPS, method="spearman")}))
row.names(correlation)<-colnames(EXPS.CORRECTED)

# plot correlations
pdf(file=paste0(PLOTS_DIR, "/1.6_SuppFig5_TCGA_Exposure_SpearmanCorrelations_NewSignatures_vs_GoldSignatures.pdf"), width = 200/25.4, height = 150/25.4)
p=Heatmap(correlation, column_title = "Spearman's rho between Exposures of original vs new signatures", cluster_rows = FALSE, cluster_columns = FALSE,
          cell_fun = function(j, i, x, y, width, height, fill) {
              grid.text(sprintf("%.3f", correlation[i, j]), x, y, gp = gpar(fontsize = 6))
              })
print(p)
dev.off()


# We calculate pan-cancer prevalence we counted in how many samples each signature activity
# was above its specific threshold
dtPrev = as.data.frame(do.call(rbind, lapply(1:length(allSigs), function(thisSig){
    acts=EXPS.CORRECTED[,thisSig]
    nZero=length(acts[acts==0])
    prev=(length(acts)-nZero)/length(acts)
    out=c(allSigs[thisSig],nZero,length(acts)-nZero,prev)
    return(out)
})))
colnames(dtPrev)=c("Sig","zeros","non-zeros","prevalence")
dtPrev$prevalence=as.numeric(dtPrev$prevalence)

# Plot prevalence
dtPrev$Sig=factor(dtPrev$Sig, levels = allSigs)
pOut = ggplot(dtPrev, aes(x = Sig, y = prevalence)) + 
    geom_crossbar(aes(ymin = prevalence, ymax = prevalence), width = 0.6) + 
    scale_y_continuous(labels = scales::percent, limits=c(0,1)) + labs(x = "Signatures", y = "Prevalence") +
    theme(axis.line = element_line(size = 1), axis.ticks = element_line(size = 1),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

png(paste0(PLOTS_DIR, "/1.6_NewSignatures_Prevalence_CorrectedExposures.png"), width = 120, height = 60, units = "mm", res = 300)
print(pOut)
dev.off()


##################### END OF THE SCRIPT ########################################