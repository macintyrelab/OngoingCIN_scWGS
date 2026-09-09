################################################################################
## This WORKFLOW documents the derivation of genome-level signatures using the new 
## changepoint encoding
## Started: 19/01/23
## Updated: 27/07/25
##
## Procedure: 
## Here I'm going to apply the original approach (5 features) but using the new
## changepoint encoding to quantify our 17 genome-level signatures (Drews et al. 2022)
## in the TCGA SNP6 profiles.
## I quantify exposures using YAPSA and the original definitions. Then, I compare 
## sample exposures with the original exposures to see if there are significant differences. 
##
## Main conclusion:
## The modification of the changepoint encoding do not impact on the ability to detect
## our compendium of CIN signatures. 
##
## I use YAPSA to quantify the 17 signatures with the new changepoint encoding 
## to generate our set of signature activities for this new feature encoding.
## These activities will be used as gold-standard for all projects applying the new
## feature encoding.
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
# Create directories for outputs
INTER_DIR=file.path(NEWSIGS_DIR,"intermediate")
dir.create(INTER_DIR, showWarnings = FALSE, recursive = TRUE)
OUTPUT_DIR=file.path(NEWSIGS_DIR,"outputs")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
PLOTS_DIR=file.path(NEWSIGS_DIR,"plots")
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)


### LOAD INPUT DATA
# All data needed for running this script can be found downloaded from https://github.com/markowetzlab/Drews2022_CIN_Compendium
# Please, download all files and include them in a /data folder within the BASE path 

# Segment tables of samples with detectable CIN
TCGA.TABLE=readRDS(paste0(INPUT_DIR, "/tcga/0_TCGA_Segments_dCIN.rds"))
TCGA.TABLE$sample=factor(TCGA.TABLE$sample)
# Mixture models
MODELS=readRDS(paste0(INPUT_DIR, "/Mixmodels_merged_components.rds"))
# Pan-cancer signatures definitions
GOLD.CNSIGS.DEFS=readRDS(paste0(INPUT_DIR, "/tcga/Signature_Compendium_v5_Cosine-0.74_Signatures_NAMESAPRIL21.rds"))


### FUNCTIONS
source(file.path(BASE, "scripts/main_functions.R"))
source(file.path(BASE, "scripts/newCNSigs_functions.R"))


################################################################################
################# Chapter 1: Fit mixture models ################################
### Step 1.1: Extract copy-number features
lTCGA.TABLE=split(TCGA.TABLE, as.character(TCGA.TABLE$sample))
RMNORM=TRUE #remove normals
CORES=1

lECNF=extractCopynumberFeatures.newchp(lTCGA.TABLE, cores = CORES, prePath = paste0(INPUT_DIR,"/"), rmNorm = RMNORM)
saveRDS(lECNF, paste0(INTER_DIR, "/0.1_TCGA_ECNF_newCHP.rds"))

### Step 1.2: Fit and choose optimal mixture models
# I'm going to use the same components than in the original approach


################################################################################
############ Chapter 2: Generate input matrix for NMF ##########################

### Step 2.1: Calculate probability of each segment to belong to each category & create count matrix
lECNF=readRDS(paste0(INTER_DIR, "/0.1_TCGA_ECNF_newCHP.rds"))
UNINFPRIOR=TRUE
FEATS=names(MODELS)

lMats = lapply(FEATS, function(thisFeature) {
    
    print(thisFeature)
    thisEcnf = lECNF[[thisFeature]]
    thisModel = MODELS[[thisFeature]]
    
    dat = as.numeric(thisEcnf[,2] )
    # We want a posterior, hence likelihood (density) times prior (weight)
    if( thisFeature %in% c("bp10MB","bpchrarm","osCN") ) {
        # Poisson model
        print("Poisson-based posterior")
        if(UNINFPRIOR){
            postDatUnscaled = sapply(1:nrow(thisModel), function(x) dpois(x = dat, lambda = thisModel[[x,"Mean"]]) )
        } else {
            postDatUnscaled = sapply(1:nrow(thisModel), function(x) dpois(x = dat, lambda = thisModel[[x,"Mean"]]) * thisModel[[x, "Weight"]] )
        }
        
    } else {
        # Gaussian model
        print("Gauss-based posterior")
        if(UNINFPRIOR){
            postDatUnscaled = sapply(1:nrow(thisModel), function(x) dnorm(x = dat, mean = thisModel[[x,"Mean"]], sd = thisModel[[x,"SD"]]) )
        } else {
            postDatUnscaled = sapply(1:nrow(thisModel), function(x) dnorm(x = dat, mean = thisModel[[x,"Mean"]], sd = thisModel[[x,"SD"]]) * thisModel[[x, "Weight"]] )
        }
    }
    
    # Normalise densities to probabilities
    postDatScaled = data.frame(postDatUnscaled / rowSums(postDatUnscaled))
    postDatScaled$Sample = thisEcnf[,1]
    matSxC = aggregate(. ~ Sample, postDatScaled, sum)
    rownames(matSxC) = matSxC$Sample
    matSxC$Sample = NULL
    matSxC = as.matrix(matSxC)
    
    # Should be sorted but just to be sure
    matSxC = matSxC[ ,order(thisModel[["Mean"]])]
    colnames(matSxC) = paste0(thisFeature, 1:ncol(matSxC))
    
    return(matSxC)
})
names(lMats)=FEATS
# Create input sum-of-posterior matrix
SxC=do.call(cbind, lMats)
saveRDS(SxC, file=paste0(INTER_DIR, "/0.2_TCGA_SxC_newCHP.rds"))


################################################################################
############ Chapter 3: Quantify signatures using the new approach ##############
# I applied YAPSA to quantify signatures after applying the new changepoint encoding
# and also to derive de novo genome-level signatures in the modified 5 feature space.
# This analysis is done by running the "00a_Quantify_TCGA_OriginalCINSigs_newCHP.R" script.
# Results show that we recapitulate all 17 signatures with this new changepoint encoding,
# so it's not impacting on our signature compendium


################################################################################
######## Chapter 4: Quantify signatures detected with the new encoding #########
set.seed(7777)
V = t(readRDS(paste0(INTER_DIR, "/0.2_TCGA_SxC_newCHP.rds")))
W = GOLD.CNSIGS.DEFS

# Sanity check signature matrix
if(nrow(W) < ncol(W)) W = t(W)

# Check order of components and fix if necessary
if(!identical(rownames(W), rownames(V))) {
    W = W[match(rownames(V), rownames(W)), ]
}

## YAPSA needs:
# Left matrix W        sigCatalogue        signature (rows) by components (cols)  <= HAVE
# Right matrix H       expCatalogue        samples (rows) by signatures (cols)    <= WANT
# Full matrix V        mutCatalogue        samples (rows) by components (cols)    <= HAVE

EXPS=YAPSA::LCD(V,YAPSA:::normalize_df_per_dim(W,2))
EXPS=as.matrix(t(apply(EXPS,2,function(x){x/sum(x)}))) #normalised per sample
saveRDS(EXPS, file=paste0(INTER_DIR, "/0.4_Signature_17CINSigs_Activities_TCGA_newCHP.rds"))


################################################################################
############ Chapter 5: Calculate the signature-specific threshold #############
# To trust in small signature activities, we defined signature-specific thresholds 
# based on the 1,000 simulations. As in Drews et al. (2022), I fitted a single
# Gaussian distribution to the non-zero values and took the quantile 0.95 as a threshold

### Step 5.1: Simulate feature distributions with noise
### Step 5.2: Derive SxC matrices
### Step 5.3: Calculate signature activities
### Step 5.4: Visualize results & save activity thresholds
# These four steps are done in the cluster by running the `00b_Threshold_17CINSigActivities_newCHP.R`
# This script has been slightly modified from the one used in the Nature paper 

### Step 5.5: Apply signature-specific thresholds
# We set to zero those activity values lower than the signature-specific threshold
EXPS=readRDS(file=paste0(INTER_DIR, "/0.4_Signature_17CINSigs_Activities_TCGA_newCHP.rds"))
threshs=fread(file.path(NEWSIGS_DIR, "MCsimulations/MC3_Boxplots_per_sig_17CINSigs_newCHP_fullTCGA_1000sims_10pGaussian_10pSamplePoisson.txt"))

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
saveRDS(EXPS.CORRECTED,file=paste0(OUTPUT_DIR, "/0_Signature_17CINSigs_Activities_THRESH95_TCGA_newCHP.rds"))

## These EXPS are going to be used as the gold-standard for all projects using
## the new changepoint encoding!!

##################### END OF THE SCRIPT ########################################
