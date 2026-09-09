################################################################################
## This script is to derive single-event signatures in a single-cell scenario
## Started: 21/04/24
## Updated: 28/07/25
##
## Here, I use all segments to derive single-event signatures to be use for dissecting
## ongoingCIN using single-cell and/or single-cell clone data. As a reference cohort
## for deriving this new signatures I am using the TCGA cohort.
##
## To derive signatures, we use the gold-standard signature activities derived using raw copy
## number profiles and the 5 feature encoding with the new changepoint (derive in FOCAS manuscript).
## I input the sample-by-signatures and the sample-by-components matrix in YAPSA 
## to derive the definitions of the single-cell signatures. Then, YAPSA is applied 
## to quantify sample-per-signature activities using the new single-cell signatures
## definitions. Finally, I compare signatures and activities generated to the gold-standard
## signature compendium to see the signatures I can recapitulate in a clonal scenario.
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
NEWSIGS_DIR=file.path(BASE,"01_Derive_New_CINSignatures")
SC_DIR=dirname(this.path())
INPUT_DIR=file.path(BASE,"data")
INTER_DIR=file.path(SC_DIR,"intermediate")
dir.create(INTER_DIR, showWarnings = FALSE, recursive = TRUE)
OUTPUT_DIR=file.path(SC_DIR,"outputs")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
PLOTS_DIR=file.path(SC_DIR,"plots")
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

### PARAMETERS
REMOVEQUIET=TRUE # Just to check that all samples have at least 20 CNAs. If not remove.
RMNORM=TRUE
CORES=1
WIGGLE=0.1
IGNOREDELS=FALSE
SIZE=5000000
UNINFPRIOR=TRUE


### LOAD INPUT DATA
## TCGA
TCGA.TABLE=rbind(readRDS(paste0(INPUT_DIR, "/tcga/0_TCGA_Segments_dCIN.rds")))
TCGA.TABLE$sample=as.character(TCGA.TABLE$sample)
ECNF.TCGA=readRDS(paste0(INTER_DIR, "/1.1.1_TCGA_Segments_100kb_ECNF.rds"))

# Mixmodels for single-cell signatures
MIXMODELS=readRDS(paste0(INPUT_DIR,"/SC_MIXMODELS.rds"))

# Gold-standard signatures using new changepoint encoding
DEFS.gold=readRDS(file=file.path(INPUT_DIR,"tcga/Signature_Compendium_v5_Cosine-0.74_Signatures_NAMESAPRIL21.rds"))
EXPS.gold=readRDS(file=file.path(NEWSIGS_DIR,"outputs/0_Signature_17CINSigs_Activities_THRESH95_TCGA_newCHP.rds"))


################################################################################
#### 1. GENERATE SAMPLE-BY-COMPONENT MATRIX ####################################
# We generate the input matrix for YAPSA using the 3-feature single-cell encoding space 

### 1.1. Selecting samples with CIN
# This is needed for trusting in the deconvolution procedure
samps=row.names(EXPS.gold) #6335 samples

### 1.2. Generating input matrix
FEATS=names(MIXMODELS)

lMats=lapply(FEATS, function(thisFeature) {
    print(thisFeature)
    thisEcnf=ECNF.TCGA[[thisFeature]]
    thisModel=MIXMODELS[[thisFeature]]
    
    # Select samples with >= NEVENTS
    thisEcnf=thisEcnf[thisEcnf$ID%in%samps,]
    
    # We want a posterior, hence likelihood (density) times prior (weight)
    dat=as.numeric(thisEcnf[,2])
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
    # Make sure that focal amplifications with huge number of copies have weight in the chp=8 
    # Make sure that windows with a lot of segments have weight in the bpw=8
    id = which(rowSums(is.na(postDatScaled)) == ncol(postDatScaled))
    if (thisFeature == "changepoint"){
        if(length(id)!=0){postDatScaled[id,]=matrix(rep(c(rep(0,7),1), length(id)), ncol = 8, nrow =length(id), byrow = T)}
    }
    if (thisFeature == "bpwindow"){
        if(length(id)!=0){postDatScaled[id,]=matrix(rep(c(rep(0,7),1), length(id)), ncol = 8, nrow =length(id), byrow = T)}
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
saveRDS(SxC, file=paste0(INTER_DIR, "/1.2.1_TCGA_SxC_100kb.rds"))
write.table(t(SxC), file=paste0(INTER_DIR,"/1.2.1_TCGA_SxC_100kb.txt"), quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")

# Plot input matrix
SxC=t(apply(SxC,1,function(x){x/sum(x)}))
png(file=paste0(PLOTS_DIR, "/1.2.1_Heatmap_TCGA_SxC_100kb.png"), width = 200/25.4, height = 150/25.4, units="in", res = 300)
p=Heatmap(SxC, column_title = "Sample-per-component matrix", cluster_rows = TRUE, cluster_columns = FALSE,
          row_names_gp = gpar(fontsize = 5), column_names_gp = gpar(fontsize = 6))
print(p)
dev.off()


#### 2. APPLY YAPSA FOR OBTAINING THE SIGNATURE-BY-COMPONENT MATRIX ############
# Since we have the exposures of each TCGA sample to the gold-standard CX signatures, 
# here we used the YAPSA package to get signature definitions for the new single-event 
# feature encoding

## YAPSA needs:
# Left matrix W        sigCatalogue        signature (rows) by components (cols)  <= WANT
# Right matrix H       expCatalogue        samples (rows) by signatures (cols)    <= HAVE
# Full matrix V        mutCatalogue        samples (rows) by components (cols)    <= HAVE

set.seed(7777)
SxC=readRDS(file=paste0(INTER_DIR, "/1.2.1_TCGA_SxC_100kb.rds"))
DEFS=YAPSA::LCD(SxC,YAPSA:::normalize_df_per_dim(EXPS.gold,1))
DEFS=as.matrix(t(apply(DEFS,1,function(x){x/sum(x)})))
saveRDS(DEFS, file=paste0(OUTPUT_DIR,"/2_New_scCINSignatures_Definitions.rds"))

# plot definitions
pdf(file=paste0(PLOTS_DIR, "/1.2.2_scCINSignatures_Definitions_YAPSA.pdf"), width = 200/25.4, height = 150/25.4)
Heatmap(DEFS, column_title = "scCINSigs Definitions - 3 features with 37 components (20+8+9)", cluster_rows = FALSE, cluster_columns = FALSE)
dev.off()


#### 3. RE-QUANTIFY ACTIVITIES USING NEW DEFINITIONS ###########################
# After deriving the single-cell signature definitions, we quantify their activity 
# in each TCGA sample

## YAPSA needs:
# Left matrix W        sigCatalogue        signature (rows) by components (cols)  <= HAVE
# Right matrix H       expCatalogue        samples (rows) by signatures (cols)    <= WANT
# Full matrix V        mutCatalogue        samples (rows) by components (cols)    <= HAVE

SxC=readRDS(file=paste0(INTER_DIR, "/1.2.1_TCGA_SxC_100kb.rds"))
DEFS=readRDS(file=paste0(OUTPUT_DIR,"/2_New_scCINSignatures_Definitions.rds"))
EXPS=YAPSA::LCD(t(SxC),YAPSA:::normalize_df_per_dim(t(DEFS),2))
EXPS=as.matrix(t(apply(EXPS,2,function(x){x/sum(x)}))) #normalised per sample
saveRDS(EXPS, file=paste0(OUTPUT_DIR, "/2_New_scCINSignatures_Activities.rds"))

#### 4. SIGNATURE-SPECIFIC THRESHOLD ###########################################
## Calculate the signature-specific threshold 
# To trust in small signature activities, we defined signature-specific thresholds 
# based on the 1,000 simulations. As in Drews et al. (2022), I fitted a single
# Gaussian distribution to the non-zero values and took the quantile 0.95 as a threshold

## Step 4.1: Simulate feature distributions with noise
## Step 4.2: Derive SxC matrices
## Step 4.3: Calculate signature activities
## Step 4.4: Visualize results & save activity thresholds
# These four steps are done in the cluster by running the `02a_Threshold_scCINSigActivities.R`
# This script has been slightly modified from the one used in the Nature paper 

## Step 4.5: Apply signature-specific thresholds
# We set to zero those activity values lower than the signature-specific threshold
EXPS=readRDS(file=paste0(OUTPUT_DIR, "/2_New_scCINSignatures_Activities.rds"))
threshs=fread(file.path(SC_DIR, "MCsimulations/MC3_Boxplots_per_sig_scCINSigs_fullTCGA_1000sims_10pGaussian_10pSamplePoisson.txt"))

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
saveRDS(EXPS.CORRECTED,file=paste0(OUTPUT_DIR, "/2_New_scCINSignatures_Activities_THRES95.rds"))


#### 5. ASSESSMENT OF NEW SINGLE-CELL SIGNATURES TO CAPTURE AETIOLOGIES ########
# Now, we perform a correlation between the focas activities and the activity 
# levels of the new single-cell signatures. This will allow us to confirm that
# the new signature encoding (with non-integer changepoints) collects the 
# different genomic footprints leave by each CIN type
EXPS=readRDS(file=paste0(OUTPUT_DIR, "/2_New_scCINSignatures_Activities_THRES95.rds"))

EXPS.gold=EXPS.gold[row.names(EXPS.gold)%in%row.names(EXPS),]
EXPS=EXPS[order(row.names(EXPS.gold)),]
correlation<-do.call(rbind,lapply(1:ncol(EXPS), function(thisSig){
    corr=cor(EXPS[,thisSig], EXPS.gold, method="spearman")}))
row.names(correlation)<-colnames(EXPS)

# plot correlations
pdf(file=paste0(PLOTS_DIR, "/1.2.5_SuppNoteFig3_SpearmanCorrelations_Activities_gold_vs_sc.pdf"), width = 200/25.4, height = 150/25.4)
p=Heatmap(correlation, column_title = "Spearman's rho between activities using gold vs scell encoding", cluster_rows = FALSE, cluster_columns = FALSE,
          cell_fun = function(j, i, x, y, width, height, fill) {
              grid.text(sprintf("%.3f", correlation[i, j]), x, y, gp = gpar(fontsize = 6))
          })
print(p)
dev.off()


## Main conclusions!
# We observe that CX7 is lost. It's expected since CX7 was mainly defined by subclonal
# alterations, which are not possible in a single-cell scenario
# We also see that CX6 may represent early WGD. It's a signature defined by long segments
# that has high activity in samples with CX4. This was expected since the original 
# definition of CX6 (interpretation matrix) shows that long, <1 & 2-copy changes which
# is in line with what it's observed in samples with early WGD.


## Important note!
# Our aim now is to map copy number events to signatures. For that, we need that the
# signatures are defined in the 3 feature space, but this is not the case for CX1, CX6 and
# CX12 in the changepoint feature (all components with 0 values)
# Therefore, we want “fill the gaps” and ensure that each signature has weights on each feature.
# To do so, we used our previous methodology to generate the "signature interpretation matrix".


#### 6. DERIVATION OF SINGLE-CELL SIGNATURE INTERPRETATION MATRIX ##############
# In our previous work we established that shrinkage to zero weights of certain features
# is a common phenomenon when using NMF (Drews et al. Nature 2022). It also establishes
# that the missing weights are most certainly weak or moderate signals.
# NMF-inherent feature shrinkage explains why not all copy number signatures are defined
# across all 3 features (i.e. 0 values in all segsize components for CX1/6/12)


## Step 6.1: Prepare data
# Normalized SxC
SxC=readRDS(file=paste0(INTER_DIR, "/1.2.1_TCGA_SxC_100kb.rds"))
SxC[,c(1:20)] <- t(apply(SxC[,c(1:20)],1,function(x){x/sum(x)}))
SxC[,c(21:28)] <- t(apply(SxC[,c(21:28)],1,function(x){x/sum(x)}))
SxC[,c(29:37)] <- t(apply(SxC[,c(29:37)],1,function(x){x/sum(x)}))
COMPS <- colnames(SxC)

# Activities per sample
EXPS=readRDS(file=paste0(OUTPUT_DIR, "/2_New_scSignatures_Activities_THRES95.rds"))

# Definitions
DEFS=readRDS(file=paste0(OUTPUT_DIR,"/2_New_scCINSignatures_Definitions.rds"))


## Step 6.2: Linear correlations
# To check which component weights are correlated with signature activities
SIGS <- colnames(EXPS)

dfExps <- reshape2::melt(EXPS)
colnames(dfExps)<-c("sample", "sig", "act")
dfSxC <- reshape2::melt(t(SxC))
colnames(dfSxC)<-c("component", "sample", "vals")
dfLM <- left_join(dfExps,dfSxC, by="sample")

# Linear correlation
res <- dfLM[dfLM$act!=0,] %>% group_by(sig,component) %>% do(rstatix::tidy(lm(act ~ vals, data = .)))
res$beta <- res$term=="vals"
res <- res[res$beta==TRUE,] #remove intercept values
res$p.adj <- p.adjust(res$p.value, "BH")


## Step 6.3: Mix information of input matrix & activity matrix
SxC=readRDS(file=paste0(INTER_DIR, "/1.2.1_TCGA_SxC_100kb.rds"))
dfSxC <- reshape2::melt(t(SxC))
colnames(dfSxC)<-c("component", "sample", "vals")
dfStrength <- left_join(dfExps,dfSxC, by="sample")

# Multiply component weights and activity level per sample
dfStrength$strength <- dfStrength$vals*dfStrength$act

# Sum strength of each component per signature
dfStrength <- dfStrength %>% group_by(sig,component) %>% summarise(sum=sum(strength))
dfStrength <- left_join(dfStrength, res, by=c("sig", "component"))

# Normalize values to one
dfStrength[,"sum"] <- apply(dfStrength[,"sum"],2,function(x){x/sum(x)})


## Step 6.4: Add information of signature definitions for generating the interpretation matrix
dfDefs <- reshape2::melt(t(DEFS))
colnames(dfDefs)<-c("component", "sig", "weights")
dfInt <- left_join(dfStrength, dfDefs, by=c("sig", "component"))

# Normalize values to one
dfInt[,"weights"] <- apply(dfInt[,"weights"],2,function(x){x/sum(x)})

# Sum values from both matrices
dfInt$allSum <- dfInt$sum + dfInt$weights


## Step 6.5: Select positively correlated components
# Select components with significant (pvalue<0.05) and positive (beta>0) correlations
dfInt$selected <- dfInt$estimate>=0 & dfInt$p.adj<0.05
dfInt$allSum[dfInt$selected==FALSE]<-0

mInt<- as.matrix(do.call(rbind,lapply(SIGS, function(sig){
    s <- dfInt$allSum[dfInt$sig==sig]
    s <- s/sum(s)
    return(s)
})))
row.names(mInt)<-SIGS
colnames(mInt)<-colnames(DEFS)

pdf(file=paste0(PLOTS_DIR, "/1.2.6_scCINSignatures_Interpretation_Matrix.pdf"), width = 200/25.4, height = 150/25.4)
Heatmap(mInt, column_title = "scCINSigs Interpretations - 3 features with 37 components (20+8+9)", cluster_rows = FALSE, cluster_columns = FALSE)
dev.off()

saveRDS(mInt, file=paste0(OUTPUT_DIR,"/2_New_scCINSignatures_Interpretation.rds"))

## Step 6.6: Fill the gaps of the definition matrix
DEFS.CORRECTED=DEFS
DEFS.CORRECTED[1,21:28] <- mInt[1,21:28]
DEFS.CORRECTED[6,21:28] <- mInt[6,21:28]
DEFS.CORRECTED[12,21:28] <- mInt[12,21:28]

DEFS.CORRECTED=as.matrix(t(apply(DEFS.CORRECTED,1,function(x){x/sum(x)})))
saveRDS(DEFS.CORRECTED, file=paste0(OUTPUT_DIR,"/2_New_scCINSignatures_Definitions_Corrected.rds"))

# plot definitions
pdf(file=paste0(PLOTS_DIR, "/1.2.2_scCINSignatures_Definitions_Corrected_YAPSA.pdf"), width = 200/25.4, height = 150/25.4)
Heatmap(DEFS.CORRECTED, column_title = "scCINSigs Corrected Definitions - 3 features with 37 components (20+8+9)", cluster_rows = FALSE, cluster_columns = FALSE)
dev.off()


############ END SCRIPT ! ######################################################




