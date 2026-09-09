################################################################################
## Script for testing the new feature encoding using simulated data
## Started: 29/10/2022
## Updated: 28/07/2025
##
## This data was generated as part of our previous work (Drews et al 2022)
## We used this simulated samples for showing that our original feature encoding
## was able to capture different types of CIN. In this case, I'm going to use
## this simulated data to show that the substitution of 3 original features (bparm,
## bp10Mb, oscillating) by the new feature (bpwindow) is not affecting the ability
## to detect CIN types operating in a sample. 
################################################################################

rm(list=ls(all=TRUE))

## LIBRARIES
library(this.path)
library(data.table)
library(reshape2)
library(YAPSA)
library(ggplot2)
library(patchwork)
library(ggthemes)
library(lemon)
library(tidyverse)
library(foreach)
library(doMC)
library(igraph)
library(qgraph)
library(ComplexHeatmap)
library(ggplot2)
library(stringr)
library(Polychrome)

## PLOTS THEME
theme_set(theme_tufte(base_size = 6, base_family = "ArialMT"))
theme_update(text = element_text(size = 6),
             axis.text = element_text(size = 6),
             axis.line = element_line(size = 0.5), 
             axis.ticks = element_line(size = 0.5),
             axis.ticks.length = unit(.1, "cm"))


## PATHS
BASE=dirname(dirname(this.path()))
SIMDIR=dirname(this.path())
INPUT_DIR=file.path(BASE,"data")
# Create directories for outputs
OUTPUT_DIR=file.path(SIMDIR,"outputs")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
PLOTS_DIR=file.path(SIMDIR,"plots")
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(PLOTS_DIR,"nmf"), showWarnings = FALSE, recursive = TRUE)

### FUNCTIONS
source(file.path(BASE, "/scripts/newCNSigs_functions.R"))


## LOAD INPUT FILES => Download from figshare
lSEGTABLES=readRDS(paste0(INPUT_DIR, "/simulations/Out_2_Segmentation_tables_20each_N240.rds"))

# Signatures and Activities were derived with the SignatureDiscovery repo using the original approach
# The ideal comparison will be with definitions & exposures derived using the new 
# changepoint feature encoding. Since I have confirmed that the new encoding has not a significant
# impact on exposures, I'm just going to do the comparison with gold-standard matrices
# derived in Drews et al. 2022

# Signatures
SIMULATED.SIGNATURES=readRDS(paste0(INPUT_DIR, "/simulations/Signatures_XK5TRD_normalised.rds"))
SIGNAMES = c("W3"="HRDLST", "W5"="ecDNA", "W1"="CHR", "W4"="WGDearly", "W2"="WGDlate")

# Activities
SIMULATED.ACTIVITIES=readRDS(paste0(INPUT_DIR, "/simulations/Out_5_Activities_20each_N240.rds"))

# Overview of simulations
SIM.OVERVIEW=readRDS(paste0(INPUT_DIR, "/simulations/Out_2_Overview_simulation_20each_N240.rds"))

#################################################################################
### 0. PROCESS COPY NUMBER PROFILES #############################################
## Set normal segments 
SEGTABLES=data.table::rbindlist(lSEGTABLES)
SEGTABLES$segVal[SEGTABLES$segVal>(2-0.1) & SEGTABLES$segVal<(2+0.1)] = 2

# Smooth contiguous segments with equal segVal
df=idSmoothingTargets(SEGTABLES, WIGGLE = 0, colNameSegVal = "segVal", colNameChr = "chromosome", IGNOREDELS = FALSE)
lRaw=split(df, df$sample)
lSmooth=smoothSegments(lRaw, CORES = 1, SMOOTHINGFACTOR = 0, colNameMerge = "segVal", colNameChr = "chromosome",
                       colNameStart = "start", colNameEnd = "end", IGNOREDELS = FALSE, asDf = FALSE)


#################################################################################
### 1. EXTRACT COPY-NUMBER FEATURE ##############################################
lSEGTABLES=lSmooth
RMNORM=TRUE #remove normals
CORES=1
SIZE=5000000

SIM.ECNF=extractThreeCopynumberFeatures.new(lSEGTABLES, cores=CORES, 
                                            prePath=paste0(BASE,"/data"), 
                                            rmNorm=RMNORM, size=SIZE)

### 2. CALCULATE SxC MATRIX #####################################################
MODELS=readRDS(paste0(BASE,"/data/NEW_MIXMODELS.rds"))
FEATS=c("segsize","changepoint","bpwindow")
UNINFPRIOR=TRUE

lMats=lapply(FEATS, function(thisFeature) {
    print(thisFeature)
    thisEcnf=SIM.ECNF[[thisFeature]]
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
    id = which(rowSums(is.na(postDatScaled)) == ncol(postDatScaled))
    if(length(id)!=0){postDatScaled[id,]=c(rep(0,9),1)}
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

# Create input matrix
SxC=do.call(cbind, lMats)
saveRDS(SxC, file=paste0(INPUT_DIR, "/simulations/1.1_SxC_NewFeatures_20each_N240.rds"))


#### 3. EXTRACT SIGNATURES ######################################################

### Step 3.1: Apply BayesNMF
# I run the '01a_BrunetNMF_new_signature_derivation_simulations.R' script in the cluster
# I tested from K=3 to K=6, and I run the nmf algorithm 1000 times

### Step 3.2: Evaluate different solutions
lNMF=readRDS(paste0(OUTPUT_DIR,"/1.2_lExposures_allKs_20each_N240.rds"))
allKs=1:length(lNMF)

for(thisK in allKs){
    nmf=lNMF[[thisK]]
    
    
    ## Signature definitions
    defs=t(nmf[["defs"]])
    row.names(defs)=paste0("W",c(1:nrow(defs)))
    K=nrow(defs)
    
    pdf(file = paste0(PLOTS_DIR, "/nmf/1.3_CINSimulations_BrunetNMF_NewDefinitions_K", K,".pdf"), width = 200/25.4, height = 150/25.4)
    print(Heatmap(defs, column_title = "New CINSigs definitions derived from simulated profiles -- BrunetNMF",
                  cluster_rows = FALSE, cluster_columns = FALSE))
    dev.off()
    
    ## Exposures
    exps=t(nmf[["exps"]])
    colnames(exps)=paste0("W",c(1:ncol(exps)))
    
    # Compare activity levels with gold-standard approach
    correlation<-do.call(rbind,lapply(1:ncol(exps), function(thisSig){
        corr=cor(exps[,thisSig], SIMULATED.ACTIVITIES, method="pearson")}))
    row.names(correlation)<-colnames(exps)
    
    # Plot correlations
    pdf(file=paste0(SIMDIR, "/plots/nmf/1.3_PearsonCorrelations_Exposures_CINSimulations_K",K,".pdf"), width = 200/25.4, height = 150/25.4)
    print(Heatmap(correlation, column_title = "Pearson's r between exposures of original vs new signatures in simulated profiles", cluster_rows = FALSE, cluster_columns = FALSE,
                  cell_fun = function(j, i, x, y, width, height, fill) {
                      grid.text(sprintf("%.3f", correlation[i, j]), x, y, gp = gpar(fontsize = 6))
                  }))
    dev.off()
}

# Solution K=5 capture all 5 CIN types simulated! 

### 4. EVALUATE THE BEST SOLUTION K=5 ##########################################

## Evaluate presence of mutational process in each sample
thisK=5
WNAMES=c("WGDlate", "ecDNA", "WGDearly", "CHR", "HRDLST")

# Activities
EXPS=t(lNMF[[allKs[3]]]$exps)
mExps=scale(EXPS, center=TRUE, scale=TRUE)
colnames(mExps)=WNAMES

dfExps=reshape2::melt(mExps)
colnames(dfExps) = c("Samples", "Signature", "Activity")
dfExps$Samples = as.character(dfExps$Samples)

# Match samples order
thisOrder=SIM.OVERVIEW$Samples[order(SIM.OVERVIEW$Signature)]
dfExps$Samples=factor(dfExps$Samples, levels=rev(thisOrder))
dfExps$Simulation=SIM.OVERVIEW$Signature[match(as.character(dfExps$Samples),SIM.OVERVIEW$Samples)]

# Plot heatmap
p1 = ggplot(dfExps, aes(x = 1, y = Samples, fill = Simulation)) + geom_tile() + labs(x = "Simulation") +
    theme(axis.text = element_blank(), axis.ticks = element_blank()) + scale_fill_manual(values = unname(glasbey.colors(18)))
p2 = ggplot(dfExps, aes(x = Signature, y = Samples, fill = Activity)) + geom_tile() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.title.y = element_blank(), axis.line.y = element_blank()) + scale_fill_gradient2(low = "blue", mid = "white", midpoint = 0, high = "red")
pOut1 = p1 + p2 + plot_layout(nrow = 1, widths = c(0.05, 0.95), guides = "collect")

ggsave(paste0(SIMDIR, "/plots/1.4_Heatmap_Annotation_CINSimulations_K",thisK,".png"), 
       pOut1, width = 16, height = 16, units = "cm", bg = "white")


# Simplify names: remove for WGD the last entry - not dominant signature anyways
lGroup = lapply(unique(dfExps$Simulation), function(thisGroup) {
    
    # print(thisGroup)
    dfGroup = dfExps[dfExps$Simulation == thisGroup, ]
    
    # Trim down three element names because only two are dominant
    splitGroup = str_split(thisGroup, "_")[[1]]
    newGroupName = ifelse(length(splitGroup) == 3, paste0(splitGroup[1:2], collapse = "_"), thisGroup)
    
    # Replace Chr with CHR
    newGroupName = gsub("Chr", "CHR", newGroupName)
    
    dfGroup$ShortSim = newGroupName
    return(dfGroup)
})

dfShortGroup = do.call(rbind, lGroup)
dfShortGroup$Simulation = dfShortGroup$ShortSim
dfShortGroup$ShortSim = NULL

# Labeling top XX dominant samples for each each signature
topXX = c("HRDLST" = 100, "ecDNA" = 100, "CHR" = 100, "WGDearly" = 96, "WGDlate" = 24)

# Label top XX for each signature
allSigs = unique(dfShortGroup$Signature)
lSig = lapply(allSigs, function(thisSig) {
    
    dfSig=dfShortGroup[dfShortGroup$Signature == thisSig, ]
    dfSig$TopXX = "No"
    
    # Actually label samples
    numTop=topXX[names(topXX) == thisSig]
    dfSig$TopXX[order(dfSig$Activity, decreasing = TRUE)[1:numTop]] = "Yes"
    
    return(dfSig)
})
dfSig = do.call(rbind, lSig)


## Show of the top 100 samples how many have been predicted right.
# Correct HRDLST to LST
dfSig$Signature = factor(as.character(dfSig$Signature), levels = allSigs,
                         labels = c("WGDlate", "ecDNA", "WGDearly", "CHR", "LST"))

## Loop over each signature and get all samples with "Yes" or "No" (top XX samples)
allSigs = levels(dfSig$Signature)
lTab = lapply(allSigs, function(thisSig) {
    
    dfThisSig=dfSig[dfSig$Signature == thisSig, ]
    dfThisSig$RightGroup=ifelse(grepl(thisSig, dfThisSig$Simulation), "RightGroup", "WrongGroup")
    
    dfThisSig$RightGroup = factor(dfThisSig$RightGroup, levels = c("RightGroup", "WrongGroup"))
    dfThisSig$TopXX = factor(dfThisSig$TopXX, levels = c("No", "Yes"))
    
    # Call melt explicitly as reshape2 is now deprecated.
    dfSummary = reshape2::melt(table(dfThisSig$TopXX, dfThisSig$RightGroup))
    dfSummary$Signature = thisSig 
    
    return(dfSummary)
    
})
dfTab = do.call(rbind, lTab)
colnames(dfTab) = c("DominantSample", "Simulation", "N", "Signature")

dfRight = dfTab[dfTab$Simulation=="RightGroup",]
p3 = ggplot(dfRight, aes(x = Signature, y = N, fill = DominantSample)) + 
    geom_bar(stat = "identity", position = "fill") + 
    labs(y = "Proportion of correctly identified samples", fill = "Identified\nsamples with\ndominant signature?") +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = c("Yes" = "#1b9e77", "No" = "#d95f02" ))
p4 = ggplot(dfRight, aes(x = Signature, y = N, fill = DominantSample)) + 
    geom_bar(stat = "identity") + 
    labs(y = "Sample number", fill = "Identified\nsamples with\ndominant signature?") +
    scale_fill_manual(values = c("Yes" = "#1b9e77", "No" = "#d95f02" ))

## Calculate average accuracy
# Number of times a dominant signature was modelled
allDoms = sum(dfRight$N)

# Number of times we got it right
allRights = sum(dfRight$N[ dfRight$DominantSample == "Yes" ])
allRights/allDoms #70.47%

pOut3 = p3 + p4 + plot_layout(guides = "collect")

ggsave(paste0(SIMDIR, "/plots/1.4_Barplot_CorrectIndentification_CINSimulations_K",thisK,".svg"), 
       pOut3, width = 120, height = 60, units = "mm", bg = "white")


## Compare definitions with expected definitions
# Derived definitions
DEFS=t(lNMF[[allKs[3]]]$defs)
row.names(DEFS)=WNAMES

# Expected definitions
WGDlate=c(c(rep(0,4),rep(1,18)), #segsize
          c(rep(0,3),rep(1,3),rep(0,4)), #chp
          c(0,rep(1,5),rep(0,3))) #bpwindow

WGDearly=c(c(rep(0,15),rep(1,7)), #segsize
          c(rep(0,4),rep(1,2),rep(0,4)), #chp
          c(1,rep(0,8))) #bpwindow

CHR=c(c(rep(0,15),rep(1,7)), #segsize
      c(rep(0,2),rep(1,2),rep(0,6)), #chp
      c(1,rep(0,8))) #bpwindow

ecDNA=c(c(rep(0,4),rep(1,3),rep(0,15)), #segsize
        c(rep(0,8),rep(1,2)), #chp
        c(rep(0,2),rep(1,7))) #bpwindow

LST=c(c(rep(0,8),rep(1,7),rep(0,7)), #segsize
      c(rep(0,2),rep(1,2),rep(0,6)), #chp
      c(0,rep(1,5),rep(0,3))) #bpwindow

GOLD.DEFS=rbind(WGDlate,WGDearly,CHR,ecDNA,LST)
GOLD.DEFS=t(apply(GOLD.DEFS, 1, function(x) x/sum(x)))
    
# Cosine similarities
cosine<-do.call(rbind,lapply(1:nrow(DEFS), function(thisSig){
    cos=lsa::cosine(DEFS[thisSig,], t(GOLD.DEFS))}))
row.names(cosine)=paste0("SIM-",WNAMES)

# Plot definitions
scaleFUN <- function(x) sprintf("%.2f", x)

dfDefs = melt(t(DEFS))
colnames(dfDefs)<-c("component", "sig", "weight")

dfDefs$feature=gsub('[0-9]+', '', dfDefs$component)
dfDefs$feature=factor(dfDefs$feature, levels = c("segsize", "changepoint","bpwindow"))
dfDefs$sig=factor(dfDefs$sig, levels = c("WGDlate","WGDearly","CHR","ecDNA","HRDLST"))

p0 = ggplot(dfDefs, aes(x = component, y = weight, fill = feature))+
  geom_bar(stat = "identity") + facet_wrap(. ~ sig, ncol = 1, scales = "free_y") + 
  # coord_capped_cart(left = "both", bottom = "both") +
  labs(x = "Feature components", y = "Identified weight") +
  scale_y_continuous(labels=scaleFUN) +
  scale_fill_manual(values = c("#ffbe0b", "#fb5607", "#38ec83")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none") 

ggsave(paste0(SIMDIR, "/plots/1.4_Barplot_Definitions_CINSimulations_K",thisK,".svg"), 
       p0, width = 160, height = 100, units = "mm", bg = "white")


# Plot gold-standard definitions
dfGDefs = melt(t(GOLD.DEFS))
colnames(dfGDefs)<-c("component", "sig", "weight")
dfGDefs$component=dfDefs$component
dfGDefs$feature=dfDefs$feature
dfGDefs$feature=factor(dfGDefs$feature, levels = c("segsize", "changepoint","bpwindow"))

p0 = ggplot(dfGDefs, aes(x = component, y = weight, fill = feature))+
  geom_bar(stat = "identity") + facet_wrap(. ~ sig, ncol = 1, scales = "free_y") + 
  # coord_capped_cart(left = "both", bottom = "both") +
  labs(x = "Feature components", y = "Identified weight") +
  scale_y_continuous(labels=scaleFUN) +
  scale_fill_manual(values = c(alpha("#ffbe0b",0.4), alpha("#fb5607",0.4), alpha("#38ec83",0.4))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none") 

ggsave(paste0(SIMDIR, "/plots/1.4_Barplot_Definitions_CINSimulations_gold.svg"), 
       p0, width = 160, height = 100, units = "mm", bg = "white")

