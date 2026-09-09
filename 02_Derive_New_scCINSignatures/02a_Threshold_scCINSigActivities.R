## Script to add noise to feature data and compare signature activities
# conda activate toposigs
# sbatch -J MCsims --mem=100G -t1440 -o log.txt -e error.txt --wrap "Rscript 02a_Threshold_scCINSigActivities.R"

rm(list=ls(all=TRUE))

## CHECK PACKAGES ARE INSTALLED 
listOfAllPackages = c("data.table","this.path","ggplot2","ggthemes",
                      "mclust","reshape2","devtools","BiocManager")
for(thisPackage in listOfAllPackages) {
    if(thisPackage %in% rownames(installed.packages()) == FALSE){
        print(paste0("Installing ", thisPackage))
        install.packages(thisPackage, repos="https://cran.ma.imperial.ac.uk/")
    }
}

thisPackage="patchwork"
if(thisPackage %in% rownames(installed.packages()) == FALSE){
    print(paste0("Installing ", thisPackage))
    devtools::install_github("thomasp85/patchwork")
}

thisPackage="YAPSA"
if(thisPackage %in% rownames(installed.packages()) == FALSE){
    print(paste0("Installing ", thisPackage))
    BiocManager::install("YAPSA")
}

## LOAD PACKAGES
suppressPackageStartupMessages(library(data.table))
suppressWarnings(library(YAPSA))
suppressWarnings(library(reshape2))
suppressWarnings(library(ggplot2))
suppressWarnings(library(ggthemes))
suppressWarnings(library(mclust))
suppressWarnings(library(this.path))
suppressWarnings(library(patchwork))

## Plot theme
theme_set(theme_light(base_size = 12, base_family = "Arial"))


## PATHS
BASE=dirname(dirname(this.path()))
SC_DIR=dirname(this.path())
OUTPATH=file.path(SC_DIR, "MCsimulations")
SAVEINTERMEDIATEFILES=TRUE
NAMEID="_scCINSigs_fullTCGA_1000sims_10pGaussian_10pSamplePoisson"

## Step 1
ECNF=file.path(SC_DIR, "intermediate/1.1.1_TCGA_Segments_100kb_ECNF.rds")
NUMSIMS=1e3
## Step 2
MODELS=file.path(BASE, "data/SC_MIXMODELS.rds")
UNINFPRIOR=TRUE
## Step 3
SIGNATUREFILE=file.path(SC_DIR, "outputs/2_New_scCINSignatures_Definitions.rds")
## Step 4
ACTIVITIES=file.path(SC_DIR, "outputs/2_New_scCINSignatures_Activities.rds")

################################################################################
## FUNCTIONS

### Step 1: Create credible distributions of features with noise
# 1. Create simulated values
addNoiseToFeatures = function(lOriECNF, allFeatures = c("segsize", "changepoint", "bpwindow"), 
                              SDPROP = 20, FINALWIDTH = 0.1, RANGECNAS = 0.1, NUMSIMS = 1) {
  
  lSim = lapply(1:NUMSIMS, function(thisN) {
    if(thisN %% 10 == 0) {print(thisN)}
    lFeatures = sapply(allFeatures, function(thisFeat) {
      dfFeat = lOriECNF[[thisFeat]]
      
      if(thisFeat %in% c("segsize")) {
        # Pass through option
        if(FINALWIDTH != 0) {
          dfFeat = addGaussianNoise(dfFeat, sdProp = SDPROP, finalWidth = FINALWIDTH)  
        }
      } else { 
        # Pass through option
        if(RANGECNAS != 0) {
          dfFeat = addSampleNoise(dfFeat, RANGECNAS = RANGECNAS)
        }
      } 
      
      return(dfFeat)
      
    }, simplify = FALSE, USE.NAMES = TRUE)
    
  })
  return(lSim)
}

# 2. Segment size => Gaussian noise
addGaussianNoise = function(dfFeat, sdProp = 20, finalWidth = NA) {
  
  ## This function adds noise to a column called "value".
  ## The level of noise is determined by the variable "sdProp" - sd = value / sdProp
  ## If wished the width of noise can be limited by "finalWidth" to a determined maximum width around
  ## the original value (to avoid large outliers).
  dfFeat$value2 = rnorm(n = length(dfFeat$value), mean = dfFeat$value, sd = dfFeat$value/sdProp)
  
  # Remove outliers if wished
  if(! is.na(finalWidth)) {
    
    dfFeat$diff = dfFeat$value / dfFeat$value2
    dfFeat$redraw = abs(1-dfFeat$diff) > finalWidth
    dfFeat$redraw[is.na(dfFeat$redraw)] = FALSE #omit artefacts => diff=NA when value2=0
    while(sum(dfFeat$redraw, na.rm=TRUE) > 0) {
      
      dfFeat$value2[dfFeat$redraw] = rnorm(n = length(dfFeat$value[dfFeat$redraw]), 
                                           mean = dfFeat$value[dfFeat$redraw], 
                                           sd = dfFeat$value[dfFeat$redraw]/sdProp)
      dfFeat$diff = dfFeat$value / dfFeat$value2
      dfFeat$redraw = abs(1-dfFeat$diff) > finalWidth
      dfFeat$redraw[is.na(dfFeat$redraw)] = FALSE #omit artefacts => diff=NA when value2=0
    }
  }
  
  # Clean up and return
  dfFeat$diff = NULL
  dfFeat$redraw = NULL
  dfFeat$value = dfFeat$value2
  dfFeat$value2 = NULL
  
  return(dfFeat)  
}

# 3. Changepoint & BPwindow => add noise
addSampleNoise = function(dfFeat, RANGECNAS = 0.1) {
  
  # Split into sample-specific lists for easier looping
  dfFeat$ID = factor(dfFeat$ID, levels = unique(dfFeat$ID))
  lFeat = split(dfFeat, dfFeat$ID)
  lSim = lapply(lFeat, function(thisSample) {
    
    numCNAs = sum(thisSample$value)
    tabSample = table(thisSample$value)
    
    # Shortcut for osCN feature where vectors of only zeros is possible
    if(numCNAs == 0) { return(thisSample) }
    
    # Extreme case if only one non-zero entry is present (apparently sample function fails)
    # Implement rescue to original values if only one value is present (independent of how many)
    if(length(tabSample) == 1) { return(thisSample) }
    
    # Draw new values with probabilities of existing sample
    vSim = sample(x = as.numeric(names(tabSample)), size = sum(tabSample), 
                  prob = tabSample/sum(tabSample), replace = TRUE)
    currentDraw = sum(vSim)
    
    # Check if this draw is above the user-defined limit. If so, enter while loop and draw until resolved
    while( abs(1-(currentDraw/numCNAs)) > RANGECNAS) {
      vSim = sample(x = as.numeric(names(tabSample)), size = sum(tabSample), 
                    prob = tabSample/sum(tabSample), replace = TRUE)
      currentDraw = sum(vSim)
    }
    
    thisSample$value = vSim
    return(thisSample)
  })
  
  # Merge list and prepare for return
  dfSim = do.call(rbind, lSim)
  rownames(dfSim) = NULL
  dfSim$ID = as.character(dfSim$ID)
  
  return(dfSim)
}

### Step 2: Derive SxC matrix
deriveSxCMatrices = function(lSimulation, allModels = NULL, allFeatures = names(allModels), UNINFPRIOR = TRUE) {
  
  if(is.null(allModels)) {
    stop("No models supplied.")
  }
  
  lMatrices = lapply(lSimulation, function(thisSimulation) {
    
    lMats = lapply(allFeatures, function(thisFeature) {
      
      thisEcnf = thisSimulation[[ thisFeature ]]
      thisModel = allModels[[ thisFeature ]]
      
      dat = as.numeric(thisEcnf[,2])
      # We want a posterior, hence likelihood (density) times prior (weight)
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
      # Normalise densities to probabilities
      postDatScaled = data.frame(postDatUnscaled/rowSums(postDatUnscaled))
      # Make sure that focal amplifications with huge number of copies have weight in the chp=8 
      id = which(rowSums(is.na(postDatScaled)) == ncol(postDatScaled))
      if (thisFeature == "changepoint"){
        if(length(id)!=0){postDatScaled[id,]=matrix(rep(c(rep(0,7),1), length(id)), ncol = 8, nrow =length(id), byrow = T)}
      }
      # Make sure that windows with a lot of segments have weight in the bpw=8
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
    
    allMats = do.call(cbind, lMats)
    return(allMats)
  })
}

### Step 3: Calculate signature activities
calculateSignatureActivities = function(lMatrices, W = NULL) {
  
  if(is.null(W)) {
    stop("No signatures supplied.")
  }
  
  lSignatures = lapply(lMatrices, function(V) {
    
    # Sanity check signature matrix
    if(nrow(W) < ncol(W)) W = t(W)
    # Sanity check input matrix
    if(nrow(V) > ncol(V)) V = t(V)
    
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
    Hraw = YAPSA::LCD(in_mutation_catalogue_df = V, 
                      in_signatures_df = YAPSA:::normalize_df_per_dim(W,2),
                      in_per_sample_cutoff = 0)
    H = as.matrix(t(apply(Hraw,2,function(x){x/sum(x)})))
    return(H)
  })
  
  return(lSignatures)
}

### Step 4: Plot & save results
# For identifying thresholds
q25 = function(x){quantile(x, probs = 0.25)}
q75 = function(x){quantile(x, probs = 0.75)}
q95 = function(x){quantile(x, probs = 0.95)}

# Plot boxplot for each signature and sample
plotAllSigs = function(lSignatures, dtOri, minmaxNorm = FALSE, orderOri = FALSE, DOTSIZE = 0.25,
                       TITLE = "TCGA samples", FILEOUT = NULL) {
  
  # Check if output path exists
  if(is.null(FILEOUT)){
    warning("No output path for threshold text file supplied. Defaulting to './'!")
    FILEOUT = "./Thresholds.txt"
  }
  
  # Bring list of signature activities into one data table
  lMelt = list()
  allSigs = length(lSignatures)
  for(i in 1:allSigs) {
    dtSignatures = data.table(melt(lSignatures[[i]]))
    dtSignatures$iter = i
    lMelt[[length(lMelt)+1]] = dtSignatures
  }
  
  dtSigs = rbindlist(lMelt)
  
  # Plot boxplot for each sample and signature
  allSigs = levels(dtSigs$Var2)
  lPlots = lapply(allSigs, function(thisSig) {
    
    dtCS2 = dtSigs[dtSigs$Var2 == thisSig, ]
    
    # Prepare original values
    dtOriCS2 = dtOri[dtOri$Var2 == thisSig, ]
    
    # Sort samples by original signature value or by median simulation value
    if(orderOri) {
      newOrder = as.character(dtOriCS2$Var1[order(dtOriCS2$value, decreasing = TRUE)])
    } else {
      dtOrder = aggregate(value ~ Var1, dtCS2, median)
      newOrder = as.character(dtOrder$Var1[order(dtOrder$value, decreasing = TRUE)])
    }
    
    # Reorder factors
    dtOriCS2$Var1 = factor(as.character(dtOriCS2$Var1), levels = newOrder)
    dtCS2$Var1 = factor(as.character(dtCS2$Var1), levels = newOrder)
    
    # Reorder the names themselves (needed for identifying the correct threshold later)
    dtOriCS2 = dtOriCS2[match(newOrder, dtOri$Var1), ]
    
    ## Four methods for identifying thresholds
    # Method 1: First time an IQR crosses zero
    dfQ25 = aggregate(value ~ Var1, dtCS2, q25)
    # Order of samples should match the sorting from above aka the factor levels
    if(! identical(as.character(dfQ25$Var1), levels(dfQ25$Var1))) {dfQ25 = dfQ25[match(newOrder, dfQ25$Var1), ]}
    # Identify first index, corresponding sample and original value when simulation reached zero
    firstZeroIndex = rle(dfQ25$value <= 0)$lengths[1]+1
    thresh1 = dtOriCS2$value[ firstZeroIndex ]
    
    
    # Method 2: 95% quantile of all true zero samples
    allZeroSamples = as.character(dtOriCS2$Var1[dtOriCS2$value == 0])
    
    dfQ75 = aggregate(value ~ Var1, dtCS2[ dtCS2$Var1 %in% allZeroSamples, ], q75)
    thresh2 = q95(dfQ75$value)
    
    
    # Method 3 & 4: Use values from true zero samples to fit a Gaussian mixture model and identify
    # first samples that don't fit into the Gaussian. 
    # Threshold 3: <5%
    # Threshold 4: <1%
    dat = dtCS2$value[dtCS2$Var1 %in% as.character(dtOriCS2$Var1[dtOriCS2$value == 0])]
    gModel = Mclust(dat, modelNames = "V", G = 1, verbose = FALSE)
    modelMean = gModel$parameters$mean
    modelVar = gModel$parameters$variance$sigmasq
    thresh3 = qnorm(0.95, mean = modelMean, sd = sqrt(modelVar), lower.tail = TRUE)
    thresh4 = qnorm(0.99, mean = modelMean, sd = sqrt(modelVar), lower.tail = TRUE)
    
    # Plot IQR
    p1 = ggplot(dtCS2, aes(x = Var1, y = value)) + geom_boxplot(outlier.shape = NA, coef = 0) +
      geom_point(data = dtOriCS2, aes(x = Var1, y = value), colour = "red", size = DOTSIZE) +
      # geom_point(aes(y = ori), colour = "red", size = DOTSIZE) +
      geom_hline(yintercept = thresh1, size = 1, linetype = "dotted", colour = "red") +
      geom_vline(xintercept = firstZeroIndex, size = 1, linetype = "dotted", colour = "red") +
      geom_hline(yintercept = thresh2, size = 1, linetype = "dashed", colour = "blue") +
      geom_hline(yintercept = thresh3, size = 1, linetype = "dotdash", colour = "green") +
      geom_hline(yintercept = thresh4, size = 1, linetype = "dotdash", colour = "green") +
      annotate("text", x = 400, y = thresh1*1.5, size = 6, label = as.character(signif(thresh1, digits = 2))) +
      annotate("text", x = 400, y = 0, size = 6, label = as.character(signif(thresh2, digits = 2))) +
      annotate("text", x = 6000, y = 0, size = 6, label = as.character(signif(thresh3, digits = 2))) +
      annotate("text", x = 6000, y = thresh4*1.5, size = 6, label = as.character(signif(thresh4, digits = 2))) +
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank()) +
      ggtitle(TITLE) + ylab(paste(thisSig)) +
      scale_x_discrete(breaks = NULL) + scale_y_continuous(minor_breaks = seq(0, 1, 0.05))
    
    return(list(p1, c(thresh1, thresh2, thresh3, thresh4)))
  })
  
  names(lPlots) = allSigs
  
  # Split plots and thresholds
  lOutPlot = list()
  lOutThresh = list()
  for(thisSig in allSigs) {
    lOutPlot[[length(lOutPlot)+1]] = lPlots[[thisSig]][[1]]
    lOutThresh[[length(lOutThresh)+1]] = c(thisSig, signif(lPlots[[thisSig]][[2]], digits = 4))
  } 
  
  # Combine thresholds and save output file
  dtOutThresh = data.table(do.call(rbind, lOutThresh))
  colnames(dtOutThresh) = c("Sig", "Thresh_Q25ZeroHit", "Thresh_Q95ofZeroQ75", 
                            "Thresh_ZeroGMM_0.05", "Thresh_ZeroGMM_0.01")
  write.table(dtOutThresh, FILEOUT, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  return(lOutPlot)
}


################################################################################
## Step 0: LOAD DATA
dir.create(OUTPATH, recursive = TRUE, showWarnings = FALSE)
lOriECNF = readRDS(ECNF)
allModels = readRDS(MODELS)
W = readRDS(SIGNATUREFILE)
dtOri = data.table(melt(readRDS(ACTIVITIES)))


## Step 1: Simulate feature distributions with noise
print("Start simulating data...")
FEATS=c("changepoint", "segsize","bpwindow")
lSimulation = addNoiseToFeatures(lOriECNF, allFeatures = FEATS, 
                                 SDPROP = 20, FINALWIDTH = 0.1, RANGECNAS = 0.1, NUMSIMS = NUMSIMS)


## Step 2: Derive SxC matrices
print("Start deriving SxC matrices data...")
lMatrices = deriveSxCMatrices(lSimulation, allModels = allModels, 
                              allFeatures = names(allModels), UNINFPRIOR = TRUE)
if(SAVEINTERMEDIATEFILES){
  saveRDS(lMatrices, file.path(OUTPATH, paste0("MC1_SxC_Matrices", NAMEID, ".rds")))
}


## Step 3: Calculate signature activities
print("Start calculating signature activities...")
lSignatures = calculateSignatureActivities(lMatrices, W = W)
if(SAVEINTERMEDIATEFILES) {
  saveRDS(lSignatures, file.path(OUTPATH, paste0("MC2_Activities", NAMEID, ".rds")))
}


# Step 4: Visualise results
# Plot boxplot for each signature and sample
lPlots = plotAllSigs(lSignatures, dtOri, orderOri = TRUE, DOTSIZE = 0.25, TITLE = "1000 runs full TCGA", 
                     FILEOUT = file.path(OUTPATH, paste0("MC3_Boxplots_per_sig", NAMEID, ".txt")))

for(i in 1:length(lPlots)){
  png(file.path(OUTPATH, paste0("MC3_Boxplots_per_sig", NAMEID,"_CX", i,".png")), width = 200, height = 160, units = "mm", res = 300) 
  print(lPlots[[i]])
  dev.off()
}

