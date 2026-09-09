################################################################################
############ Chapter 3: Derive signatures using the new approach ###############
### Call signatures using BrunetNMF
# sbatch -J bNMF -c 5 --mem=100G -t1440 -o log.txt -e error.txt --wrap "Rscript 01a_BrunetNMF_new_signature_derivation.R"


### LIBRARIES
rm(list=ls(all=TRUE))

## CHECK PACKAGES ARE INSTALLED 
# Install Biobase to be able to install R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if ("Biobase" %in% rownames(installed.packages()) == FALSE){
    print("Installing Biobase")
    BiocManager::install(version = "3.16")
    BiocManager::install("Biobase")
}

# Install packages
listOfAllPackages = c("NMF", "data.table","doMC","foreach", "argparse", 
                      "this.path", "ggplot2", "RColorBrewer")
for(thisPackage in listOfAllPackages) {
    if(thisPackage %in% rownames(installed.packages()) == FALSE){
        print(paste0("Installing ", thisPackage))
        install.packages(thisPackage, repos="https://cran.ma.imperial.ac.uk/", dependencies=TRUE)
    }
}


## LOAD PACKAGES
suppressPackageStartupMessages(library(data.table))
suppressWarnings(library(doMC))
suppressWarnings(library(foreach))
suppressWarnings(library(this.path))
suppressWarnings(library(ggplot2))
suppressWarnings(library(RColorBrewer))
suppressWarnings(library(NMF))
registerDoMC(5)
options(scipen = 999)

### PATHS
BASE=dirname(dirname(this.path()))
SIMDIR=dirname(this.path())
INPUT_DIR=file.path(BASE,"data")
# Create directories for outputs
OUTPUT_DIR=file.path(SIMDIR,"outputs")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
PLOTS_DIR=file.path(SIMDIR,"plots")
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)

################################################################################
####### Apply BrunetNMF to derive signatures using the new approach ############
SxC=readRDS(file=paste0(INPUT_DIR, "/simulations/1.1_SxC_NewFeatures_20each_N240.rds"))

### Step 1: Apply BrunetNMF
# Test different K
allKs=c(3,4,5,6)
NMFALG="brunet"
SEED=77777
NRUN=1000
NCORES=5

print("---- Starting Step 1: Apply BrunetNMF ----")
lNMF <- foreach(thisK=allKs)%dopar%{
    print(thisK)
    # Run the NMF.
    sigs=NMF::nmf(t(SxC), rank=thisK, seed=SEED, nrun=NRUN, method=NMFALG, .opt=paste0("p",NCORES))
    
    # Get W matrix => definitions
    sigmat=basis(sigs)
    
    # Get H matrix => activities
    expmat=coef(sigs)
    
    # Normalize matrices
    sigmat=apply(sigmat,2,function(x){x/sum(x)})
    expmat=apply(expmat,2,function(x){x/sum(x)})
    
    out=list()
    out[["defs"]]=sigmat
    out[["exps"]]=expmat
    return(out)
}
names(lNMF)=allKs
saveRDS(lNMF,paste0(OUTPUT_DIR,"/1.2_lExposures_allKs_20each_N240.rds"))

### Step 2: Decide optimal k solution
# Apart from manually check the definitions of all solutions, I also run this analysis 
# to select the best solution

print("---- Starting Step 2: Check different k solutions ----")
estim.r=NMF::nmfEstimateRank(t(SxC), allKs[1]:allKs[length(allKs)], seed=SEED, nrun=NRUN,
                             verbose=F,method=NMFALG, .opt=paste0("p",NCORES))
estim.r.random=NMF::nmfEstimateRank(randomize(t(SxC)), allKs[1]:allKs[length(allKs)], 
                                    seed=SEED, nrun=NRUN,
                                    verbose=F,method=NMFALG, .opt=paste0("p",NCORES))
p=plot(estim.r,estim.r.random, 
        what = c("cophenetic", "dispersion","sparseness", "silhouette"),xname="Observed",yname="Randomised",main="")+
    theme(axis.text=element_text(size=5),axis.title=element_text(size=5),
          strip.text.x = element_text(size = 5),
          strip.text.y = element_text(size = 5),
          legend.text = element_text(size = 5),
          legend.title = element_text(size = 7))
ggsave(paste0(PLOTS_DIR,"/1.2_Plot_simulations_statistics_brunetNMF.png"), p, type="cairo", 
       width=6, height=4, units="in", 
       dpi=300)
