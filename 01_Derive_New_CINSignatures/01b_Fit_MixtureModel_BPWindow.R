# ## Script for fitting mixture model for BPwindow feature

# To be sure this is properly loaded
library(this.path)
# Poisson
library(flexmix)
# Gaussian
library(mclust)
library(data.table)
library(ggplot2)
library(cowplot)
# Parallel 
library(foreach)
library(doMC)

## PATHS
BASE=dirname(dirname(this.path()))
NEWSIGS_DIR=dirname(this.path())

FEATURE="nsegments"
MINCOMP=1
MAXCOMP=10
FITDIST="pois"
PRIOR=0.001
NITER=1000
DECISION="BIC"
SEED=77777

# Functions
source(file.path(BASE, "scripts/newCNSigs_functions.R"))
source(file.path(BASE, "scripts/main_functions.R"))

# Load data
TCGA.ECNF=readRDS(paste0(NEWSIGS_DIR,"/intermediate/1.1_TCGA_ThreeECNF.rds"))
bpw=TCGA.ECNF[['bpwindow']]

################################################################################

# Fit model
dat=bpw$value

#### Fixed categories for dat<=6 (0 in 10%, 1 in 30%, 2 in 25%, 3 in 15%, 4/5/6 in 5%) ####
Lambda=c(0,1,2,3,4,5,6) #component capturing those segments without or with low closer breakpoints

#### For data>6 (5%) => fit model for the tail of the distribution ####
dt=dat[dat>6] #remove zeros

fittedDist = fitComponent2(dt, dist=FITDIST, min_prior=PRIOR, min_comp=MINCOMP, 
                           max_comp=MAXCOMP, niter=NITER, model_selection=DECISION, seed=SEED)

# Convert cluster components to a table for easier validation
if(FITDIST == "pois") {
    # Poisson distribution -> flexmix object
    listComp = unlist(fittedDist@components)
    allComps = sapply(listComp, function(thisComp) {unlist(thisComp@parameters)})
    dtComp = data.table(Lambda = sort(allComps))
} else {
    # Normal distribution -> Mclust object
    dtComp = data.table(Mean = fittedDist$parameters$mean,
                        SD = sqrt(fittedDist$parameters$variance$sigmasq))
    dtComp = dtComp[order(Mean),]
}
dtComp=rbind(as.data.frame(Lambda),dtComp)

################################################################################
# Save
saveRDS(dtComp, file = paste0(NEWSIGS_DIR, "/intermediate/MIXMODEL_BPWINDOW.rds"))

