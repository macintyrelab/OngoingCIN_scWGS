################################################################################
## This script is to explore number of cells for unique detection
## Started: 24/10/25
## Updated: 22/06/26
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


cores = 10
registerDoMC(cores = cores)


### PATHS
BASE=dirname(dirname(this.path()))
MISC_DIR=dirname(this.path())
SCELL_DIR=file.path(BASE,"04_Process_scWGS")
# Create directories for outputs
OUTPUT_DIR=file.path(MISC_DIR,"outputs")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
PLOTS_DIR=file.path(MISC_DIR,"plots")
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
source(file.path(BASE, "scripts/uniqueEvents_functions.R"))


### PARAMETERS
THRESH=0.90 #% for defining a shared event as recent
SIGN.THRESH=0.05 #for defining significant unique segments

REMOVEQUIET=TRUE 
WINDOW=100e3 
RMNORM=TRUE
CORES=1
WIGGLE=0.1
IGNOREDELS=FALSE
SIZE=5000000
UNINFPRIOR=TRUE
RESOLUTION="100kb"
BIN.SIZE=as.numeric(substr(RESOLUTION,1,nchar(RESOLUTION)-2))*1000


### LOAD INPUT DATA
# # Raw data
# sc.files=dir(paste0(BASE,"/data/minussi_sc"),full=T,pattern="absoluteCN")
# sc.files=c(readRDS(sc.files[4]),readRDS(sc.files[7]),readRDS(sc.files[8]))

# Profiles 
SCELLS.CN=readRDS(paste0(SCELL_DIR, "/intermediate/1.0_Minussi_Powered_Segments.rds")) 
SAMPS=unique(SCELLS.CN$model)

# Bin tables
BIN.TABLE=readRDS(paste0(SCELL_DIR, "/intermediate/1.1.1_Minussi_Bins.rds")) 
SEG.TABLE=readRDS(paste0(SCELL_DIR, "/intermediate/1.1.1_Minussi_Segments.rds")) 

# Feature components
MIXMODELS=readRDS(paste0(BASE,"/data/SC_MIXMODELS.rds"))

# Null distribution matrix => for unique calling
lDistributions=readRDS(paste0(SCELL_DIR, "/intermediate/confidence/NullDistribution_Confidence_Segsize.rds"))


################################################################################
#### UNIQUE EVENTS DETECTION (GROUND TRUE) #####################################

# I run the statistic test to call real unique events in these samples, and 
# also to get subconal events
# see script "01a_minussi_real_unique.rds"


################################################################################
#### SUBSAMPLING TEST ##########################################################

# We run in the HPC cluster a subsampling analysis where an increased number of cells
# is selected at each iteration. This is run 20 times
# see script "01b_unique_fp_detection.rds"

# We rerun in the HPC cluster a subsampling analysis to check how many cells of the 
# whole population carry the false calls. This is run 20 times
# see script "01c_unique_fp_cells_in_population.rds"


################################################################################
#### PLOTTING SUBSAMPLING TEST #################################################

## 1) Plotting false unique detected => Figure 1c
dfTests=rbind(readRDS(file.path(OUTPUT_DIR,paste0("stats_ncells_unique_fp_detection_TN1.rds"))),
              readRDS(file.path(OUTPUT_DIR,paste0("stats_ncells_unique_fp_detection_TN2.rds"))),
              readRDS(file.path(OUTPUT_DIR,paste0("stats_ncells_unique_fp_detection_TN5.rds"))),
              readRDS(file.path(OUTPUT_DIR,paste0("stats_ncells_unique_fp_detection_MDAMB231_popp31.rds"))),
              readRDS(file.path(OUTPUT_DIR,paste0("stats_ncells_unique_fp_detection_MDAMB231c8.rds"))))


# compute ratio of total cells
for(thisModel in SAMPS){
    ntotal=max(dfTests$ncells[dfTests$model==thisModel])
    dfTests$tcr[dfTests$model==thisModel]=dfTests$ncells[dfTests$model==thisModel]/ntotal
}

# compute ratio of fp
dfTests$fpr=dfTests$nfalse/dfTests$nunique
dfTests$fpr[is.na(dfTests$fpr)]=0



pdf(paste0(PLOTS_DIR, "/RateFalse_vs_RateCells.pdf"), width = 100/25.4, height = 60/25.4)
p = ggplot(dfTests, aes(x = tcr, y = fpr)) + 
    geom_point(aes(color = model, fill=model), size=0.3, alpha = 0.2) +
    geom_smooth(aes(color = model), method = "lm", size=1, fullrange = TRUE) +
    scale_x_continuous(breaks=seq(0,1,0.1)) +
    scale_y_continuous(breaks=seq(0,1,0.1)) +
    scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#931bd8","#1bd893")) +
    labs(title="",
         subtitle="20 iterations",
         x = "Fraction of cells subsampled", y = "FPR of unique events")
print(p)
dev.off()


## 2) Plotting fraction of cells carrying false unique => Added in SuppNote 1
dfTests=rbind(readRDS(file.path(OUTPUT_DIR,paste0("stats_false_unique_cells_in_pop_TN1.rds"))),
              readRDS(file.path(OUTPUT_DIR,paste0("stats_false_unique_cells_in_pop_TN2.rds"))),
              readRDS(file.path(OUTPUT_DIR,paste0("stats_false_unique_cells_in_pop_TN5.rds"))),
              readRDS(file.path(OUTPUT_DIR,paste0("stats_false_unique_cells_in_pop_MDAMB231_popp31.rds"))),
              readRDS(file.path(OUTPUT_DIR,paste0("stats_false_unique_cells_in_pop_MDAMB231c8.rds"))))

# Drop events that did not map to the final population (NA) if any
dfTests <- dfTests[!is.na(dfTests$n_cells_in_pop), ]

# Sampling fraction
dfTests$tcr <- dfTests$ncells / dfTests$ntotal

# Population fraction carrying the false unique event
dfTests$fpop <- dfTests$n_cells_in_pop / dfTests$ntotal


pdf(paste0(PLOTS_DIR, "/FracPopFalseUnique_vs_RateCells.pdf"), width = 100/25.4, height = 60/25.4)
p = ggplot(dfTests, aes(x = tcr, y = fpop)) +
    geom_point(aes(color = model, fill = model), size=0.3, alpha = 0.2) +
    geom_smooth(aes(color = model), method = "loess", size=1, fullrange = TRUE, se = FALSE) +
    ylim(0,1) +
    scale_x_continuous(breaks=seq(0,1,0.1)) +
    scale_y_continuous(breaks=seq(0,1,0.1)) +
    scale_color_manual(values=c("#D81B60","#FFC107","#1E88E5","#931bd8","#1bd893")) +
    labs(title="",
         subtitle="20 iterations",
         x = "Fraction of cells subsampled", y = "Fraction of cells\ncarrying false unique event")
print(p)
dev.off()

