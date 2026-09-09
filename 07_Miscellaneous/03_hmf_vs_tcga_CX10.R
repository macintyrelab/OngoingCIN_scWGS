################################################################################
## This script is to explore CX10 activities in TCGA and HMF cohorts
## Started: 27/10/25
## Updated: 23/04/26
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


### PATHS
BASE=dirname(dirname(this.path()))
MISC_DIR=dirname(this.path())
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



### LOAD INPUT DATA
# Signature activities
TCGA.EXPS=readRDS(file.path(BASE, 'data/tcga/Activities_NewCINSigs_TCGA_yapsa_k17_THRESH95.rds'))
HMF.EXPS=readRDS(file.path(BASE, 'data/hmf/Activities_NewCINSigs_HMF_yapsa_k17_THRESH95.rds')) # This requires data access permission from HMF

# Metadata
HMF.META=fread(file.path(BASE, 'data/hmf/hmf_metadata.tsv')) # This requires data access permission from HMF
HMF.META$patientIdentifier <- gsub('TI*V*$', '', HMF.META$sampleId, perl=TRUE)
HMF.TNBC=HMF.META[HMF.META$primaryTumorSubType=='Triple negative',]

TCGA.META=readRDS(file.path(BASE,'data/tcga/TCGA_clinical_data.RDS'))
TCGA.TNBC=TCGA.META[TCGA.META$type=="TNBC" & !is.na(TCGA.META$type),]

# Copy number profiles
HMF.CN.TNBC=readRDS(file.path(BASE, 'data/hmf/hmf_tnbc_segtable.rds')) # This requires data access permission from HMF
TCGA.CN.TNBC=readRDS(file.path(BASE, 'data/tcga/tcga_tnbc_segtable.rds'))


################################################################################
#### EXPLORE SIGNATURE CX10 ACTIVITIES IN BOTH COHORTS #########################

dtExps=rbind(TCGA.EXPS[row.names(TCGA.EXPS)%in%TCGA.TNBC$bcr_patient_barcode,],
             HMF.EXPS[row.names(HMF.EXPS)%in%HMF.TNBC$sampleId,])

# Scale activitites to be comparable
scExps=as.data.frame(scale(dtExps))

# Add cohort information
scExps$sample=row.names(scExps)
scExps$cohort="TCGA"
scExps$cohort[scExps$sample%in%HMF.TNBC$sampleId]="HMF"
table(scExps$cohort)

scExps$cohort=factor(scExps$cohort,levels = c("TCGA","HMF"))
effsize::cohen.d(scExps$CX10[scExps$cohort=="HMF"], scExps$CX10[scExps$cohort=="TCGA"], paired = FALSE) #cohens'd = 0.3094811

# Compare CX10 activities
pdf(file=paste0(PLOTS_DIR, "/BoxPlot_scCX10_TCGA_vs_HMF_TNBC.pdf"), width = 45/25.4, height = 40/25.4)
p=ggplot(scExps, aes(x = cohort, y = CX10)) + 
    geom_boxplot(outlier.shape = NA) +
    stat_compare_means(size=1.5)+
    geom_point(pch = 21, color = "grey20", fill = "grey20", position = position_jitter(0.2), alpha = 0.8, size = 1) +
    xlab("") + ylab("Scaled CX10 activities") +
    labs(title="",
         subtitle = "TCGA-TNBC vs HMF-TNBC")+
    theme(text = element_text(size = 6),
          axis.text = element_text(size = 6),
          axis.line = element_line(size = 0.5), 
          axis.ticks = element_line(size = 0.5),
          axis.ticks.length = unit(.1, "cm"))
print(p)
dev.off()


## Is CX10 correlated with number of CNAs or ploidy status?

# nCNAs
dCIN=rbind(as.data.frame(table(HMF.CN.TNBC$sample[HMF.CN.TNBC$segVal!=2])),
           as.data.frame(table(TCGA.CN.TNBC$sample[TCGA.CN.TNBC$segVal!=2])))
colnames(dCIN)=c("sample","nCNAs")
dCIN=left_join(dCIN,scExps[,c(10,18:19)],by="sample")
dCIN=dCIN[!is.na(dCIN$cohort),]

pdf(paste0(PLOTS_DIR, "/ScatterPlot_CX10_vs_nCNAs.pdf"), width = 80/25.4, height = 50/25.4)
p=ggplot(dCIN, aes(x=nCNAs, y=CX10)) + 
    geom_point(aes(color=cohort),size=1) +
    # facet_wrap(~cohort, ncol = 2, scales = "free_x") +
    geom_smooth(method = lm, formula = y ~ x, se = TRUE, col = "darkred") +
    ggpubr::stat_cor(method = "kendall", label.x = 0, label.y = 3.5, size = 2) +
    labs(title="") +
    xlab("Number of CNAs") +
    ylab("Scaled CX10 activities") +
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 5),
          axis.text.y = element_text(size = 5))
print(p)
dev.off()

# Ploidy
HMF.ploidy=as.data.frame(do.call(rbind,lapply(unique(HMF.CN.TNBC$sample),function(thisSamp){
    segTab=HMF.CN.TNBC[HMF.CN.TNBC$sample==thisSamp,]
    segLen<-(as.numeric(segTab$end)-as.numeric(segTab$start))
    ploidy<-sum((segLen/sum(segLen))*as.numeric(segTab$segVal))
    out<-c(thisSamp,ploidy)
    return(out)
})))
colnames(HMF.ploidy)=c("sample","ploidy")

TCGA.ploidy=as.data.frame(do.call(rbind,lapply(as.character(unique(TCGA.CN.TNBC$sample)),function(thisSamp){
    segTab=TCGA.CN.TNBC[TCGA.CN.TNBC$sample==thisSamp,]
    segLen<-(as.numeric(segTab$end)-as.numeric(segTab$start))
    ploidy<-sum((segLen/sum(segLen))*as.numeric(segTab$segVal))
    out<-c(thisSamp,ploidy)
    return(out)
})))
colnames(TCGA.ploidy)=c("sample","ploidy")

dfPloidy=rbind(HMF.ploidy,TCGA.ploidy)
dfPloidy=left_join(dfPloidy,scExps[,c(10,18:19)],by="sample")
dfPloidy=dfPloidy[!is.na(dfPloidy$cohort),]
dfPloidy$ploidy=as.numeric(dfPloidy$ploidy)

pdf(paste0(PLOTS_DIR, "/ScatterPlot_CX10_vs_Ploidy.pdf"), width = 80/25.4, height = 50/25.4)
p=ggplot(dfPloidy, aes(x=ploidy, y=CX10)) + 
    geom_point(aes(color=cohort),size=1) +
    # facet_wrap(~cohort, ncol = 2, scales = "free_x") +
    geom_smooth(method = lm, formula = y ~ x, se = TRUE, col = "darkred") +
    ggpubr::stat_cor(method = "kendall", label.x = 2, label.y = 3.5, size = 2) +
    labs(title="") +
    xlab("Ploidy") +
    ylab("Scaled CX10 activities") +
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 5),
          axis.text.y = element_text(size = 5))
print(p)
dev.off()



######## END OF THE SCRIPT! ####################################################

