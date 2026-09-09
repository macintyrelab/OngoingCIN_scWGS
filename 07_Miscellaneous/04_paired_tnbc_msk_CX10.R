################################################################################
## This script is to explore CX10 activities in paired data from the MSK cohort
## Started: 07/11/25
## Updated: 24/03/26
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



### LOAD INPUT DATA
# Signature activities --> FACETS-derived profiles from MSK-IMPACT data & signature quantification
# Data access needs to be requested
EXPS=readRDS(file.path(BASE, 'data/mskcc/tnbc_msk_cx10.rds'))
EXPS$scCX10=as.numeric(scale(EXPS$CX10))
EXPS$sample_type=factor(EXPS$sample_type, levels = c("Primary","Metastasis"))

# Other genomic metrics
MSK.metrics=fread(file.path(BASE, 'data/mskcc/TNBC_CX10_ploidy_CNA.tsv'))


################################################################################
#### EXPLORE SIGNATURE CX10 ACTIVITIES IN PAIRED DATA ##########################

## Prepare data
Prims=EXPS[EXPS$sample_type=="Primary",]
colnames(Prims)[1]="CX10_prims"
colnames(Prims)[4]="scCX10_prims"
Mets=EXPS[EXPS$sample_type=="Metastasis",]
colnames(Mets)[1]="CX10_mets"
colnames(Mets)[4]="scCX10_mets"

dfExps=left_join(Mets[,c(1,3,4)],Prims[,c(1,3,4)],by="patient")
dfExps$absence=dfExps$CX10_mets==0&dfExps$CX10_prims==0
# res=wilcox.test(dfExps$scCX10_prims, dfExps$scCX10_mets, paired = TRUE, alternative = "less")
# res$p.value #0.02582311
res=wilcox.test(dfExps$scCX10_prims[dfExps$absence==F], dfExps$scCX10_mets[dfExps$absence==F], paired = TRUE, alternative = "less")
res$p.value #0.02471924
effsize::cohen.d(dfExps$scCX10_mets[dfExps$absence==F], dfExps$scCX10_prims[dfExps$absence==F], paired = TRUE) #cohens'd = 1.010515

# Compare CX10 activities between primaries and mets
pdf(file=paste0(PLOTS_DIR, "/BoxPlot_scCX10_MSK_TNBC.pdf"), width = 45/25.4, height = 40/25.4)
p=ggplot(EXPS[EXPS$patient%in%dfExps$patient[dfExps$absence==F],], aes(x = sample_type, y = scCX10)) + 
    geom_boxplot(outlier.shape = NA) +
    geom_line(aes(group=patient), position = position_dodge(0.2), alpha = 0.8, size = 0.4) +
    geom_point(aes(fill=sample_type,group=patient), pch = 21, color = "grey20", position = position_dodge(0.2), alpha = 0.8, size = 1) +
    xlab("") + ylab("Scaled CX10 activities") +
    labs(title="",
         subtitle = "Paired wilcoxon test p=0.025")+
    theme(text = element_text(size = 6),
          axis.text = element_text(size = 6),
          axis.line = element_line(size = 0.5), 
          axis.ticks = element_line(size = 0.5),
          axis.ticks.length = unit(.1, "cm"),
          legend.position = "none")
print(p)
dev.off()


## Is CX10 correlated with number of CNAs or ploidy status?
MSK.metrics$increase=dfExps$increase[match(MSK.metrics$patient,dfExps$patient)]
MSK.metrics$sample_type=factor(MSK.metrics$sample_type,levels=c("Primary","Metastasis"))
MSK.metrics=MSK.metrics[!duplicated(MSK.metrics[,1:2]),]

pdf(paste0(PLOTS_DIR, "/BoxPlots_msk_CX10_vs_nCNAs_vs_ploidy.pdf"), width = 90/25.4, height = 50/25.4)
p1=ggplot(MSK.metrics[MSK.metrics$patient%in%dfExps$patient[dfExps$absence==F],], aes(x = sample_type, y = ploidy)) + 
    geom_boxplot(outlier.shape = NA) +
    geom_line(aes(group=patient), position = position_dodge(0.2), alpha = 0.8, size = 0.4) +
    geom_point(aes(fill=increase,group=patient), pch = 21, color = "grey20", position = position_dodge(0.2), alpha = 0.8, size = 1) +
    xlab("") + ylab("Ploidy") +
    labs(title="",
         subtitle = "")+
    theme(text = element_text(size = 6),
          axis.text = element_text(size = 6),
          axis.line = element_line(size = 0.5), 
          axis.ticks = element_line(size = 0.5),
          axis.ticks.length = unit(.1, "cm"),
          legend.position = "none")
p2=ggplot(MSK.metrics[MSK.metrics$patient%in%dfExps$patient[dfExps$absence==F],], aes(x = sample_type, y = n_CNA)) + 
    geom_boxplot(outlier.shape = NA) +
    geom_line(aes(group=patient), position = position_dodge(0.2), alpha = 0.8, size = 0.4) +
    geom_point(aes(fill=increase,group=patient), pch = 21, color = "grey20", position = position_dodge(0.2), alpha = 0.8, size = 1) +
    xlab("") + ylab("Number of CNAs") +
    labs(title="",
         subtitle = "")+
    theme(text = element_text(size = 6),
          axis.text = element_text(size = 6),
          axis.line = element_line(size = 0.5), 
          axis.ticks = element_line(size = 0.5),
          axis.ticks.length = unit(.1, "cm"),
          legend.position = "none")
p=ggarrange(p1,p2,ncol=2,hjust=1)
print(p)
dev.off()

MSK.metrics$scCX10_prims=dfExps$scCX10_prims[match(MSK.metrics$patient,dfExps$patient)]
MSK.metrics$scCX10_mets=dfExps$scCX10_mets[match(MSK.metrics$patient,dfExps$patient)]
mets=MSK.metrics[MSK.metrics$sample_type=="Metastasis",]
prims=MSK.metrics[MSK.metrics$sample_type!="Metastasis",]

pdf(paste0(PLOTS_DIR, "/ScatterPlots_msk_CX10_vs_nCNAs_vs_ploidy.pdf"), width = 100/25.4, height = 100/25.4)
p1=ggplot(prims[prims$patient%in%dfExps$patient[dfExps$absence==F],], aes(x=ploidy, y=scCX10_prims)) + 
    geom_point(aes(color = increase),size=1) +
    geom_smooth(method = lm, formula = y ~ x, se = TRUE, col = "darkred") +
    ggpubr::stat_cor(method = "kendall", label.x = 2, label.y = 2, size = 2) +
    labs(title="", subtitle="Primary") +
    xlab("Ploidy") +
    ylab("Scaled CX10 activities") +
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 5),
          axis.text.y = element_text(size = 5),
          legend.position = "none")
p2=ggplot(mets[mets$patient%in%dfExps$patient[dfExps$absence==F],], aes(x=ploidy, y=scCX10_mets)) + 
    geom_point(aes(color = increase),size=1) +
    geom_smooth(method = lm, formula = y ~ x, se = TRUE, col = "darkred") +
    ggpubr::stat_cor(method = "kendall", label.x = 2, label.y = 4, size = 2) +
    labs(title="", subtitle="Metastasis") +
    xlab("Ploidy") +
    ylab("Scaled CX10 activities") +
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 5),
          axis.text.y = element_text(size = 5),
          legend.position = "none")
p3=ggplot(prims[prims$patient%in%dfExps$patient[dfExps$absence==F],], aes(x=n_CNA, y=scCX10_prims)) + 
    geom_point(aes(color = increase),size=1) +
    geom_smooth(method = lm, formula = y ~ x, se = TRUE, col = "darkred") +
    ggpubr::stat_cor(method = "kendall", label.x = 20, label.y = 2, size = 2) +
    labs(title="", subtitle="Primary") +
    xlab("Number of CNAs") +
    ylab("Scaled CX10 activities") +
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 5),
          axis.text.y = element_text(size = 5),
          legend.position = "none")
p4=ggplot(mets[mets$patient%in%dfExps$patient[dfExps$absence==F],], aes(x=n_CNA, y=scCX10_mets)) + 
    geom_point(aes(color = increase),size=1) +
    geom_smooth(method = lm, formula = y ~ x, se = TRUE, col = "darkred") +
    ggpubr::stat_cor(method = "kendall", label.x = 20, label.y = 4, size = 2) +
    labs(title="", subtitle="Metastasis") +
    xlab("Number of CNAs") +
    ylab("Scaled CX10 activities") +
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 5),
          axis.text.y = element_text(size = 5),
          legend.position = "none")
p=ggarrange(p1,p2,p3,p4,ncol=2,nrow=2,hjust=T)
print(p)
dev.off()


############ END OF THE SCRIPT ! ###############################################
