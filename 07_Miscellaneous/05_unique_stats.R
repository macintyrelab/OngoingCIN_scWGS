################################################################################
## This script is for testing the confidence of unique events
## Started: 23/05/2023
## Updated: 09/04/2026
##
## For calling ongoingCIN we are going to use events that are unique to each cell.
## This is under the assumption that those genomic alterations shared across cells 
## are as result of a previous clonal expansion and therefore historical.
## 
## However, it's needed to be sure that those unique events are really unique. 
## Are these events really unique? Our approach compare raw absolute copy number
## values of bins covering the segment to the value in the rest of the cells. 
## Here I show some examples of false unique and true unique events
#################################################################################


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
MISC_DIR=dirname(this.path())
SCELL_DIR=file.path(BASE,"04_Process_scWGS")
IND_DIR=file.path(BASE,"05_Induced_CIN")
# Create directories for outputs
OUTPUT_DIR=file.path(MISC_DIR,"outputs")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
PLOTS_DIR=file.path(MISC_DIR,"plots")
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(PLOTS_DIR,"unique"), showWarnings = FALSE, recursive = TRUE)

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
# Unique events from in vitro BRCA1/2 models
dfUnique=readRDS(file.path(IND_DIR,"outputs/list_unique_per_model_BRCA_KO.rds"))
# Bins value
lSCELLS.BIN=readRDS(paste0(SCELL_DIR, "/intermediate/3.1.1_BRCA_KO_Bins.rds")) 
# Segments values
lSCELLS.CN=readRDS(paste0(SCELL_DIR, "/intermediate/3.1.1_BRCA_KO_Segments.rds"))
# Null distribution matrix
lDistributions=readRDS(paste0(SCELL_DIR, "/intermediate/confidence/NullDistribution_Confidence_Segsize.rds"))
dists=do.call(c,lDistributions)
quantile(dists)


################################################################################
#### PLOT CONFIDENCE VALUES OF UNIQUE EVENTS
t.test(dfUnique$minchange[dfUnique$significant=="yes"],dfUnique$minchange[dfUnique$significant=="no"])$p.value #7.588971e-33

#add null distribution
dtDists=as.data.frame(cbind(minchange=dists,significant="null"))
dtDists$minchange=as.numeric(dtDists$minchange)

dfPlot=rbind(dfUnique[,c(8,14)],dtDists)
dfPlot$significant=factor(dfPlot$significant,levels = c("null","no","yes"))

# Plotting distribution of segVals in all single clones
pdf(file=paste0(PLOTS_DIR, "/Distribution_confidence_unique_true_vs_false.pdf"), width = 100/25.4, height = 60/25.4)
p=ggplot(dfPlot) +
    geom_density(aes(x = minchange, color = significant, fill = significant)) +
    # facet_wrap(~segsize, scales = "free_y") +
    labs(title="Distribution of confidence values of unique events",subtitle=paste0("T-test p-value = 7.59e-33"))+
    scale_color_manual(values = alpha(c("grey40","blue","red"), 1))+ 
    scale_fill_manual(values = alpha(c("grey40","blue","red"), 0.5)) +
    scale_x_continuous(breaks=seq(0, 10, 1))
print(p)
dev.off()


#### PLOT SOME UNIQUE EVENTS ######

thisModel="scBRCA2_DLP22"
# for(thisModel in unique(false$model)){
    unique = dfUnique[dfUnique$model==thisModel,]
    bins = lSCELLS.BIN[[thisModel]]
    
    ## False unique events
    false=unique[unique$significant=="no",]
    
    # Overlaps with all cells
    unique.gr <- GRanges(seqnames = false$chr, IRanges(start = false$start, end = false$end))
    bins.gr = GRanges(seqnames = bins$chromosome, IRanges(start = bins$start, end = bins$end))
    overlaps <- findOverlaps(unique.gr, bins.gr)
    
    # Plot
    for (i in 1:nrow(false)){
        idxs=overlaps@to[overlaps@from == i]
        dat=bins[idxs,]
        dat$bin_id=rep(1:table(dat$sample)[1],length(unique(dat$sample)))
        sample=false$sample[i]
        othersamps=unique(dat$sample[dat$sample!=sample])
        
        # Plot
        xrange <- c(0,table(dat$sample)[1])
        yrange <- c(0,round(max(dat$binVal)+1))
        xat <- c(seq(0,table(dat$sample)[1],by=1))
        yat <- c(seq(0,round(max(dat$binVal)+1),by=1))

        pdf(paste0(PLOTS_DIR,"/unique/false_unique_",i,".pdf"),height = 80/25.4, width = 150/25.4)
        graphics::par(mfrow = c(1,1), xpd = FALSE, mar=c(4, 4, 4, 4), oma=c(0,1,0,2))
        plot(xrange, yrange, type="n", xaxt='n', yaxt='n', xlab="Number of bins", ylab="Copy number", cex.lab=1, bty='n')
        axis(side=1,at=xat,labels=xat,las=1,cex.axis=1)
        axis(side=2,at=yat,labels=yat,las=1,cex.axis=1)
        for(s in othersamps) {
            lines(x = dat$bin_id[dat$sample==s], y = dat$segVal[dat$sample==s], cex=1, lwd=2, col = "grey40")
            points(x = dat$bin_id[dat$sample==s], y = dat$binVal[dat$sample==s], pch=16, cex=1, col= alpha("grey40", 0.6)) 
        }
        lines(x = dat$bin_id[dat$sample==sample], y = dat$segVal[dat$sample==sample], cex=1, lwd=2, col = "blue")
        points(x = dat$bin_id[dat$sample==sample], y = dat$binVal[dat$sample==sample], pch=16, cex=1, col= alpha("blue", 0.6))
        dev.off()
        
        pdf(paste0(PLOTS_DIR,"/unique/false_unique_rounded_",i,".pdf"),height = 80/25.4, width = 150/25.4)
        graphics::par(mfrow = c(1,1), xpd = FALSE, mar=c(4, 4, 4, 4), oma=c(0,1,0,2))
        plot(xrange, yrange, type="n", xaxt='n', yaxt='n', xlab="Number of bins", ylab="Copy number", cex.lab=1, bty='n')
        axis(side=1,at=xat,labels=xat,las=1,cex.axis=1)
        axis(side=2,at=yat,labels=yat,las=1,cex.axis=1)
        for(s in othersamps) {
            lines(x = dat$bin_id[dat$sample==s], y = round(dat$segVal[dat$sample==s]), cex=1, lwd=2, col = "grey40")
            points(x = dat$bin_id[dat$sample==s], y = dat$binVal[dat$sample==s], pch=16, cex=1, col= alpha("grey40", 0.6)) 
        }
        lines(x = dat$bin_id[dat$sample==sample], y = round(dat$segVal[dat$sample==sample]), cex=1, lwd=2, col = "blue")
        points(x = dat$bin_id[dat$sample==sample], y = dat$binVal[dat$sample==sample], pch=16, cex=1, col= alpha("blue", 0.6))
        dev.off()
    }
    
    ## True unique events
    true=unique[unique$significant=="yes",]
    
    # Overlaps with all cells
    unique.gr <- GRanges(seqnames = true$chr, IRanges(start = true$start, end = true$end))
    bins.gr = GRanges(seqnames = bins$chromosome, IRanges(start = bins$start, end = bins$end))
    overlaps <- findOverlaps(unique.gr, bins.gr)
    
    # Plot
    for (i in 1:nrow(true)){
        idxs=overlaps@to[overlaps@from == i]
        dat=bins[idxs,]
        dat$bin_id=rep(1:table(dat$sample)[1],length(unique(dat$sample)))
        sample=true$sample[i]
        othersamps=unique(dat$sample[dat$sample!=sample])
        
        # Plot
        xrange <- c(0,table(dat$sample)[1])
        yrange <- c(0,round(max(dat$binVal)+1))
        xat <- c(seq(0,table(dat$sample)[1],by=1))
        yat <- c(seq(0,round(max(dat$binVal)+1),by=1))
        
        pdf(paste0(PLOTS_DIR,"/unique/true_unique_",i,".pdf"),height = 80/25.4, width = 150/25.4)
        graphics::par(mfrow = c(1,1), xpd = FALSE, mar=c(4, 4, 4, 4), oma=c(0,1,0,2))
        plot(xrange, yrange, type="n", xaxt='n', yaxt='n', xlab="Number of bins", ylab="Copy number", cex.lab=1, bty='n')
        axis(side=1,at=xat,labels=xat,las=1,cex.axis=1)
        axis(side=2,at=yat,labels=yat,las=1,cex.axis=1)
        for(s in othersamps) {
            lines(x = dat$bin_id[dat$sample==s], y = dat$segVal[dat$sample==s], cex=1, lwd=2, col = "grey40")
            points(x = dat$bin_id[dat$sample==s], y = dat$binVal[dat$sample==s], pch=16, cex=1, col= alpha("grey40", 0.6)) 
        }
        lines(x = dat$bin_id[dat$sample==sample], y = dat$segVal[dat$sample==sample], cex=1, lwd=2, col = "red")
        points(x = dat$bin_id[dat$sample==sample], y = dat$binVal[dat$sample==sample], pch=16, cex=1, col= alpha("red", 0.6))
        dev.off()
        
        pdf(paste0(PLOTS_DIR,"/unique/true_unique_rounded_",i,".pdf"),height = 80/25.4, width = 150/25.4)
        graphics::par(mfrow = c(1,1), xpd = FALSE, mar=c(4, 4, 4, 4), oma=c(0,1,0,2))
        plot(xrange, yrange, type="n", xaxt='n', yaxt='n', xlab="Number of bins", ylab="Copy number", cex.lab=1, bty='n')
        axis(side=1,at=xat,labels=xat,las=1,cex.axis=1)
        axis(side=2,at=yat,labels=yat,las=1,cex.axis=1)
        for(s in othersamps) {
            lines(x = dat$bin_id[dat$sample==s], y = round(dat$segVal[dat$sample==s]), cex=1, lwd=2, col = "grey40")
            points(x = dat$bin_id[dat$sample==s], y = dat$binVal[dat$sample==s], pch=16, cex=1, col= alpha("grey40", 0.6)) 
        }
        lines(x = dat$bin_id[dat$sample==sample], y = round(dat$segVal[dat$sample==sample]), cex=1, lwd=2, col = "red")
        points(x = dat$bin_id[dat$sample==sample], y = dat$binVal[dat$sample==sample], pch=16, cex=1, col= alpha("red", 0.6))
        dev.off()
    }
    
# }





