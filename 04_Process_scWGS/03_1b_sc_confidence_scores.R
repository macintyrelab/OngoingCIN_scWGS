################################################################################
## This script is to compute confidence scores for each segment
## Run in HPC cluster as follow
## conda activate swgs-relcn
## sbatch -J scores -c10 -t1440 --mem=100G -o log.txt -e error.txt --wrap "Rscript 03_1b_sc_confidence_scores.R --cntable 3.1.1_BRCA_KO_Segments.rds --bintable 3.1.1_BRCA_KO_Bins.rds --sample 1"
## IMP! Modify path & files names prior to run the .sh script
################################################################################

### LIBRARIES
rm(list=ls(all=TRUE))

library(argparse)
library(this.path)
library(data.table)
library(QDNAseq)
library(GenomicRanges)
library(foreach)
library(doMC)

cores = 10
registerDoMC(cores = cores)

## ARGUMENTS
parser <- argparse::ArgumentParser(description='get confidence scores')
parser$add_argument('--cntable', type = "character", required = TRUE)
parser$add_argument('--bintable', type = "character", required = TRUE)
parser$add_argument('--sample', type = "numeric", required = TRUE)
args <- parser$parse_args()

### PATHS
BASE=dirname(this.path())
INPUT_DIR=file.path(BASE,"intermediate")
OUTPUT_DIR=file.path(BASE,"intermediate/confidence")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)


### LOAD INPUT DATA
# SEGMENTS
lcntables=args$cntable
lcntables=readRDS(file.path(INPUT_DIR, lcntables))

# BINS
bintable=args$bintable
lbintables=readRDS(file.path(INPUT_DIR, bintable))


################################################################################
#### Null distribution of confidence fittings ##################################

thisTab=args$sample

## Get confidence matrix per segment
print("Confidence values per segment...")
print(paste0("Sample ID: ",names(lcntables)[[thisTab]]))
dir.create(file.path(OUTPUT_DIR,names(lcntables)[[thisTab]]), showWarnings = FALSE, recursive = TRUE)

# Per cell
cntable=lcntables[[thisTab]]
bintable=lbintables[[thisTab]]
cells=unique(cntable$sample)

# confidence matrix
foreach(thisCell=cells) %dopar% {
    print(paste0("Cell ID: ",thisCell))
    
    tab=cntable[cntable$sample==thisCell,]
    
    # get confidence scores
    confis = as.data.frame(do.call(rbind,lapply(1:nrow(tab), function(thisSeg){
        print(paste0("Segment ID: ",thisSeg))
        # get genomic location of segment
        segs=tab[thisSeg,]
        segments.gr = GRanges(seqnames = segs$chromosome, IRanges(start = segs$start, end = segs$end))
        
        # get raw value in specific cell
        bin=bintable[bintable$sample%in%segs$sample & bintable$chromosome%in%segs$chromosome,]
        bins.gr=GRanges(seqnames = bin$chromosome, IRanges(start = bin$start, end = bin$end))
        overlaps <- findOverlaps(segments.gr, bins.gr)
        rawVal=mean(bin[overlaps@to,4])
        
        # get raw values in the other cells
        bin=bintable[!bintable$sample%in%segs$sample & bintable$chromosome%in%segs$chromosome,]
        bins.gr=GRanges(seqnames = bin$chromosome, IRanges(start = bin$start, end = bin$end))
        overlaps <- findOverlaps(segments.gr, bins.gr)
        oths=bin[overlaps@to,]
        oths <- do.call(c,lapply(unique(oths$sample), function(thisOth){
            bins=mean(oths[oths$sample==thisOth,4])
            return(bins)
        }))
        
        # confidence score => mode of difference with value in other cells
        confi=abs(rawVal-round(oths))
        min=min(confi)
        confi=as.numeric(names(sort(-table(confi))[1]))
        
        #output
        segs$confidence=confi
        segs$minchange=min
        return(segs)
    })))
    colnames(confis)=c("chromosome","start","end","segVal","sample","confidence","minchange")
    saveRDS(confis, paste0(file.path(OUTPUT_DIR,names(lcntables)[[thisTab]]),"/",thisCell,"_Segments_Confidence.rds"))
}

confidence.files=dir(file.path(OUTPUT_DIR,names(lcntables)[[thisTab]]),full=T,pattern="Confidence")
Confidence=as.data.frame(do.call(rbind,lapply(confidence.files,function(thisFile){
    segs=readRDS(thisFile)
    return(segs)
})))
saveRDS(Confidence, paste0(OUTPUT_DIR,"/",names(lcntables)[[thisTab]],"_Segments_Confidence.rds"))

print("....done!")



