library("NMF")
library("quadprog")
#library("ComplexHeatmap")
#library("circlize")
# for the manual parallelisation and also used by NMF for their internal parallelisation
library("foreach")
# needed for faster nmf runs based on the developer's own recommendations:
# https://cran.r-project.org/web/packages/NMF/vignettes/NMF-vignette.pdf
# page 16 ff.
library("doParallel")
# Libraries "bigmemory" and "synchronicity" are recommended for non-Windows machines
# for remodelling data structures
library("reshape2")
# for plotting
library("ggplot2")
#library("Cairo")
# Correlation of signatures
library("nnls")
library("corrplot")


locationOfThisScript = function() # Function LocationOfThisScript returns the location of this .R script (may be needed to source other files in same dir)
{
    this.file = NULL
    # This file may be 'sourced'
    for (i in -(1:sys.nframe())) {
        if (identical(sys.function(i), base::source)) this.file = (normalizePath(sys.frame(i)$ofile))
    }
    
    if (!is.null(this.file)) return(dirname(this.file))
    
    # But it may also be called from the command line
    cmd.args = commandArgs(trailingOnly = FALSE)
    cmd.args.trailing = commandArgs(trailingOnly = TRUE)
    cmd.args = cmd.args[seq.int(from=1, length.out=length(cmd.args) - length(cmd.args.trailing))]
    res = gsub("^(?:--file=(.*)|.*)$", "\\1", cmd.args)
    
    # If multiple --file arguments are given, R uses the last one
    res = tail(res[res != ""], 1)
    if (0 < length(res)) return(dirname(res))
    
    # Both are not the case. Maybe we are in an R GUI?
    return(NULL)
}
this_path<-locationOfThisScript()
source(paste(this_path,"helper_functions.R",sep="/"))


#### Extract features ####

extractCopynumberFeatures.2024<-function(CN_data, cores = 1, prePath="data/", allowedError = 0.1, rmNorm = FALSE)
{
    #get chromosome lengths
    chrlen<-read.table(paste0(prePath, "hg19.chrom.sizes.txt"),sep="\t",stringsAsFactors = F)[1:24,]
    
    #get centromere locations
    gaps<-read.table(paste0(prePath, "gap_hg19.txt"),sep="\t",header=F,stringsAsFactors = F)
    centromeres<-gaps[gaps[,8]=="centromere",]
    
    if(cores > 1) {
        require(foreach)
        require(doMC)

        registerDoMC(cores)

        temp_list = foreach(i=1:6) %dopar% {
            if(i == 1){
                list(segsize = getSegsize(CN_data, rmNorm = rmNorm) )
            } else if (i == 2) {
                list(bp10MB = getBPnum(CN_data,chrlen) )
            } else if (i == 3) {
                list(osCN = getOscilation.2024(CN_data,chrlen) )
            } else if (i == 4) {
                list(bpchrarm = getCentromereDistCounts(CN_data,centromeres,chrlen) )
            } else if (i == 5) {
                list(changepoint = getChangepointCN(CN_data, allowedError, rmNorm = rmNorm) )
            } else {
                list(copynumber = getCN(CN_data, rmNorm = rmNorm) )
            }
        
        }
        
        # Another failsafe that the outcome is definitely numeric
        # temp_list = unlist( temp_list, recursive = FALSE )
        outList = lapply(temp_list, function(thisDF) { 
            thisDF[,2] = as.numeric(thisDF[,2])
            return(thisDF)
        })
        return( outList )
        
    } else {  
        
        segsize<-getSegsize(CN_data, rmNorm = rmNorm)
        bp10MB<-getBPnum(CN_data,chrlen)
        osCN<-getOscilation.2024(CN_data,chrlen)
        bpchrarm<-getCentromereDistCounts(CN_data,centromeres,chrlen)
        changepoint<-getChangepointCN(CN_data,allowedError,rmNorm)
        copynumber<-getCN(CN_data, rmNorm = rmNorm)
        
        temp_list = list(segsize=segsize,bp10MB=bp10MB,osCN=osCN,bpchrarm=bpchrarm,changepoint=changepoint,copynumber=copynumber)
        # temp_list = unlist( temp_list, recursive = FALSE )
        outList = lapply(temp_list, function(thisDF) { 
            thisDF[,2] = as.numeric(thisDF[,2])
            return(thisDF)
        })
        return( outList )
        
    }

}

extractCopynumberFeatures.newchp<-function(CN_data, cores = 1, prePath="data/", allowedError = 0.1, rmNorm = FALSE)
{
    #get chromosome lengths
    chrlen<-read.table(paste0(prePath, "hg19.chrom.sizes.txt"),sep="\t",stringsAsFactors = F)[1:24,]
    
    #get centromere locations
    gaps<-read.table(paste0(prePath, "gap_hg19.txt"),sep="\t",header=F,stringsAsFactors = F)
    centromeres<-gaps[gaps[,8]=="centromere",]
    
    if(cores > 1) {
        require(foreach)
        require(doMC)
        
        registerDoMC(cores)
        
        temp_list = foreach(i=1:6) %dopar% {
            if(i == 1){
                list(segsize = getSegsize(CN_data, rmNorm = rmNorm) )
            } else if (i == 2) {
                list(bp10MB = getBPnum(CN_data,chrlen) )
            } else if (i == 3) {
                list(osCN = getOscilation.2024(CN_data,chrlen) )
            } else if (i == 4) {
                list(bpchrarm = getCentromereDistCounts(CN_data,centromeres,chrlen) )
            } else if (i == 5) {
                list(changepoint = getChangepointCN.new(CN_data, allowedError, rmNorm = rmNorm) )
            } else {
                list(copynumber = getCN(CN_data, rmNorm = rmNorm) )
            }
            
        }
        
        # Another failsafe that the outcome is definitely numeric
        # temp_list = unlist( temp_list, recursive = FALSE )
        outList = lapply(temp_list, function(thisDF) { 
            thisDF[,2] = as.numeric(thisDF[,2])
            return(thisDF)
        })
        return( outList )
        
    } else {  
        
        segsize<-getSegsize(CN_data, rmNorm = rmNorm)
        bp10MB<-getBPnum(CN_data,chrlen)
        osCN<-getOscilation.2024(CN_data,chrlen)
        bpchrarm<-getCentromereDistCounts(CN_data,centromeres,chrlen)
        changepoint<-getChangepointCN.new(CN_data,allowedError,rmNorm)
        copynumber<-getCN(CN_data, rmNorm = rmNorm)
        
        temp_list = list(segsize=segsize,bp10MB=bp10MB,osCN=osCN,bpchrarm=bpchrarm,changepoint=changepoint,copynumber=copynumber)
        # temp_list = unlist( temp_list, recursive = FALSE )
        outList = lapply(temp_list, function(thisDF) { 
            thisDF[,2] = as.numeric(thisDF[,2])
            return(thisDF)
        })
        return( outList )
        
    }
    
}

