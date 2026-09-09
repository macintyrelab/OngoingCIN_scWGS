
#### FUNCTION FOR PLOTTING HEATMAP OF SINGLE-CELL COPY NUMBER PROFILES #########

#' @title Heatmap representing copy number profiles from multiple cells/samples
#' @description 
#' @name plotSCHeatmap
#' @param cntable segment table with copy numbers. From QDNAseq objects, you can get this table using getSegTable() function
#' @param bin.size resolution in bp 
#' @param path path where you have saved the "uniqueEvents_functions.R" file, which contains functions needed
#' @param cluster add TRUE or FALSE to cluster or not the cells by similarity in the heatmap
#' @param cell_names add TRUE or FALSE to add or not the name of the cells in the heatmap
#' @param format add "png" or "pdf" 
#' @param outdir add path where you want to save the heatmap
#' @param outname name of the output file
#' @param dims add something in case you want to change plot high based on number of cells (default in NULL)
#'  
#' @return a heatmap with all copy number profiles
#' @export
#' @examples

plotSCHeatmap=function(cntable, bin.size, path, cluster, cell_names, format, outdir, outname, dims=NULL){
    
    ## 1) Unify copy number segments across cells
    ## Unified matrix
    # Define parameters
    window=as.numeric(bin.size*2)
    
    # Unify segments
    source(paste0(path,"/uniqueEvents_functions.R")) #upload functions
    UNIFY.MAT<-unifySegments.multiple(profiles=cntable,window=window,bin.size=bin.size)
    
    # Scale chromosome sizes => Number of heatmap columns proportional to chr size
    unify=chrScale(UNIFY.MAT,window)
    
    # Round copy number values
    unify[,4:ncol(unify)]=sapply(4:ncol(unify),function(x)round(unify[,x]))
    
    ## 2) Heatmap with copy number values
    # Set annotation
    ann=data.frame(unify$chromosome)
    colnames(ann)="chr"
    ann$chr=as.numeric(ann$chr)
    cols=c("#E046DC","#83A263","#64E399","#7E40E1","#E0AD4F","#E57155","#7BA9D9","#B56DDA","#D99CD8",
           "#878F8E","#DA9799","#C7E1E3","#CDE389","#6DC4DC","#E9D9B3","#85E658","#7C88DB","#DFC9E0",
           "#E55A9E","#6CE8DD","#D7E149","#B0E3C0")
    cols=stats::setNames(cols,unique(ann$chr))
    colAnn=ComplexHeatmap::HeatmapAnnotation(chr=ann$chr,col=list(chr=cols))
    
    # Plot as discrete value
    unify=t(unify[,4:ncol(unify)])
    colMain=c("0"="#3182BD","1"="#9ECAE1","2"="#E5E5E5","3"="#FDCC8A","4"="#FC8D59","5"="#E34A33",
              "6"="#B30000","7"="#980043","8"="#DD1C77","9"="#DF65B0","10"="#C994C7","11"="#D4B9DA",
              "12"="#BEAEEA","13"="#cee393","14"="#9bb060","15"="#616e3c")
    unify[unify>14]=15
    
    #Dimensions
    if(is.null(dims)){heigth=10}
    if(!is.null(dims)){heigth=nrow(unify)/2}
    
    if(format == "png"){
        png(file=paste0(outdir, "/Heatmap_unify_profiles_",outname,".png"), width = 14, height = heigth, units = "in", res = 300)
        print(Heatmap(unify, col = colMain, column_title = paste0("Copy number profiles of ", length(unique(cntable$sample))," cells unified"), 
                      cluster_rows = cluster, cluster_columns = FALSE,
                      top_annotation = colAnn, show_column_names = FALSE, 
                      show_row_names = cell_names, row_dend_width = unit(2, "cm")))
        dev.off()
    }
    
    if(format == "pdf"){
        pdf(file=paste0(outdir, "/Heatmap_unify_profiles_",outname,".pdf"), width = 14, height = heigth)
        print(Heatmap(unify, col = colMain, column_title = paste0("Copy number profiles of ", length(unique(cntable$sample))," cells unified"), 
                      cluster_rows = cluster, cluster_columns = FALSE,
                      top_annotation = colAnn, show_column_names = FALSE, 
                      show_row_names = cell_names, row_dend_width = unit(2, "cm")))
        dev.off()
    }
}


#' @title Function for scaling chromosome sizes in heatmap
#' @description 
#' @name chrScale
#' @param unify unified matrix with copy number profiles
#' @param window window used to split genome (bin.size*2)
#'  
#' @return a dataframe with data ready for heatmap
#' @export
#' @examples

chrScale=function(unify, window){
    # Get chromosome sizes & split into fixed windows
    chr_sizes=CNpare:::chr_sizes
    coordinates=as.data.frame(do.call(rbind,lapply(chr_sizes$chr, function(chr){
        start=seq(1,chr_sizes$length[chr_sizes$chr==chr],by=window)
        end=c(start[2:length(start)]-1,chr_sizes$length[chr_sizes$chr==chr])
        return(cbind(rep(chr,length(start)),cbind(start,end)))
    })))
    colnames(coordinates)=colnames(unify)[1:3]
    coordinates[,2:3]=sapply(2:3, function(x)as.numeric(coordinates[,x]))
    
    # Map unified segments with windows
    window.gr=GenomicRanges::GRanges(seqnames=coordinates$chromosome, IRanges(start=coordinates$start, end=coordinates$end))
    mat.gr=GenomicRanges::GRanges(seqnames=unify$chromosome, IRanges(start=unify$start, end=unify$end))
    overlaps=GenomicRanges::findOverlaps(window.gr, mat.gr)
    
    # Get copy number values
    df=matrix(data = NA, nrow = nrow(coordinates), ncol=ncol(unify)-3)
    colnames(df)=colnames(unify[4:ncol(unify)])
    for (i in 1:nrow(df)){
        idx = overlaps@to[overlaps@from == i]
        if(length(idx)!=0){
            df[i, ] <- colMeans(unify[idx, 4:ncol(unify)])
        }
    }
    df=cbind(coordinates,df)
    return(df[!is.na(df[,4]),])
}


