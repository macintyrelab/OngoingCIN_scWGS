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

#### MAIN FUNCTIONS #############################################################

extractThreeCopynumberFeatures.sc<-function(CN_data, cores, prePath, allowedError=0.1, rmNorm=TRUE, size=6000000, res=100000)
{
    #get chromosome lengths
    chrlen<-read.table(paste0(prePath, "/hg19.chrom.sizes.txt"),sep="\t",stringsAsFactors = F)[1:24,]
    
    if(cores > 1) {
        require(foreach)
        require(doMC)
        registerDoMC(cores)
        
        temp_list=foreach(i=1:3)%dopar%{
            if(i==1){
                list(segsize=getSegsize.sc(CN_data, rmNorm=rmNorm, res=res))
            }else if(i==2){
                list(changepoint=getChangepointCN.sc(CN_data, allowedError=allowedError, rmNorm=rmNorm))
            }else{
                list(bpwindow=getBPwindow(CN_data, rmNorm=rmNorm, size=size))
            }
        }
        outList=lapply(temp_list, function(thisDF){ 
            thisDF[,2] = as.numeric(thisDF[,2])
            return(thisDF)
        })
        return(outList)
        
    } else {
        #SEGSIZE
        segsize<-getSegsize.sc(CN_data, rmNorm=rmNorm, res=res)
        #CHANGEPOINT
        changepoint<-getChangepointCN.sc(CN_data, allowedError, rmNorm=rmNorm)
        #BPWINDOW
        bpwindow=getBPwindow(CN_data, rmNorm=rmNorm, size=size)
        
        temp_list=list(segsize=segsize,changepoint=changepoint,bpwindow=bpwindow)
        outList=lapply(temp_list, function(thisDF) { 
            thisDF[,2]=as.numeric(thisDF[,2])
            return(thisDF)
        })
        return(outList)
    }
}

extractThreeCopynumberFeatures.new<-function(CN_data, cores, prePath, allowedError=0.1, rmNorm=TRUE, size=6000000)
{
    #get chromosome lengths
    chrlen<-read.table(paste0(prePath, "/hg19.chrom.sizes.txt"),sep="\t",stringsAsFactors = F)[1:24,]
    
    if(cores > 1) {
        require(foreach)
        require(doMC)
        registerDoMC(cores)
        
        temp_list=foreach(i=1:3)%dopar%{
            if(i==1){
                list(segsize=getSegsize(CN_data, rmNorm=rmNorm))
            }else if(i==2){
                list(changepoint=getChangepointCN.new(CN_data, allowedError=allowedError, rmNorm=rmNorm))
            }else{
                list(bpwindow=getBPwindow(CN_data, rmNorm=rmNorm, size=size))
            }
        }
        outList=lapply(temp_list, function(thisDF){ 
            thisDF[,2] = as.numeric(thisDF[,2])
            return(thisDF)
        })
        return(outList)
        
    } else {
        #SEGSIZE
        segsize<-getSegsize(CN_data, rmNorm=rmNorm)
        #CHANGEPOINT
        changepoint<-getChangepointCN.new(CN_data, allowedError, rmNorm=rmNorm)
        #BPWINDOW
        bpwindow=getBPwindow(CN_data, rmNorm=rmNorm, size=size)
        
        temp_list=list(segsize=segsize,changepoint=changepoint,bpwindow=bpwindow)
        outList=lapply(temp_list, function(thisDF) { 
            thisDF[,2]=as.numeric(thisDF[,2])
            return(thisDF)
        })
        return(outList)
    }
}



##### HELPER FUNCTIONS #########################################################
getSegsize.sc<-function(abs_profiles, rmNorm = FALSE, res = 100000)
{
    out<-c()
    samps<-getSampNames(abs_profiles)
    for(i in samps)
    {
        if(class(abs_profiles)=="QDNAseqCopyNumbers")
        {
            segTab<-getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
        }
        else
        {
            segTab<-abs_profiles[[i]]
            colnames(segTab)[4]<-"segVal"
        }
        # First tap.
        segTab$segVal = as.numeric(segTab$segVal)
        # If wished, don't consider normal segments
        if(rmNorm) { segTab = segTab[ segTab$segVal != 2, ] }
        # Avoiding potential artefact
        segTab$segVal[segTab$segVal<0]<-0
        seglen<-segTab$end-segTab$start
        seglen<-seglen[seglen>0]
        # match single-cell resolution
        seglen[seglen<res]=res
        # Double tap.
        out<-rbind(out,cbind(ID=rep(i,length(seglen)),value=as.numeric(seglen)))
    }
    rownames(out)<-NULL
    data.frame(out,stringsAsFactors = F)
}

getSegsize<-function(abs_profiles, rmNorm = FALSE)
{
    out<-c()
    samps<-getSampNames(abs_profiles)
    for(i in samps)
    {
        if(class(abs_profiles)=="QDNAseqCopyNumbers")
        {
            segTab<-getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
        }
        else
        {
            segTab<-abs_profiles[[i]]
            colnames(segTab)[4]<-"segVal"
        }
        # First tap.
        segTab$segVal = as.numeric(segTab$segVal)
        # If wished, don't consider normal segments
        if(rmNorm) { segTab = segTab[ segTab$segVal != 2, ] }
        # Avoiding potential artefact
        segTab$segVal[segTab$segVal<0]<-0
        seglen<-segTab$end-segTab$start
        seglen<-seglen[seglen>0]
        # Double tap.
        out<-rbind(out,cbind(ID=rep(i,length(seglen)),value=as.numeric(seglen)))
    }
    rownames(out)<-NULL
    data.frame(out,stringsAsFactors = F)
}

getChangepointCN.sc<-function(abs_profiles, allowedError = 0.1, rmNorm = FALSE)
{
    out<-c()
    samps<-getSampNames(abs_profiles)
    for(i in samps)
    {
        if(class(abs_profiles)=="QDNAseqCopyNumbers") {
            segTab<-getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
        } else {
            segTab<-abs_profiles[[i]]
            colnames(segTab)[4]<-"segVal"
        }
        segTab$segVal = as.numeric(segTab$segVal)
        # Avoiding potential artefact
        segTab$segVal[segTab$segVal<0]<-0 #make sure any segVal is negtive
        segTab<-segTab[segTab$start!=segTab$end,] #make sure that there are not any segment with length=0
        
        # Compute by chromosome
        chrs<-unique(segTab$chromosome)
        allcp<-c()
        for(c in chrs) {
            tab=segTab[segTab$chromosome==c,]
            
            ## Changepoint for events from right
            r.chp=c()
            currseg<-as.numeric(tab$segVal)
            if(length(currseg)>0){
                firstSeg = abs(2 - currseg[1])
                # As we look only at the left side of a CNA, we might miss a changepoint at the beginning of the p-arm
                # That's why we check manually but only regard this value if it is higher than an allowed error rate.
                if(firstSeg <= allowedError) {
                    theseChanges = currseg[-1]-currseg[-length(currseg)]
                    if(rmNorm) { theseChanges = theseChanges[ currseg[-1] != 2 ] }
                    r.chp<-c(r.chp, theseChanges)
                } else {
                    theseChanges = c( firstSeg, currseg[-1]-currseg[-length(currseg)] )
                    if(rmNorm) { theseChanges = theseChanges[ currseg != 2 ] }
                    r.chp<-c(r.chp, theseChanges)
                }
            }
            
            ## Changepoint for events from left
            l.chp=c()
            currseg<-rev(as.numeric(tab$segVal)) #reverse to compute the change point from the chromosome end to the centromere
            if(length(currseg)>0){
                firstSeg = abs(2 - currseg[1])
                if(firstSeg <= allowedError) {
                    theseChanges = currseg[-1]-currseg[-length(currseg)]
                    if(rmNorm) { theseChanges = theseChanges[ currseg[-1] != 2 ] }
                    l.chp<-c(l.chp, theseChanges)
                } else {
                    theseChanges = c( firstSeg, currseg[-1]-currseg[-length(currseg)] )
                    if(rmNorm) { theseChanges = theseChanges[ currseg != 2 ] }
                    l.chp<-c(l.chp, theseChanges)
                }
            }
            
            ## Select correct changepoint 
            # Huge negative copy number changes are due to the previous focal amplification
            # Getting the absolute copychange inflates the signal of replication stress signatures
            b.chp=cbind(rev(currseg[currseg!=2]),r.chp,rev(l.chp))
            chp=ifelse(b.chp[,2]>=-3,b.chp[,2],b.chp[,3])
            # Correct cases with focal amplifications in tandem by subtracting normal segVal
            chp=ifelse(chp>=-3,abs(chp),abs(b.chp[,1]-2)) 
            
            ## Adapt to single-cell scenario
            # Avoid changepoints lower than 1
            chp[chp<1]=1
            # Round chp to the closest integer
            chp=round(chp)
            allcp<-c(allcp,unname(chp))
        }
        # Make sure it's really numeric
        out<-rbind(out,cbind(ID=rep(i,length(allcp)),value=as.numeric(allcp)))
    }
    rownames(out)<-NULL
    data.frame(out,stringsAsFactors = F)
}

getChangepointCN.new<-function(abs_profiles, allowedError = 0.1, rmNorm = FALSE)
{
    out<-c()
    samps<-getSampNames(abs_profiles)
    for(i in samps)
    {
        if(class(abs_profiles)=="QDNAseqCopyNumbers") {
            segTab<-getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
        } else {
            segTab<-abs_profiles[[i]]
            colnames(segTab)[4]<-"segVal"
        }
        segTab$segVal = as.numeric(segTab$segVal)
        # Avoiding potential artefact
        segTab$segVal[segTab$segVal<0]<-0 #make sure any segVal is negtive
        segTab<-segTab[segTab$start!=segTab$end,] #make sure that there are not any segment with length=0
        
        # Compute by chromosome
        chrs<-unique(segTab$chromosome)
        allcp<-c()
        for(c in chrs) {
            tab=segTab[segTab$chromosome==c,]
            
            ## Changepoint for events from right
            r.chp=c()
            currseg<-as.numeric(tab$segVal)
            if(length(currseg)>0){
                firstSeg = abs(2 - currseg[1])
                # As we look only at the left side of a CNA, we might miss a changepoint at the beginning of the p-arm
                # That's why we check manually but only regard this value if it is higher than an allowed error rate.
                if(firstSeg <= allowedError) {
                    theseChanges = currseg[-1]-currseg[-length(currseg)]
                    if(rmNorm) { theseChanges = theseChanges[ currseg[-1] != 2 ] }
                    r.chp<-c(r.chp, theseChanges)
                } else {
                    theseChanges = c( firstSeg, currseg[-1]-currseg[-length(currseg)] )
                    if(rmNorm) { theseChanges = theseChanges[ currseg != 2 ] }
                    r.chp<-c(r.chp, theseChanges)
                }
            }
            
            ## Changepoint for events from left
            l.chp=c()
            currseg<-rev(as.numeric(tab$segVal)) #reverse to compute the change point from the chromosome end to the centromere
            if(length(currseg)>0){
                firstSeg = abs(2 - currseg[1])
                if(firstSeg <= allowedError) {
                    theseChanges = currseg[-1]-currseg[-length(currseg)]
                    if(rmNorm) { theseChanges = theseChanges[ currseg[-1] != 2 ] }
                    l.chp<-c(l.chp, theseChanges)
                } else {
                    theseChanges = c( firstSeg, currseg[-1]-currseg[-length(currseg)] )
                    if(rmNorm) { theseChanges = theseChanges[ currseg != 2 ] }
                    l.chp<-c(l.chp, theseChanges)
                }
            }
            
            ## Select correct changepoint 
            # Huge negative copy number changes are due to the previous focal amplification
            # Getting the absolute copychange inflates the signal of replication stress signatures
            b.chp=cbind(rev(currseg[currseg!=2]),r.chp,rev(l.chp))
            chp=ifelse(b.chp[,2]>=-3,b.chp[,2],b.chp[,3])
            # Correct cases with focal amplifications in tandem by subtracting normal segVal
            chp=ifelse(chp>=-3,abs(chp),abs(b.chp[,1]-2)) 
            allcp<-c(allcp,unname(chp))
        }
        # Make sure it's really numeric
        out<-rbind(out,cbind(ID=rep(i,length(allcp)),value=as.numeric(allcp)))
    }
    rownames(out)<-NULL
    data.frame(out,stringsAsFactors = F)
}

getBPwindow<-function(abs_profiles, rmNorm = FALSE, size = 3000000)
{
    out<-c()
    samps<-getSampNames(abs_profiles)
    for(i in samps)
    {
        if(class(abs_profiles)=="QDNAseqCopyNumbers"){
            segTab<-getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
        }else{
            segTab<-abs_profiles[[i]]
            colnames(segTab)[4]<-"segVal"
        }
        segTab$segVal = as.numeric(segTab$segVal)
        
        #Avoiding potential artefact
        segTab$segVal[segTab$segVal<0]<-0 #make sure any segVal is negative
        segTab<-segTab[segTab$start!=segTab$end,] #make sure that there are not any segment with length=0
        
        # Get the number of events contiguous to the segment in a window
        chrs<-unique(segTab$chromosome)
        allnum<-unlist(lapply(chrs,function(c){
            chrTab<-segTab[segTab$chromosome==c,]
            events<-1:nrow(chrTab)
            #If wished, don't consider normal segments
            if(rmNorm){events<-which(chrTab$segVal!=2)}
            #If no event in the chromosome
            if(length(events)==0){
                allnum=NULL
            }
            #If only one segment in the chromosome, number of contiguous events is 0
            if(length(events)==1){
                allnum=0
            }
            #If more than one segment in the chromosome => take the maximum number of the right and left hand-side
            if(length(events)>1){
                allnum<-unlist(lapply(events,function(e){
                    wup=chrTab$end[e]+size
                    nup=(chrTab$start<=wup&chrTab$end>chrTab$end[e])
                    nup=length(nup[nup==TRUE])
                    wdown=chrTab$start[e]-size
                    ndown=(chrTab$end>=wdown&chrTab$start<chrTab$start[e])
                    ndown=length(ndown[ndown==TRUE])
                    thisNum=max(nup,ndown)
                    return(thisNum)
                }))
            }
            return(allnum)
        }))
        out<-rbind(out,cbind(ID=rep(i,length(allnum)),value=allnum))
    }
    rownames(out)<-NULL
    out<-data.frame(out,stringsAsFactors = F)
    # Make sure it's really numeric
    out$value=as.numeric(out$value)
    return(out)
}


getKLD=function(matrix1, matrix2){
    kldiv=lapply(1:nrow(matrix1),function(i){
        x=rbind(matrix1[i,],matrix2[i,])
        kld=suppressMessages(philentropy::KL(x, unit="log")) #equal to my calculations
        return(kld)
    })
    kldiv=sum(do.call(rbind,kldiv)[,1])
    return(kldiv)
}

getCosSignatures=function(signatures){
    if(nrow(signatures)==1){
        cos=1
    }else{
        cos=mean(unlist(lapply(1:nrow(signatures),function(i){
            x=signatures[i,]
            oth=signatures[-i,]
            if(is.null(nrow(oth))){
                sim=lsa::cosine(x, oth)[1,1]
            }else{
                sim=unlist(lapply(1:nrow(oth), function(thisSig){cos=lsa::cosine(x, oth[thisSig,])}))
            }
            return(sim)
        })))
    }
    return(cos)
}

cosineSimForMatrices <- function(X, corr=FALSE){
    if(corr){ X = apply(X, 2, function(x){ x-mean(x) }) }
    denom = solve(diag(sqrt(diag(t(X)%*%X))))
    return( denom%*%(t(X)%*%X)%*%denom )
}

getSampNames<-function(abs_profiles)
{
    if(class(abs_profiles)=="QDNAseqCopyNumbers")
    {
        samps<-colnames(abs_profiles)
    }
    else
    {
        samps<-names(abs_profiles)
    }
    samps
}


getSegTable<-function(x)
{
    dat<-x
    sn<-Biobase::assayDataElement(dat,"segmented")
    fd <- Biobase::fData(dat)
    fd$use -> use
    fdfiltfull<-fd[use,]
    sn<-sn[use,]
    segTable<-c()
    for(c in unique(fdfiltfull$chromosome))
    {
        snfilt<-sn[fdfiltfull$chromosome==c]
        fdfilt<-fdfiltfull[fdfiltfull$chromosome==c,]
        sn.rle<-rle(snfilt)
        starts <- cumsum(c(1, sn.rle$lengths[-length(sn.rle$lengths)]))
        ends <- cumsum(sn.rle$lengths)
        lapply(1:length(sn.rle$lengths), function(s) {
            from <- fdfilt$start[starts[s]]
            to <- fdfilt$end[ends[s]]
            segValue <- sn.rle$value[s]
            c(fdfilt$chromosome[starts[s]], from, to, segValue)
        }) -> segtmp
        segTableRaw <- data.frame(matrix(unlist(segtmp), ncol=4, byrow=T),stringsAsFactors=F)
        segTable<-rbind(segTable,segTableRaw)
    }
    colnames(segTable) <- c("chromosome", "start", "end", "segVal")
    segTable
}


getPloidy<-function(abs_profiles)
{
    out<-c()
    samps<-getSampNames(abs_profiles)
    for(i in samps)
    {
        if(class(abs_profiles)=="QDNAseqCopyNumbers")
        {
            segTab<-getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
        }
        else
        {
            segTab<-abs_profiles[[i]]
            colnames(segTab)[4]<-"segVal"
        }
        segLen<-(as.numeric(segTab$end)-as.numeric(segTab$start))
        ploidy<-sum((segLen/sum(segLen))*as.numeric(segTab$segVal))
        out<-c(out,ploidy)
    }
    data.frame(out,stringsAsFactors = F)
}

# Function for getting copy number tables from QDNAseq object
getCNTable<-function(cnobject, WIGGLE=0.1, IGNOREDELS=FALSE, REMOVEQUIET=FALSE){
    
    ## Step 1: Get segment tables
    samps=colnames(cnobject)
    CN.TABLE=as.data.frame(do.call(rbind,lapply(samps, function(thisClon){
        segTab=CNpare::getSegTable(cnobject[,thisClon])
        segTab$sample=thisClon
        segTab[,2:4]=sapply(2:4,function(x)as.numeric(segTab[,x]))
        return(segTab)
    })))
    
    ## Step 2: Smooth normal segments
    # Get normal segments
    CN.TABLE$segVal[CN.TABLE$segVal > (2-WIGGLE) & CN.TABLE$segVal < (2+WIGGLE)] = 2
    # Merge segments only when segments with same copy number value follow each other 
    df=idSmoothingTargets(CN.TABLE, WIGGLE = 0, colNameSegVal = "segVal", colNameChr = "chromosome", IGNOREDELS = IGNOREDELS)
    # Split by sample name
    lRaw=split(df, df$sample)
    # Smooth segments by taking the weighted average of the segVal and their lengths
    lSmooth=smoothSegments(lRaw, CORES = 1, SMOOTHINGFACTOR = 0, colNameMerge = "segVal", colNameChr = "chromosome",
                           colNameStart = "start", colNameEnd = "end", IGNOREDELS = IGNOREDELS, asDf = FALSE)
    
    CN.TABLE=data.frame(rbindlist(lSmooth))
    return(CN.TABLE)
}

idSmoothingTargets = function(dfAllSegs, WIGGLE, colNameSegVal, colNameChr, IGNOREDELS = TRUE) {
    
    ### Check column name
    testSegVal = dfAllSegs[[colNameSegVal]][1]
    testChr = dfAllSegs[[colNameChr]][1]
    if(!is.numeric(testSegVal)){stop("Segment Value column has no numeric value in it. Supplied correct column name? Forgot conversion?")}
    if(is.null(testSegVal)){stop("Chromosome column has no numeric value in it. Supplied correct column name?")}
    
    # take differences to segment down below
    dfAllSegs$diffs = c(abs(dfAllSegs[[colNameSegVal]][1:(nrow(dfAllSegs)-1)] - dfAllSegs[[colNameSegVal]][2:nrow(dfAllSegs)]), WIGGLE+1)
    # set TRUE if difference to next segment is smaller than the user supplied cutoff
    dfAllSegs$smooth = dfAllSegs$diffs <= WIGGLE
    # set all segments which are last in a chromosome to FALSE. This also prevents leaking to other samples and cohorts.
    dfAllSegs$smooth[cumsum(rle(as.character(dfAllSegs[[colNameChr]]))$lengths)] = FALSE
    
    # Ignore deletions if wished
    if(IGNOREDELS) {dfAllSegs$smooth[dfAllSegs[[colNameSegVal]] == 0] = FALSE}
    
    return(dfAllSegs)
}


smoothSegments = function(lRaw, CORES, SMOOTHINGFACTOR, colNameMerge, colNameChr, colNameStart, colNameEnd,
                          IGNOREDELS = TRUE, asDf = FALSE) {
    
    ### Check column names
    test = lRaw[[1]]
    testMerge = test[[colNameMerge]][1]
    testChr = test[[colNameChr]][1]
    testStart = test[[colNameStart]][1]
    testEnd = test[[colNameEnd]][1]
    if(! is.numeric(testMerge)) { stop("Merge column has no numeric value in it. Supplied correct column name?")}
    if(is.null(testChr)) { stop("Chromosome column has no numeric value in it. Supplied correct column name?")}
    if(! is.numeric(testStart)) { stop("Start column has no numeric value in it. Supplied correct column name?")}
    if(! is.numeric(testEnd)) { stop("End column has no numeric value in it. Supplied correct column name?")}
    
    # Add diff column to names we want to keep when merging (comes from function "idSmoothingTargets").
    colNameMerge = c(colNameMerge, "diffs")
    
    registerDoMC(CORES)
    lSmooth = foreach(thisSample = lRaw, .final = function(x) setNames(x, names(lRaw)) ) %dopar% {
        
        thisOut = thisSample
        stillSmoothing = sum(thisOut$smooth)
        while( stillSmoothing > 0 ) {
            # For the while loop:
            # Read lines from thisSample and change in thisOut. Hence for a new iteration I need to sync the two.
            thisSample = thisOut
            
            rleRaw = rle(thisSample$smooth)
            # This takes the indeces of the FALSE chains and adds 1. This should give you the next segment which is TRUE.
            # Two challenges:
            # 1) Last segment always FALSE (see above), hence removal of the last number as this would indicate to a segment outside the df.
            # 2) If it starts with a TRUE segment, this would not be found when looking at the FALSE chains. Hence, adding index 1 manually if chain starts with TRUE.
            indRaw = cumsum(rleRaw$lengths)[ ! rleRaw$values ] + 1
            indRaw = indRaw[ -length(indRaw) ]
            if( rleRaw$values[1] ) { indRaw = c(1, indRaw) }
            
            # loop over start indices of TRUE chains.
            for(i in indRaw) {
                # detect length of segments to smooth. add 1 as the last segment has a FALSE value in it but still belongs to this chain.
                endOfStreak = i + rle(thisSample$smooth[i:nrow(thisSample)])$lengths[1]
                # extract reads
                dfMerge = thisSample[i:endOfStreak,]
                
                # too stupid to make this work with data.table
                newElement = as.data.frame( dfMerge[1,] )
                # Get new end and check first wether valid number.
                newEnd = dfMerge[nrow(dfMerge),][[colNameEnd]]
                if(! is.null(newEnd)) {
                    newElement[[colNameEnd]] = newEnd
                } else {
                    stop("New end coordinate is null. Supplied correct column name?")
                }
                ## Column "segVal" will be dealt with in a minute. Column "diffs" later when running again idSmoothingTargets.
                
                # Merge cn specifically by taking the length of the elements into consideration
                widthWeights = dfMerge[[colNameEnd]] - dfMerge[[colNameStart]]
                # save those cases where length of the segment is zero in both segs to merge (that happens when SV breaks are included in the release)
                if(widthWeights[1]==0 && widthWeights[2]==0){
                    newElement[[colNameMerge[1]]] = mean(dfMerge[[colNameMerge[1]]])
                } else {
                    newElement[[colNameMerge[1]]] = weighted.mean(dfMerge[[colNameMerge[1]]], widthWeights)
                }
                
                # Replace all to merge segments with the new merged segment. Later delete duplicated.
                thisOut[i:endOfStreak,] = newElement
            }
            
            # as we have replaced all segments with the new mean segment, we need to remove the duplicates
            thisOut = thisOut[ ! duplicated(thisOut), ]
            # again detect segments which needs smoothing
            thisOut = idSmoothingTargets(thisOut, SMOOTHINGFACTOR, colNameSegVal = colNameMerge[[1]], colNameChr = colNameChr,
                                         IGNOREDELS = IGNOREDELS)
            stillSmoothing = sum(thisOut$smooth)
        }
        
        # after smoothing is finished, change name of cohort
        thisOut$smooth = NULL
        thisOut$diffs = NULL
        return( thisOut )
    }
    
    if( isTRUE(asDf) ) {
        dfSmooth = setDT( rbindlist( lSmooth ) )
        return( dfSmooth )
    } else {
        return( lSmooth )
    }
}



