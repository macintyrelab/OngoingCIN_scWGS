################################################################################
## This script tests the approach for mapping events to specific CIN signatures
## Started: 10/03/25
## Updated: 02/09/26
##
## Here I perform simulations to model the acquisition of copy number alterations
## by a particular mutational process. Alterations have a specific configuration
## based on the definition of the causal signature. I generate 1000 samples with 
## different sample-level signature activities with 50-to-200 alterations. 
##
################################################################################

### LIBRARIES
rm(list=ls(all=TRUE))

library(this.path)
library(tidyverse)
library(data.table)
library(YAPSA)
library(NMF)
library(ggpubr)
library(ggthemes)

## PLOTS THEME
theme_set(theme_tufte(base_size = 6, base_family = "ArialMT"))
theme_update(text = element_text(size = 6),
             axis.text = element_text(size = 6),
             axis.line = element_line(size = 0.5), 
             axis.ticks = element_line(size = 0.5),
             axis.ticks.length = unit(.1, "cm"))


## PATHS
BASE=dirname(dirname(this.path()))
SIMS_DIR=dirname(this.path())
SC_DIR=file.path(BASE,"02_Derive_New_scCINSignatures")
INPUT_DIR=file.path(BASE,"data")
# Create directories for outputs
OUTPUT_DIR=file.path(SIMS_DIR,"outputs")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
PLOTS_DIR=file.path(SIMS_DIR,"plots")
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)

## LOAD DATA
DEFS=readRDS(file=paste0(SC_DIR, "/outputs/2_New_scCINSignatures_Definitions.rds"))

MIXMODELS=readRDS(file=paste0(INPUT_DIR, "/SC_MIXMODELS.rds"))

## PARAMETERS
min.cnas=50 #50 copy number alterations in each sample as minimum
max.cnas=200 #200 copy number alterations in each sample as maximum
samps=2000 #2000 samples
mean.active=5 #mean number of active mutational processes per sample

################################################################################
#### SIMULATE INPUT MATRIX #####################################################
set.seed(77777)

### 1. Create sample-level activities
# Signature activities are drawn from a uniform distribution and renormalised to sum
# to 1. Every mutational process, and every signature within a process, is equally
# likely to be active, and at most one signature per process is active in a sample
sigs=row.names(DEFS)

#mutational process behind each signature. Signatures with an unknown aetiology are
#pooled into a single process
PROCESS=list(chr.missegregation=c("CX1","CX6","CX14"),
             impaired.HR=c("CX2","CX3","CX5"),
             WGD="CX4",
             replication.stress=c("CX8","CX9","CX11"),
             ecDNA="CX13",
             impaired.NHEJ="CX10")
PROCESS$unknown=setdiff(sigs,unlist(PROCESS))

#number of active processes per sample: Poisson(mean.active) truncated to [2,#processes]
n.active=sample(2:length(PROCESS),samps,replace=TRUE,prob=dpois(2:length(PROCESS),mean.active))

sample_activities=do.call(rbind,lapply(1:samps,function(s){
    acts=setNames(numeric(length(sigs)),sigs)
    
    #processes active in this sample, all of them equally likely
    procs=sample(names(PROCESS),n.active[s])
    
    #one signature per active process, all of them equally likely
    active.sigs=sapply(procs,function(p){
        if(length(PROCESS[[p]])==1){PROCESS[[p]]}else{sample(PROCESS[[p]],1)}})
    
    #activity of each signature drawn uniformly, then renormalised to sum to 1
    drawn=runif(length(active.sigs))
    acts[active.sigs]=drawn/sum(drawn)
    
    return(acts)
}))
sample_activities=as.data.frame(sample_activities)
row.names(sample_activities)=paste0("sample",1:nrow(sample_activities))

### 2. Expected configurations based on high weights
#segment size
defs_size=DEFS[,1:20]
defs_size=defs_size/rowSums(defs_size)

#changepoint
defs_chp=DEFS[,21:28]
#signatures with 0-weights in all chpoint components (CX1, CX6, CX12): define chpoint
nochp=names(which(rowSums(defs_chp)==0))
defs_chp=defs_chp/rowSums(defs_chp)
for(s in nochp){defs_chp[s,]=c(0.5,0.5,rep(0,6))}

#bpwindow
defs_bpw=DEFS[,29:37]
defs_bpw=defs_bpw/rowSums(defs_bpw)

### 3. Generate copy number alterations in each sample
# draw the component of each event, given the signature causing it and a definition
draw=function(sig,defs) apply(defs[sig,,drop=FALSE],1,function(p)sample(length(p),1,prob=p))

## 3.1. Single-cell scenario
cntable=rbindlist(lapply(1:samps,function(s){
    print(paste0("Sample ",s))
    
    #activities
    acts=as.numeric(sample_activities[s,])
    
    #create events
    cnas=sample(min.cnas:max.cnas,1)
    
    ##signature causing each event
    sig=sample(sigs,cnas,replace=TRUE,prob=acts)
    
    ##configuration of each event
    data.table(sample=paste0("sample",s), event=1:cnas, signature=sig,
               segsize=MIXMODELS[["segsize"]]$Mean[draw(sig,defs_size)],
               changepoint=MIXMODELS[["changepoint"]]$Lambda[draw(sig,defs_chp)],
               bpwindow=MIXMODELS[["bpwindow"]]$Lambda[draw(sig,defs_bpw)])
}))
cntable=as.data.frame(cntable)


### 4. Create ECNF 
ecnf=list()
ecnf[["segsize"]]=as.data.frame(cbind(ID=cntable$sample,value=cntable$segsize))
ecnf[["changepoint"]]=as.data.frame(cbind(ID=cntable$sample,value=cntable$changepoint))
ecnf[["bpwindow"]]=as.data.frame(cbind(ID=cntable$sample,value=cntable$bpwindow))


### 5. Create ExC matrix
FEATS=names(MIXMODELS)

lMats=lapply(FEATS, function(thisFeature) {
    print(thisFeature)
    thisEcnf=ecnf[[thisFeature]]
    thisModel=MIXMODELS[[thisFeature]]
    
    # We want a posterior, hence likelihood (density) times prior (weight)
    dat=as.numeric(thisEcnf[,2])
    print("Calculating probabilities...")
    if(thisFeature!="segsize"){
        # Poisson model
        print("Poisson-based posterior")
        postDatUnscaled = sapply(1:nrow(thisModel), function(x) dpois(x = dat, lambda = thisModel[[x,"Lambda"]]))
    }
    else{
        # Gaussian model
        print("Gauss-based posterior")
        postDatUnscaled = sapply(1:nrow(thisModel), function(x) dnorm(x = dat, mean = thisModel[[x,"Mean"]], sd = thisModel[[x,"SD"]]))
    }
    # Normalize densities to probabilities
    postDatScaled = data.frame(postDatUnscaled/rowSums(postDatUnscaled))
    # Make sure that focal amplifications with huge number of copies have weight in the chp=8 
    # Make sure that windows with a lot of segments have weight in the bpw=8
    id = which(rowSums(is.na(postDatScaled)) == ncol(postDatScaled))
    if (thisFeature == "changepoint"){
        if(length(id)!=0){postDatScaled[id,]=matrix(rep(c(rep(0,7),1), length(id)), ncol = 8, nrow =length(id), byrow = T)}
    }
    if (thisFeature == "bpwindow"){
        if(length(id)!=0){postDatScaled[id,]=matrix(rep(c(rep(0,8),1), length(id)), ncol = 9, nrow =length(id), byrow = T)}
    }
    # Create component matrix
    colnames(postDatScaled) = paste0(thisFeature, 1:ncol(postDatScaled))
    postDatScaled=as.matrix(postDatScaled)
    
    return(postDatScaled)
})
names(lMats)=FEATS
ExC=do.call(cbind, lMats)


### 6. Map events to signatures
scdefs=cbind(defs_size,cbind(defs_chp,defs_bpw))

# score each event against every signature definition, keeping only the
# signatures that are active in the sample the event belongs to
active=as.matrix((sample_activities[cntable$sample,]>0)*1)
sig.scores=(ExC%*%t(scdefs))*active
sig.scores=sig.scores/rowSums(sig.scores)

events=cbind(cntable,as.data.frame(sig.scores))
events$map.sig=colnames(sig.scores)[max.col(sig.scores,ties.method="first")]


### 7. Test assignment accuracy
# Signature level
events$match=events$signature==events$map.sig
events$map.sig=factor(events$map.sig,levels = row.names(scdefs))
events$signature=factor(events$signature,levels = row.names(scdefs))
events$match=factor(events$match, levels = c("TRUE","FALSE"))

stats = events %>% group_by(map.sig,match, .drop=F) %>% summarise(n=n())
stats$freq = do.call(c,lapply(sigs,function(x){
    d=stats[stats$map.sig==x,]
    freq=d$n/sum(d$n)
    return(freq)
}))

# overall statistics
sum(stats$n[stats$match=="TRUE"]) #188352
sum(stats$n) #243184
mean(stats$freq[stats$match=="TRUE"], na.rm = TRUE) #77.15% of events correctly mapped!

# per signature
stats.sig = stats[stats$match=="TRUE",c(1,4)]

# per sample
stats.samp = events %>% group_by(sample,match,.drop=F) %>% summarise(n=n())
stats.samp$freq = do.call(c,lapply(unique(stats.samp$sample),function(x){
    d=stats.samp[stats.samp$sample==x,]
    freq=d$n/sum(d$n)
    return(freq)
}))
mean(stats.samp$freq[stats.samp$match==T]) #77.29%
sd(stats.samp$freq[stats.samp$match==T]) #12.90%
mean(stats.samp$freq[stats.samp$match==F]) #22.71%
sd(stats.samp$freq[stats.samp$match==F]) #12.90%


### 8. Plot mapping frequencies 

# plot per signature
df = events %>% group_by(map.sig,signature, .drop=F) %>% summarise(n=n())
df$freq = do.call(c,lapply(sigs,function(x){
    d=df[df$map.sig==x,]
    freq=(d$n/sum(d$n))*100
    return(freq)
}))

df=reshape::cast(df, signature ~ map.sig, mean, value = "freq")
row.names(df)=df[,1]
df = as.data.frame(df[,-1])

pdf(file=paste0(PLOTS_DIR, "/2.6_SuppFig_Assignment_Accuracy_Simulations_signatures.pdf"), height = 100/25.4, width = 150/25.4)
p=ComplexHeatmap::Heatmap(df, column_title = "Rate of copy number events assigned to signatures", cluster_rows = FALSE, cluster_columns = FALSE,
          cell_fun = function(j, i, x, y, width, height, fill) {
              grid.text(sprintf("%.1f", df[i, j]), x, y, gp = gpar(fontsize = 6))
          })
print(p)
dev.off()

