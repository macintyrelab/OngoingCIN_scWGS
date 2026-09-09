################################################################################
## False-unique events: how many cells in the FINAL population carry each one
## Based on: 01b_unique_fp_detection.R
## Started: 19/06/26
##
## For each subsampling fraction we re-run the unique-event detection (exactly as
## in 01b), keep the events that are called "unique" in the subsample but are NOT
## truly unique in the full population (i.e. the FALSE unique events), and for
## each of them count how many cells of the final population actually carry that
## event. "Carry" = a subclonal event in data/all_events_minussi.rds that overlaps
## the false event's coordinates and has the same (rounded) copy-number value.
##
## Output (per model): a long data.frame, one row per false-unique detection:
##   model, iter, ncells, ntotal, chromosome, start, end, segVal, id, n_cells_in_pop
##
################################################################################

## strict numeric coercion
num <- function(x) suppressWarnings(as.numeric(as.character(x)))

## ---------- Libraries ----------
suppressPackageStartupMessages({
    library(this.path)
    library(data.table)
    library(GenomicRanges)
    library(foreach)
    library(doMC)
    library(dplyr)
})
options(stringsAsFactors = FALSE)

## be explicit about threads
data.table::setDTthreads(1)

## ---------- Parallel cores ----------
cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "10"))
registerDoMC(cores = cores)

## ---------- Paths ----------
BASE      <- dirname(dirname(this.path()))
MISC_DIR  <- file.path(BASE,"07_Miscellaneous")
SCELL_DIR <- file.path(BASE,"04_Process_scWGS")
INPUT_DIR <- file.path(BASE,"data")

OUTPUT_DIR <- file.path(MISC_DIR,"outputs")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

## ---------- Functions ----------
source(file.path(BASE, "scripts/uniqueEvents_functions.R"))
source(file.path(BASE, "scripts/newCNSigs_functions.R"))
source(file.path(BASE, "scripts/main_functions.R"))

## ---------- Parameters ----------
set.seed(77777)
THRESH      <- 0.90
SIGN.THRESH <- 0.01
WINDOW      <- 100e3
RMNORM      <- TRUE
CORES_FEAT  <- 1L
SIZE        <- 5e6
RESOLUTION  <- "100kb"
BIN.SIZE    <- as.numeric(sub("kb$","",RESOLUTION))*1000

## ============================================================================
##                       LOAD INPUT DATA
## ============================================================================

message("Loading precomputed tables...")

# 1. Lists of Long DataFrames
BIN.TABLE <- readRDS(paste0(INPUT_DIR, "/1.1.1_Minussi_Bins.rds"))
SEG.TABLE <- readRDS(paste0(INPUT_DIR, "/1.1.1_Minussi_Segments.rds"))

# 2. Mixture models and nulls
MIXMODELS      <- readRDS(file.path(BASE,"02_Derive_New_scCINSignatures/data/SC_MIXMODELS.rds"))
lDistributions <- readRDS(file.path(SCELL_DIR, "intermediate/confidence/NullDistribution_Confidence_Segsize.rds"))
dfTrueUnique   <- readRDS(file.path(OUTPUT_DIR,"real_unique_minussi.rds"))
dfTrueUnique$id <- paste0(dfTrueUnique$chromosome,":",
                          format(dfTrueUnique$start, scientific = FALSE, trim = TRUE),"-",
                          format(dfTrueUnique$end,   scientific = FALSE, trim = TRUE))

# 3. Per-cell classified events in the FINAL population (from colleague).
#    list per model: data.frame(chromosome,start,end,segVal,sample,class)
ALL.EVENTS <- readRDS(file.path(OUTPUT_DIR, "all_events_minussi.rds"))


## ============================================================================
##                       HELPER FUNCTIONS
## ============================================================================
mean_bin_val_over_range <- function(chr, s, e, target_sample, current_ids, bin_matrix, bins_gr) {
    
    segment.gr <- GRanges(seqnames = chr, ranges = IRanges(start = s, end = e))
    overlaps <- findOverlaps(segment.gr, bins_gr)
    idx <- subjectHits(overlaps)
    
    target_mat <- bin_matrix[idx, target_sample, drop = FALSE]
    target_vec <- suppressWarnings(as.numeric(target_mat[, 1]))
    rawV <- mean(target_vec, na.rm = TRUE)
    
    oth_ids <- intersect(current_ids, colnames(bin_matrix))
    oth_ids <- setdiff(oth_ids, target_sample)
    
    if (!length(oth_ids)) {
        return(list(rawVal = rawV, others = numeric(0)))
    }
    
    others_mat <- bin_matrix[idx, oth_ids, drop = FALSE]
    
    if (nrow(others_mat) == 0) {
        return(list(rawVal = rawV, others = numeric(0)))
    }
    
    others_mat_num <- suppressWarnings(
        matrix(as.numeric(others_mat),
               nrow = nrow(others_mat),
               ncol = ncol(others_mat))
    )
    colnames(others_mat_num) <- colnames(others_mat)
    
    if (ncol(others_mat_num) == 1) {
        others <- mean(others_mat_num[, 1], na.rm = TRUE)
        names(others) <- colnames(others_mat_num)[1]
    } else {
        others <- colMeans(others_mat_num, na.rm = TRUE)
    }
    
    list(rawVal = rawV, others = others)
}

## Build the full-population SUBCLONAL reference (GRanges) for a model.
## Step 1 of the colleague's recipe: drop "unique", keep "subclonal".
build_subclonal_ref <- function(df) {
    sub <- df[df$class == "subclonal", c("chromosome","start","end","segVal","sample")]
    if (nrow(sub) == 0) return(NULL)
    GRanges(seqnames = as.character(sub$chromosome),
            ranges   = IRanges(start = as.numeric(sub$start), end = as.numeric(sub$end)),
            segVal   = round(as.numeric(sub$segVal)),
            sample   = as.character(sub$sample))
}

## Count cells in the final population carrying a false-unique event.
## Map the event coordinates onto the subclonal events (boundaries may differ, so
## we use coordinate OVERLAP), keep the cells whose overlapping subclonal event has
## the same rounded copy number, and count the distinct cells.
## (A single false segment may overlap several subclonal events; a cell is counted
##  if ANY of its overlapping subclonal events matches the copy number.)
count_cells_in_pop <- function(chr, s, e, cn, ref_gr) {
    if (is.null(ref_gr)) return(NA_integer_)
    q  <- GRanges(seqnames = as.character(chr), ranges = IRanges(start = as.numeric(s), end = as.numeric(e)))
    ov <- findOverlaps(q, ref_gr)
    if (length(ov) == 0) return(0L)
    hit <- ref_gr[subjectHits(ov)]
    hit <- hit[hit$segVal == round(as.numeric(cn))]
    length(unique(hit$sample))
}

## ============================================================================
##                          MAIN SIMULATIONS
## ============================================================================

niters <- 1:20

MODELS_TO_RUN <- intersect(names(BIN.TABLE), names(SEG.TABLE))
# Filter for specific models if desired
# MODELS_TO_RUN <- "TN1"

for (model in MODELS_TO_RUN) {
    message("\n=== Model: ", model, " ===")
    
    ## -----------------------------------------------------------------------
    ## 1. PREPARE DATA FOR THIS MODEL (With Intersection Check)
    ## -----------------------------------------------------------------------
    
    DT_BINS <- as.data.table(BIN.TABLE[[model]])
    DT_SEGS <- as.data.table(SEG.TABLE[[model]])
    
    if(is.null(DT_BINS$binVal)) stop("BIN.TABLE missing 'binVal' column")
    
    message("  Constructing bin matrix (Long to Wide)...")
    MAT_DT <- data.table::dcast(DT_BINS, chromosome + start + end ~ sample, value.var = "binVal")
    MAT_DT <- MAT_DT[order(as.integer(gsub("X|Y", "23", chromosome)), start),]
    
    MODEL_BINS_GR <- GRanges(seqnames = MAT_DT$chromosome,
                             ranges   = IRanges(start = MAT_DT$start, end = MAT_DT$end))
    
    sample_cols <- setdiff(colnames(MAT_DT), c("chromosome", "start", "end"))
    MODEL_BIN_MAT <- as.matrix(MAT_DT[, ..sample_cols])
    
    rm(MAT_DT); gc()
    
    # --- CRITICAL FIX: Intersect Samples ---
    cells_in_bins <- colnames(MODEL_BIN_MAT)
    cells_in_segs <- unique(DT_SEGS$sample)
    
    common_cells <- intersect(cells_in_bins, cells_in_segs)
    
    message(sprintf("  Samples in Bins: %d | Samples in Segs: %d | Common: %d",
                    length(cells_in_bins), length(cells_in_segs), length(common_cells)))
    
    if(length(common_cells) < 2) {
        message("  Skipping model ", model, " (Not enough common samples).")
        next
    }
    
    # Update datasets to only include the common cells
    MODEL_BIN_MAT <- MODEL_BIN_MAT[, common_cells, drop=FALSE]
    DT_SEGS       <- DT_SEGS[sample %in% common_cells]
    all_cells     <- common_cells
    
    ## Full-population subclonal reference for this model (for cell counting)
    REF_GR <- if (model %in% names(ALL.EVENTS)) build_subclonal_ref(ALL.EVENTS[[model]]) else NULL
    if (is.null(REF_GR)) message("  WARNING: no subclonal reference for model ", model,
                                 " -> n_cells_in_pop will be NA")
    
    ## -----------------------------------------------------------------------
    ## 2. DEFINE ITERATION GRID
    ## -----------------------------------------------------------------------
    n_total <- length(all_cells)
    
    if (n_total > 150) {
        ncells <- c(seq(2, 100), seq(105, n_total, by = 5))
    } else if (n_total < 100) {
        ncells <- 2:n_total
    } else {
        ncells <- 2:100
    }
    
    ## -----------------------------------------------------------------------
    ## 3. RUN SIMULATION LOOPS  (collect FALSE-unique events per iter/ncells)
    ## -----------------------------------------------------------------------
    EMPTY <- data.table(iter = integer(), ncells = integer(),
                        chromosome = character(), start = numeric(),
                        end = numeric(), segVal = numeric(), id = character())
    
    lTests <- foreach(i = niters, .combine = "c", .multicombine = TRUE) %do% {
        message("  Iteration: ", i)
        
        res_i <- foreach(n = ncells, .combine = "c", .multicombine = TRUE) %dopar% {
            
            out <- tryCatch({
                ## A) Sample Cells
                ids     <- sample(all_cells, n)
                cntable <- as.data.frame(DT_SEGS[sample %in% ids])
                
                ## B) Confidence Calculation
                segDT <- as.data.table(cntable[, c("chromosome","start","end","segVal","sample")])
                
                segDT[, c("confidence","minchange") := {
                    mm <- mean_bin_val_over_range(
                        chr           = as.character(chromosome),
                        s             = as.numeric(start),
                        e             = as.numeric(end),
                        target_sample = as.character(sample),
                        current_ids   = ids,
                        bin_matrix    = MODEL_BIN_MAT,
                        bins_gr       = MODEL_BINS_GR
                    )
                    rawVal <- mm$rawVal
                    oths   <- mm$others
                    
                    if (!length(oths) || is.na(rawVal)) {
                        list(NA_real_, NA_real_)
                    } else {
                        oths <- oths[is.finite(oths)]
                        if (!length(oths)) {
                            list(NA_real_, NA_real_)
                        } else {
                            diffs <- abs(rawVal - round(oths))
                            diffs <- diffs[is.finite(diffs)]
                            if (!length(diffs)) {
                                list(NA_real_, NA_real_)
                            } else {
                                tab   <- table(diffs)
                                confi <- as.numeric(names(sort(tab, decreasing = TRUE))[1])
                                list(confi, min(diffs))
                            }
                        }
                    }
                }, by = seq_len(nrow(segDT))]
                
                confidence <- as.data.frame(segDT)
                
                ## C) Unify & Classify
                unify <- unifySegments.multiple(
                    profiles = cntable, window = WINDOW, bin.size = BIN.SIZE
                )
                if (ncol(unify) > 3) {
                    unify[, 4:ncol(unify)] <- lapply(unify[, 4:ncol(unify)], function(x) round(as.numeric(x)))
                }
                
                events_list <- classEvents.sc(unify, threshold = THRESH)
                events <- data.table::rbindlist(events_list, idcol = "sample_id")
                
                if (nrow(events) == 0) {
                    EMPTY
                } else {
                    
                    required_cols <- c("chromosome", "start", "end", "segVal", "sample_id", "class")
                    if (!all(required_cols %in% colnames(events))) {
                        EMPTY
                    } else {
                        
                        events <- as.data.frame(events[, ..required_cols])
                        colnames(events)[which(colnames(events) == "sample_id")] <- "sample"
                        
                        events <- events[events$segVal != 2, ]
                        unique_ev <- events[events$class == "unique", c("chromosome","start","end","segVal","sample")]
                        
                        ## D) Map Back to Original Segments
                        seg.unique <- data.frame()
                        if (nrow(unique_ev) > 0) {
                            unique.gr <- GRanges(seqnames = unique_ev$chromosome,
                                                 ranges   = IRanges(start = unique_ev$start, end = unique_ev$end))
                            segments.gr <- GRanges(seqnames = cntable$chromosome,
                                                   ranges   = IRanges(start = cntable$start, end = cntable$end))
                            
                            ov <- findOverlaps(unique.gr, segments.gr)
                            
                            seg.unique <- do.call(rbind, lapply(seq_len(nrow(unique_ev)), function(thisUnq) {
                                idxs <- subjectHits(ov)[queryHits(ov) == thisUnq]
                                if(length(idxs) == 0) return(NULL)
                                segs <- cntable[idxs, ]
                                samp <- unique_ev$sample[thisUnq]
                                segs[segs$sample == samp, ]
                            }))
                            
                            if (!is.null(seg.unique) && nrow(seg.unique) > 0) {
                                seg.unique <- seg.unique[!duplicated(seg.unique), ]
                                seg.unique <- seg.unique[order(seg.unique$chromosome, seg.unique$start, seg.unique$end), ]
                            }
                        }
                        
                        ## E) Stats & Null Distribution
                        if (!is.null(seg.unique) && nrow(seg.unique) > 0) {
                            seg.unique <- dplyr::left_join(
                                seg.unique, confidence,
                                by = c("chromosome","start","end","segVal","sample")
                            )
                            
                            length    <- as.numeric(seg.unique$end) - as.numeric(seg.unique$start)
                            thisModel <- MIXMODELS[["segsize"]]
                            exc <- sapply(seq_len(nrow(thisModel)), function(x) {
                                dnorm(length, mean = as.numeric(thisModel[x, "Mean"]), sd = as.numeric(thisModel[x, "SD"]))
                            })
                            if (is.null(dim(exc))) exc <- matrix(exc, nrow = 1)
                            rs <- rowSums(exc)
                            rs[rs == 0 | !is.finite(rs)] <- 1
                            exc <- exc / rs
                            seg.unique$segsize <- apply(exc, 1, function(r) paste0("segsize", which.max(r)))
                            
                            stats <- lapply(seq_len(nrow(seg.unique)), function(ii) {
                                comp <- seg.unique$segsize[ii]
                                minc <- seg.unique$minchange[ii]
                                
                                if (!comp %in% names(lDistributions)) return(c(t.pvalue=1, mean.null=NA, n=0))
                                
                                dist <- lDistributions[[comp]]
                                if (!length(dist) || is.na(minc)) {
                                    c(t.pvalue = 1, mean.null = NA_real_, n = length(dist))
                                } else {
                                    tp <- if (length(dist) > 1) (t.test(minc, dist, var.equal = TRUE, alternative = "greater"))$p.value else 1
                                    c(t.pvalue = tp, mean.null = mean(dist, na.rm = TRUE), n = length(dist))
                                }
                            })
                            stats <- do.call(rbind, stats)
                            seg.unique <- cbind(seg.unique, stats)
                            
                            # 1. Calculate Q-values (BH correction)
                            seg.unique$t.qvalue <- p.adjust(as.numeric(seg.unique$t.pvalue), method = "BH")
                            
                            # 2. Define Significant Events (Match Ground Truth Logic)
                            is_sig <- seg.unique$t.qvalue < 0.01 & seg.unique$minchange >= 0.7
                            
                            # 3. Filter
                            seg.unique <- seg.unique[is_sig, , drop = FALSE]
                            if (nrow(seg.unique) > 0) {
                                seg.unique$id <- paste0(seg.unique$chromosome, ":",
                                                        format(as.numeric(seg.unique$start), scientific = FALSE, trim = TRUE), "-",
                                                        format(as.numeric(seg.unique$end),   scientific = FALSE, trim = TRUE))
                            }
                        }
                        
                        ## F) Keep the FALSE unique events (significant unique, but not truly unique)
                        if (!is.null(seg.unique) && nrow(seg.unique) > 0) {
                            trueIds <- dfTrueUnique$id[dfTrueUnique$model == model]
                            fe <- seg.unique[!(seg.unique$id %in% trueIds), , drop = FALSE]
                            if (nrow(fe) > 0) {
                                data.table(iter = i, ncells = n,
                                           chromosome = as.character(fe$chromosome),
                                           start = as.numeric(fe$start),
                                           end = as.numeric(fe$end),
                                           segVal = as.numeric(fe$segVal),
                                           id = as.character(fe$id))
                            } else EMPTY
                        } else EMPTY
                    }
                }
            }, error = function(e) EMPTY)
            
            list(out)
        }
        res_i
    }
    
    false_events <- data.table::rbindlist(lTests, fill = TRUE)
    
    ## -----------------------------------------------------------------------
    ## 4. COUNT CELLS IN FINAL POPULATION CARRYING EACH FALSE EVENT
    ##    (compute once per distinct coords+CN, then map back)
    ## -----------------------------------------------------------------------
    if (nrow(false_events) > 0) {
        key <- paste(false_events$chromosome, false_events$start, false_events$end,
                     round(false_events$segVal), sep = "_")
        uk  <- !duplicated(key)
        ue  <- false_events[uk]
        message(sprintf("  Counting cells for %d distinct false events (%d total detections)...",
                        nrow(ue), nrow(false_events)))
        ucnt <- mapply(count_cells_in_pop, ue$chromosome, ue$start, ue$end, ue$segVal,
                       MoreArgs = list(ref_gr = REF_GR))
        false_events$n_cells_in_pop <- as.integer(ucnt[match(key, key[uk])])
    } else {
        false_events$n_cells_in_pop <- integer()
    }
    
    false_events$model  <- model
    false_events$ntotal <- n_total
    
    saveRDS(false_events,
            file.path(OUTPUT_DIR, paste0("stats_false_unique_cells_in_pop_", model, ".rds")))
    message(sprintf("  Saved %d false-event rows for %s.", nrow(false_events), model))
}
