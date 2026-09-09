################################################################################
## Faster script: number of cells vs unique detection (Precomputed List Input)
## Started: 24/10/25   |   Fast rewrite: 14/11/25   |   List-Input: 19/11/25
## Fixed: 23/11/25 (Added Sample Intersection Logic to fix TN5 crash)
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
  library(ggplot2)
  library(ggthemes)
  library(ggpubr)
  library(RColorBrewer)
  library(lemon)
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
PLOTS_DIR  <- file.path(MISC_DIR,"plots")
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)

## ---------- Theme ----------
theme_set(theme_tufte(base_size = 6, base_family = "ArialMT"))
theme_update(text = element_text(size = 6),
             axis.text = element_text(size = 6),
             axis.line = element_line(size = 0.5),
             axis.ticks = element_line(size = 0.5),
             axis.ticks.length = unit(.1, "cm"))

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

## ============================================================================
##                          MAIN SIMULATIONS
## ============================================================================

niters <- 1:20

MODELS_TO_RUN <- intersect(names(BIN.TABLE), names(SEG.TABLE))
# Filter for specific models if desired
# MODELS_TO_RUN <- MODELS_TO_RUN[!(MODELS_TO_RUN %in% c("TN1","TN2","MDAMB231_popp31"))]
# MODELS_TO_RUN <- "TN1"
# MODELS_TO_RUN <- "TN2"
# MODELS_TO_RUN <- "MDAMB231_popp31"

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
  ## 3. RUN SIMULATION LOOPS
  ## -----------------------------------------------------------------------
  lTests <- foreach(i = niters, .combine = "rbind") %do% {
    message("  Iteration: ", i)

    res_i <- foreach(n = ncells, .combine = "rbind") %dopar% {

      tryCatch({
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
          c(iter = i, ncells = n, nunique = 0, nreal = 0, nfalse = 0)
        } else {

          required_cols <- c("chromosome", "start", "end", "segVal", "sample_id", "class")
          if (!all(required_cols %in% colnames(events))) {
            return(c(iter = i, ncells = n, nunique = 0, nreal = 0, nfalse = 0))
          }

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
            # Ground Truth used SIGN.THRESH = 0.01. Ensure simulation uses the same.
            is_sig <- seg.unique$t.qvalue < 0.01 & seg.unique$minchange >= 0.7

            # 3. Filter
            seg.unique <- seg.unique[is_sig, , drop = FALSE]
            if (nrow(seg.unique) > 0) {
              seg.unique$id <- paste0(seg.unique$chromosome, ":",
                                      format(as.numeric(seg.unique$start), scientific = FALSE, trim = TRUE), "-",
                                      format(as.numeric(seg.unique$end),   scientific = FALSE, trim = TRUE))
            }
          }

          nunique <- if(!is.null(seg.unique)) nrow(seg.unique) else 0
          nreal   <- if (nunique > 0) sum(seg.unique$id %in% dfTrueUnique$id[dfTrueUnique$model == model]) else 0
          nfalse  <- nunique - nreal

          c(iter = i, ncells = n, nunique = nunique, nreal = nreal, nfalse = nfalse)
        }
      }, error = function(e) {
        # Silent error handling (or log if preferred)
        c(iter = i, ncells = n, nunique = -1, nreal = -1, nfalse = -1)
      })
    }
    res_i
  }

  dfTest <- as.data.frame(lTests)
  dfTest$iter    <- as.integer(dfTest$iter)
  dfTest$ncells  <- as.integer(dfTest$ncells)
  dfTest$nunique <- as.integer(dfTest$nunique)
  dfTest$nreal   <- as.integer(dfTest$nreal)
  dfTest$nfalse  <- as.integer(dfTest$nfalse)
  dfTest$model   <- model

  saveRDS(dfTest, file.path(OUTPUT_DIR, paste0("stats_ncells_unique_fp_detection_", model, ".rds")))
}
