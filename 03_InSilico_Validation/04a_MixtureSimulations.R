################################################################################
## This script simulates samples with alterations caused by two signatures. The
## contribution of each signature to the final pool of alterations is modified.
## We also add a low number of alterations caused by a random process.
##
################################################################################


library(dplyr)
library(this.path)

## PATHS
BASE <- dirname(dirname(this.path()))
SIMS_DIR <- dirname(this.path())
SC_DIR <- file.path(BASE,"02_Derive_New_scCINSignatures")
OUTPUT_DIR <- file.path(SIMDIR, "outputs", "mixture_simulations")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
PLOTS_DIR <- file.path(SIMDIR, "plots", "mixture_simulations")
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)


## LOAD DATA
DEFS <- readRDS(file = paste0(SC_DIR, "/outputs/2_New_scCINSignatures_Definitions.rds"))

## PARAMETERS
signatures <- c("CX1", "CX2", "CX13")
n_total <- 105
n_random <- 5
n_events <- n_total - n_random       # Total events per sample
n_replicates <- 10
seed <- 42

set.seed(seed)

# All pairwise combinations
sig_pairs <- combn(signatures, 2, simplify = FALSE)

# Feature column names from DEFS
feature_cols <- colnames(DEFS)
seg_cols <- feature_cols[1:22]
chp_cols <- feature_cols[23:32]
bpw_cols <- feature_cols[33:41]

# Normalize signature definitions per feature group
get_probs <- function(sig, cols) {
  v <- as.numeric(DEFS[sig, cols])
  if (sum(v) == 0) rep(1 / length(cols), length(cols))
  else v / sum(v)
}

# Sample a single event from a signature, returning a one-hot row (41 columns)
sample_one_event <- function(sig) {
  row <- rep(0L, length(feature_cols))

  row[sample(1:22,  1, prob = get_probs(sig, seg_cols))] <- 1L
  row[22 + sample(1:10, 1, prob = get_probs(sig, chp_cols))] <- 1L
  row[32 + sample(1:9,  1, prob = get_probs(sig, bpw_cols))] <- 1L

  row
}

# Sample a single random event with uniform component distributions
sample_random_event <- function() {
  row <- rep(0L, length(feature_cols))

  row[sample(1:22,  1)] <- 1L
  row[22 + sample(1:10, 1)] <- 1L
  row[32 + sample(1:9,  1)] <- 1L

  row
}

# Run simulation
all_events <- list()
all_truth <- list()

for (rep_id in 1:n_replicates) {
  for (pair in sig_pairs) {
    sigA <- pair[1]
    sigB <- pair[2]

    for (nA in (n_events-1):1) {
      nB <- n_events - nA
      sample_id <- paste0("rep", rep_id, "_", sigA, "_", sigB, "_", nA, "v", nB)

      # Generate events from each signature + random events
      rand_events <- do.call(rbind, replicate(n_random, sample_random_event(), simplify = FALSE))
      events <- rbind(
        do.call(rbind, replicate(nA, sample_one_event(sigA), simplify = FALSE)),
        do.call(rbind, replicate(nB, sample_one_event(sigB), simplify = FALSE)),
        rand_events
      )
      events_df <- as.data.frame(events)
      colnames(events_df) <- feature_cols
      events_df$sample_id <- sample_id
      events_df$source_sig <- c(rep(sigA, nA), rep(sigB, nB), rep("random", n_random))

      all_events[[sample_id]] <- events_df

      all_truth[[sample_id]] <- data.frame(
        sample_id = sample_id,
        replicate = rep_id,
        sigA = sigA,
        sigB = sigB,
        n_events_A = nA,
        n_events_B = nB,
        stringsAsFactors = FALSE
      )
    }
  }
}

# Assemble outputs
events_df <- do.call(rbind, all_events)
rownames(events_df) <- NULL

truth_df <- do.call(rbind, all_truth)
rownames(truth_df) <- NULL

# Save
saveRDS(events_df, file = file.path(OUTPUT_DIR, "mixture_samples.rds"))
saveRDS(truth_df, file = file.path(OUTPUT_DIR, "mixture_summary.rds"))

cat("Done.", length(unique(events_df$sample_id)), "samples,",
    nrow(events_df), "total events generated.\n")
cat("Output saved to:", OUTPUT_DIR, "\n")
