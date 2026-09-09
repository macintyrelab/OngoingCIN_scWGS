################################################################################
## QC analyses for the mixture simulations produced in 04a_MixtureSimulations.R
## and mapped in 04b_MixtureSimulations_Mapping.R.
##
## Goals supported:
##   I.  Sample-level activity recovery quality: per-pair CCC, R^2, RMSE and
##       bias for both raw mapped activities and thresholded activities (vs
##       real), consolidated in qc_recovery_table.csv. CCC is also annotated
##       on Fig 3 by 03_MixtureSimulations_Mapping.R.
##   II. Per-event detection rate at the boundaries (Fig 5): mean per-sample
##       fraction of events whose hard-assigned signature matches the
##       ground-truth source, comparing raw mapping with the thresholded
##       cohort, with SE across replicates.
################################################################################

### LIBRARIES
rm(list = ls(all = TRUE))

library(tidyverse)
library(ggpubr)
library(ggthemes)

## PLOTS THEME
theme_set(theme_tufte(base_size = 6, base_family = "ArialMT"))
theme_update(text = element_text(size = 6),
             axis.text = element_text(size = 6),
             axis.line = element_line(linewidth = 0.5),
             axis.ticks = element_line(linewidth = 0.5),
             axis.ticks.length = unit(.1, "cm"))

## PATHS
BASE <- dirname(dirname(this.path()))
SIMS_DIR <- dirname(this.path())
SC_DIR <- file.path(BASE,"02_Derive_New_scCINSignatures")
OUTPUT_DIR <- file.path(SIMDIR, "outputs", "mixture_simulations")
PLOTS_DIR <- file.path(SIMDIR, "plots", "mixture_simulations")

## LOAD DATA
DEFS <- readRDS(file = paste0(SC_DIR, "/outputs/2_New_scCINSignatures_Definitions.rds"))
sample_matrix <- readRDS(file = file.path(OUTPUT_DIR, "mixture_samples.rds"))
truth_df <- readRDS(file = file.path(OUTPUT_DIR, "mixture_summary.rds"))

## CCC (Lin's concordance correlation coefficient): bounded agreement metric
## in [-1, 1]; 1 = perfect agreement on y = x. Penalises both correlation loss
## and systematic offset/scale shift.
ccc <- function(x, y) {
  ok <- complete.cases(x, y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 2) return(NA_real_)
  2 * cov(x, y) / (var(x) + var(y) + (mean(x) - mean(y))^2)
}

## PARAMETERS
signatures <- rownames(DEFS)
feature_cols <- colnames(DEFS)
seg_cols <- feature_cols[1:22]
chp_cols <- feature_cols[23:32]
bpw_cols <- feature_cols[33:41]

## Normalize signature definitions per feature group (zero-sum -> uniform)
safe_normalize <- function(mat) {
  rs <- rowSums(mat)
  for (i in seq_len(nrow(mat))) {
    if (rs[i] == 0) mat[i, ] <- 1 / ncol(mat)
    else mat[i, ] <- mat[i, ] / rs[i]
  }
  mat
}

defs_seg <- safe_normalize(DEFS[signatures, seg_cols])
defs_chp <- safe_normalize(DEFS[signatures, chp_cols])
defs_bpw <- safe_normalize(DEFS[signatures, bpw_cols])
scdefs <- as.matrix(cbind(defs_seg, defs_chp, defs_bpw))

################################################################################
#### MAP EVENTS (replicates the mapping in 04b_MixtureSimulations_Mapping.R) ####

feature_data <- as.matrix(sample_matrix[, feature_cols])
storage.mode(feature_data) <- "double"

event_ids <- as.character(sample_matrix$sample_id)
event_pairs <- truth_df[match(event_ids, truth_df$sample_id), ]

raw_scores <- feature_data %*% t(scdefs)
act_mat_bin <- (raw_scores > 0) * 1
raw_masked <- raw_scores * act_mat_bin

row_sums <- rowSums(raw_masked, na.rm = TRUE)
zero_rows <- which(row_sums == 0 | is.na(row_sums))
row_sums[zero_rows] <- 1
all_scores <- raw_masked / row_sums

# Hard assignment to argmax (NA on tie)
mapped_sig <- apply(all_scores, 1, function(x) {
  if (all(is.na(x))) return(NA_character_)
  sigs <- colnames(all_scores)[which(x == max(x, na.rm = TRUE))]
  if (length(sigs) > 1) return(NA_character_)
  return(sigs)
})

# Ground-truth source signature for each event (sigA / sigB / "random")
source_sig <- as.character(sample_matrix$source_sig)

################################################################################
#### 1) PER-SAMPLE COUNTS, VALID SIGNATURES, AND ACTIVITIES ####################

# Reproduce the per-sample tables produced by script 03 so the QC numbers and
# the figure annotations stay in sync.

base_df <- data.frame(
  sample_id = event_ids,
  sigA = event_pairs$sigA,
  sigB = event_pairs$sigB,
  mapped_sig = mapped_sig,
  stringsAsFactors = FALSE
)

counts_table <- as.data.frame(table(base_df$sample_id, base_df$mapped_sig),
                              stringsAsFactors = FALSE)
colnames(counts_table) <- c("sample_id", "mapped_sig", "Freq")
counts_table <- counts_table[!is.na(counts_table$mapped_sig) &
                               counts_table$mapped_sig != "", , drop = FALSE]

counts_wide <- counts_table %>%
  pivot_wider(names_from = mapped_sig, values_from = Freq, values_fill = 0) %>%
  rename_with(.fn = function(x) paste0("n_to_", x), .cols = -sample_id)

sample_results <- base_df %>%
  group_by(sample_id, sigA, sigB) %>%
  summarise(
    n_to_A = sum(mapped_sig == sigA, na.rm = TRUE),
    n_to_B = sum(mapped_sig == sigB, na.rm = TRUE),
    n_total = n(),
    .groups = "drop"
  ) %>%
  left_join(counts_wide, by = "sample_id") %>%
  mutate(across(starts_with("n_to_"), ~ replace_na(.x, 0)))

for (sig in signatures) {
  col_name <- paste0("n_to_", sig)
  if (!col_name %in% colnames(sample_results)) sample_results[[col_name]] <- 0
}

# 5% / >=3-event threshold (AND, matching script 03 implementation)
sample_results <- sample_results %>%
  rowwise() %>%
  mutate(
    valid_sigs = list(
      signatures[sapply(signatures, function(sig) {
        n_to_sig <- get(paste0("n_to_", sig))
        (n_to_sig / n_total >= 0.05) & (n_to_sig >= 3)
      })]
    ),
    n_total_thresh = sum(sapply(signatures, function(sig) {
      if (sig %in% valid_sigs) get(paste0("n_to_", sig)) else 0
    }))
  ) %>%
  ungroup()

# Join with truth and compute raw and thresholded activities
mapped_activities <- sample_results %>%
  inner_join(truth_df, by = c("sample_id", "sigA", "sigB")) %>%
  rowwise() %>%
  mutate(
    real_activity_A = n_events_A / (n_events_A + n_events_B),
    real_activity_B = n_events_B / (n_events_A + n_events_B),
    mapped_activity_A = n_to_A / n_total,
    mapped_activity_B = n_to_B / n_total,
    thresh_activity_A = ifelse(n_total_thresh > 0,
                               ifelse(sigA %in% valid_sigs,
                                      n_to_A / n_total_thresh, 0),
                               NA_real_),
    thresh_activity_B = ifelse(n_total_thresh > 0,
                               ifelse(sigB %in% valid_sigs,
                                      n_to_B / n_total_thresh, 0),
                               NA_real_),
    minor_sig = ifelse(n_events_A <= n_events_B, sigA, sigB),
    minor_n   = pmin(n_events_A, n_events_B),
    minor_frac = minor_n / n_total,
    n_to_minor = ifelse(minor_sig == sigA, n_to_A, n_to_B),
    raw_minor_activity    = n_to_minor / n_total,
    thresh_minor_activity = ifelse(n_total_thresh > 0 &
                                     minor_sig %in% valid_sigs,
                                   n_to_minor / n_total_thresh, NA_real_),
    pair = paste0(sigA, " vs ", sigB)
  ) %>%
  ungroup()

################################################################################
#### 2) RECOVERY-QUALITY TABLE (CCC, R^2, RMSE, BIAS) ##########################

recovery_table <- mapped_activities %>%
  group_by(sigA, sigB) %>%
  summarise(
    n_samples       = n(),
    CCC_raw_A       = ccc(mapped_activity_A, real_activity_A),
    CCC_raw_B       = ccc(mapped_activity_B, real_activity_B),
    CCC_thresh_A    = ccc(thresh_activity_A, real_activity_A),
    CCC_thresh_B    = ccc(thresh_activity_B, real_activity_B),
    R2_raw_A        = cor(mapped_activity_A, real_activity_A,
                          use = "complete.obs")^2,
    R2_raw_B        = cor(mapped_activity_B, real_activity_B,
                          use = "complete.obs")^2,
    R2_thresh_A     = cor(thresh_activity_A, real_activity_A,
                          use = "complete.obs")^2,
    R2_thresh_B     = cor(thresh_activity_B, real_activity_B,
                          use = "complete.obs")^2,
    rmse_raw_A      = sqrt(mean((mapped_activity_A - real_activity_A)^2)),
    rmse_raw_B      = sqrt(mean((mapped_activity_B - real_activity_B)^2)),
    rmse_thresh_A   = sqrt(mean((thresh_activity_A - real_activity_A)^2,
                                 na.rm = TRUE)),
    rmse_thresh_B   = sqrt(mean((thresh_activity_B - real_activity_B)^2,
                                 na.rm = TRUE)),
    bias_raw_A      = mean(mapped_activity_A - real_activity_A),
    bias_raw_B      = mean(mapped_activity_B - real_activity_B),
    bias_thresh_A   = mean(thresh_activity_A - real_activity_A,
                            na.rm = TRUE),
    bias_thresh_B   = mean(thresh_activity_B - real_activity_B,
                            na.rm = TRUE),
    .groups = "drop"
  )

write.csv(recovery_table,
          file = file.path(OUTPUT_DIR, "qc_recovery_table.csv"),
          row.names = FALSE)

################################################################################
#### 3) PER-EVENT DETECTION RATE: RAW vs THRESHOLDED (Fig 5) ###################

# True per-event detection: an event is "detected" if its hard-assigned
# signature matches its ground-truth source. We then average the per-event
# correctness across all signature-derived events of each sample, and
# aggregate the resulting per-sample rates across the 20 replicates available
# at each (pair, minor_n) bin.
#
# Two flavours:
#   - Raw:         mapped_sig == source_sig (any signature-derived event).
#   - Thresholded: same, but events whose mapped signature fails the 5% AND
#                  >=3-event filter for the host sample are counted as
#                  not-detected (denominator unchanged), so the curve also
#                  reflects information lost to thresholding at the boundaries.

ev_df <- data.frame(
  sample_id  = event_ids,
  sigA       = event_pairs$sigA,
  sigB       = event_pairs$sigB,
  source_sig = source_sig,
  mapped_sig = mapped_sig,
  stringsAsFactors = FALSE
) %>%
  filter(source_sig != "random") %>%
  left_join(sample_results %>% select(sample_id, valid_sigs),
            by = "sample_id")

# Vectorised per-event correctness (avoids slow rowwise on ~290k events)
ev_df$raw_correct <- !is.na(ev_df$mapped_sig) &
  ev_df$mapped_sig == ev_df$source_sig
ev_df$thresh_correct <- ev_df$raw_correct &
  vapply(seq_len(nrow(ev_df)),
         function(i) ev_df$mapped_sig[i] %in% ev_df$valid_sigs[[i]],
         logical(1))

# Per-sample detection rate
det_per_sample <- ev_df %>%
  group_by(sample_id, sigA, sigB) %>%
  summarise(
    rate_raw    = mean(raw_correct),
    rate_thresh = mean(thresh_correct),
    .groups = "drop"
  ) %>%
  inner_join(truth_df, by = c("sample_id", "sigA", "sigB")) %>%
  mutate(
    pair       = paste0(sigA, " vs ", sigB),
    minor_n    = pmin(n_events_A, n_events_B),
    minor_frac = minor_n / 105
  )

# Aggregate across replicates per (pair, minor_n)
detect_summary <- det_per_sample %>%
  group_by(pair, minor_n, minor_frac) %>%
  summarise(
    n_rep       = n(),
    mean_raw    = mean(rate_raw),
    sd_raw      = sd(rate_raw),
    se_raw      = sd(rate_raw) / sqrt(n()),
    mean_thresh = mean(rate_thresh),
    sd_thresh   = sd(rate_thresh),
    se_thresh   = sd(rate_thresh) / sqrt(n()),
    .groups = "drop"
  )

write.csv(detect_summary,
          file = file.path(OUTPUT_DIR, "qc_detection_rate.csv"),
          row.names = FALSE)

# Long format for plotting (pair x method)
detect_long <- bind_rows(
  detect_summary %>%
    transmute(pair, minor_frac, method = "Raw",
              rate = mean_raw, se = se_raw),
  detect_summary %>%
    transmute(pair, minor_frac, method = "Thresholded",
              rate = mean_thresh, se = se_thresh)
) %>%
  mutate(method = factor(method, levels = c("Raw", "Thresholded")))

p_detection <- ggplot(detect_long,
                      aes(x = minor_frac, y = rate,
                          colour = pair, fill = pair, linetype = method)) +
  geom_vline(xintercept = 0.05, linetype = "dotted",
             colour = "grey50", linewidth = 0.3) +
  geom_ribbon(aes(ymin = pmax(rate - se, 0),
                  ymax = pmin(rate + se, 1),
                  group = interaction(pair, method)),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.4) +
  scale_x_continuous(limits = c(0, 0.5),
                     breaks = c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5),
                     labels = c("0", "5%", "10%", "20%", "30%",
                                 "40%", "50%")) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Fraction of unique events from the minor signature",
       y = "Per-event detection rate (mean per sample)",
       colour = "Pair", fill = "Pair", linetype = "Mapping",
       title = "Per-event signature detection: raw vs thresholded") +
  theme(legend.position = "bottom",
        legend.key.size = unit(3, "mm"))

pdf(file = file.path(PLOTS_DIR, "5_MixtureSimulations_DetectionRate.pdf"),
    height = 80 / 25.4, width = 110 / 25.4)
print(p_detection)
dev.off()

################################################################################
#### CONSOLE SUMMARY ###########################################################

cat("\n=== Recovery quality (R^2, RMSE, bias) per pair, raw vs thresholded ===\n")
print(recovery_table)

cat("\nFigures:\n",
    file.path(PLOTS_DIR, "5_MixtureSimulations_DetectionRate.pdf"), "\n")
cat("Tables:\n",
    file.path(OUTPUT_DIR, "qc_recovery_table.csv"), "\n",
    file.path(OUTPUT_DIR, "qc_detection_rate.csv"), "\n")
