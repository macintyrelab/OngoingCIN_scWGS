################################################################################
## This script processes the output from 04a_MixtureSimulations.R and applies
## the signature mapping from 02_CINSimulations_Testing_Mapping.R to evaluate
## how well signature activities are recovered from pairwise mixtures.
##
## For each pair of signatures, it plots CXA/(CXA+CXB) (mapped activities)
## against the event proportion (nA:nB from 1:99 to 99:1), comparing mapped
## activities to the real (expected) activities.
################################################################################

### LIBRARIES
rm(list = ls(all = TRUE))

library(tidyverse)
library(ggpubr)
library(ggthemes)
library(this.path)

## PLOTS THEME
theme_set(theme_tufte(base_size = 6, base_family = "ArialMT"))
theme_update(text = element_text(size = 6),
             axis.text = element_text(size = 6),
             axis.line = element_line(size = 0.5),
             axis.ticks = element_line(size = 0.5),
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
## and systematic offset/scale shift, so it drops when mapped/threshold lines
## diverge from the ground-truth line even if the two are linearly related.
ccc <- function(x, y) {
  ok <- complete.cases(x, y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 2) return(NA_real_)
  2 * cov(x, y) / (var(x) + var(y) + (mean(x) - mean(y))^2)
}

## PARAMETERS
# signatures <- c("CX1","CX2","CX3","CX4","CX5","CX6","CX8","CX11","CX13")
# signatures <- c("CX1", "CX2", "CX13")
signatures <- rownames(DEFS)
feature_cols <- colnames(DEFS)
seg_cols <- feature_cols[1:22]
chp_cols <- feature_cols[23:32]
bpw_cols <- feature_cols[33:41]

## Normalize signature definitions per feature group
## Handle zero-sum rows by assigning uniform weights
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
#### MAP EVENTS TO SIGNATURES ##################################################

# Per-event mapping (same dot-product scoring as script 1, section 6).
# Each event row is one-hot encoded across feature bins; score against each
# active signature, normalize, then aggregate per sample for activities.

# Extract feature matrix as numeric
feature_data <- as.matrix(sample_matrix[, feature_cols])
storage.mode(feature_data) <- "double"

# Event metadata
event_ids <- as.character(sample_matrix$sample_id)

# Match each event to its pair info via sample_id
event_pairs <- truth_df[match(event_ids, truth_df$sample_id), ]

# 1. Score all events against all signatures at once: (n_events x 41) %*% (41 x n_sigs)
raw_scores <- feature_data %*% t(scdefs)

# 2. build per-event activity matrix (events x signatures), binary mask
act_mat_bin <- (raw_scores > 0) * 1

# 3. apply mask
raw_masked <- raw_scores * act_mat_bin

# 4. row-wise normalize to probabilities, handle zero rows
row_sums <- rowSums(raw_masked, na.rm = TRUE)
zero_rows <- which(row_sums == 0 | is.na(row_sums))
row_sums[zero_rows] <- 1    # avoid division by zero (or set to NA/other fallback)

all_scores <- raw_masked / row_sums

# Hard assignment: each event to the max-scoring signature
mapped_sig <- apply(all_scores, 1, function(x) {
  if (all(is.na(x))) return(NA)
  sigs <- colnames(all_scores)[which(x == max(x, na.rm = TRUE))]
  if (length(sigs) > 1) return(NA) # tie, no assignment
  else return(sigs)
})

# Aggregate per sample: soft scores + hard assignment counts
base_df <- data.frame(
  sample_id = event_ids,
  sigA = event_pairs$sigA,
  sigB = event_pairs$sigB,
  mapped_sig = mapped_sig,
  stringsAsFactors = FALSE
)

# counts per sample per signature (table excludes NA by default)
counts_table <- as.data.frame(table(base_df$sample_id, base_df$mapped_sig), stringsAsFactors = FALSE)
colnames(counts_table) <- c("sample_id", "mapped_sig", "Freq")
counts_table <- counts_table[!is.na(counts_table$mapped_sig) & counts_table$mapped_sig != "", , drop = FALSE]

counts_wide <- counts_table %>%
  pivot_wider(names_from = mapped_sig, values_from = Freq, values_fill = 0) %>%
  rename_with(.fn = function(x) paste0("n_to_", x), .cols = -sample_id)

# Summary per sample and attach per-signature counts
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

# Add signature columns that received zero counts (if any)
for (sig in signatures) {
    col_name <- paste0("n_to_", sig)
    if (!col_name %in% colnames(sample_results)) {
        sample_results[[col_name]] <- 0
    }
}

# Thresholded activity: a signature is only valid if >= 5% of events AND >= 3 events. Otherwise, that event does not count.
# Applied to all signatures in the signatures list; n_to_{sig}_thresh columns are generated for each.
thresh_cols <- setNames(
  lapply(signatures, function(sig) paste0("n_to_", sig, "_thresh")),
  signatures
)

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

# Add n_to_{sig}_thresh columns for all signatures
for (sig in signatures) {
  col_thresh <- paste0("n_to_", sig, "_thresh")
  col_raw    <- paste0("n_to_", sig)
  sample_results[[col_thresh]] <- ifelse(
    mapply(function(vs) sig %in% vs, sample_results$valid_sigs),
    sample_results[[col_raw]],
    0L
  )
}


# Join with truth for real activities
mapped_activities <- sample_results %>%
  inner_join(truth_df, by = c("sample_id", "sigA", "sigB")) %>%
  mutate(
    event_proportion  = n_events_A,
    n_sig_events      = n_events_A + n_events_B,   # true signature-derived events (excludes random)
    real_activity     = n_events_A / n_sig_events,
    mapped_activity   = ifelse(n_sig_events > 0, n_to_A / n_total, NA_real_)
  )

# thresh_activity
mapped_activities$thresh_activity <- mapply(
  function(sA, denom, i) {
    col <- paste0("n_to_", sA, "_thresh")
    if (denom > 0) mapped_activities[[col]][i] / denom else NA_real_
  },
  mapped_activities$sigA,
  mapped_activities$n_total_thresh,
  seq_len(nrow(mapped_activities))
)

mapped_activities <- mapped_activities %>%
  arrange(sigA, sigB, n_events_A)

################################################################################
#### PLOT RESULTS ##############################################################


sig_pairs <- combn(c("CX1", "CX2", "CX13"), 2, simplify = FALSE)

plots <- lapply(sig_pairs, function(pair) {
  sigA <- pair[1]
  sigB <- pair[2]

  df <- mapped_activities %>%
    filter(sigA == !!sigA, sigB == !!sigB)

  # Mean and SD per event proportion across replicates
  df_summary <- df %>%
    group_by(event_proportion) %>%
    summarise(
      mean_mapped  = mean(mapped_activity),
      sd_mapped    = sd(mapped_activity),
      mean_thresh  = mean(thresh_activity, na.rm = TRUE),
      sd_thresh    = sd(thresh_activity, na.rm = TRUE),
      real_activity = first(real_activity),
      .groups = "drop"
    )

  # CCC, RMSE and mean bias of mapped and thresholded activity vs real
  # activity (computed on per-replicate values, not just means)
  ccc_mapped  <- ccc(df$mapped_activity, df$real_activity)
  ccc_thresh  <- ccc(df$thresh_activity, df$real_activity)
  rmse_mapped <- sqrt(mean((df$mapped_activity - df$real_activity)^2,
                           na.rm = TRUE))
  rmse_thresh <- sqrt(mean((df$thresh_activity - df$real_activity)^2,
                           na.rm = TRUE))
  bias_mapped <- mean(df$mapped_activity - df$real_activity, na.rm = TRUE)
  bias_thresh <- mean(df$thresh_activity - df$real_activity, na.rm = TRUE)
  metrics_label <- sprintf(
    paste0("Mapped:   CCC = %.3f, RMSE = %.3f, bias = %+.3f\n",
           "Threshold: CCC = %.3f, RMSE = %.3f, bias = %+.3f"),
    ccc_mapped, rmse_mapped, bias_mapped,
    ccc_thresh, rmse_thresh, bias_thresh)

  ggplot() +
    # Shaded region between real and mapped average lines
    geom_ribbon(data = df_summary,
                aes(x = event_proportion,
                    ymin = pmin(mean_thresh, mean_mapped),
                    ymax = pmax(mean_thresh, mean_mapped)),
                fill = "grey70", alpha = 0.4) +
    # SD ribbon around mapped mean
    geom_ribbon(data = df_summary,
                aes(x = event_proportion,
                    ymin = mean_mapped - sd_mapped,
                    ymax = mean_mapped + sd_mapped),
                fill = "steelblue", alpha = 0.2) +
    # SD ribbon around thresholded mean
    geom_ribbon(data = df_summary,
                aes(x = event_proportion,
                    ymin = mean_thresh - sd_thresh,
                    ymax = mean_thresh + sd_thresh),
                fill = "firebrick", alpha = 0.15) +
    # Scatter points - unthresholded
    geom_point(data = df,
               aes(x = event_proportion, y = mapped_activity),
               color = "steelblue", alpha = 0.3, size = 0.3) +
    # Scatter points - thresholded
    geom_point(data = df,
               aes(x = event_proportion, y = thresh_activity),
               color = "firebrick", alpha = 0.3, size = 0.3) +
    # Lines with color aesthetic for legend
    geom_line(data = df_summary,
              aes(x = event_proportion, y = real_activity, color = "Real activity"),
              linetype = "dashed", linewidth = 0.4) +
    geom_line(data = df_summary,
              aes(x = event_proportion, y = mean_mapped, color = "Mapped"),
              linewidth = 0.4) +
    geom_line(data = df_summary,
              aes(x = event_proportion, y = mean_thresh, color = "Mapped (applying thresholds)"),
              linewidth = 0.4) +
    scale_color_manual(
      values = c("Real activity" = "black",
                 "Mapped" = "steelblue",
                 "Mapped (applying thresholds)" = "firebrick"),
      breaks = c("Real proportion", "Mapped", "Mapped (applying thresholds)")
    ) +
    labs(
      x = paste0("Event proportion (", sigA, ":", sigB, ")"),
      y = bquote(N[.(sigA)] / N[TOTAL]),
      title = paste0(sigA, " vs ", sigB),
      color = NULL
    ) +
    scale_x_continuous(breaks = c(1, 25, 50, 75, 99),
                       labels = c("1:99", "25:75", "50:50", "75:25", "99:1")) +
    scale_y_continuous(limits = c(0, 1)) +
    annotate("text", x = 99, y = 0.05, hjust = 1, vjust = 0,
             label = metrics_label, size = 1.8, family = "ArialMT") +
    theme(legend.position = "bottom",
          legend.key.size = unit(3, "mm"))
})

# Combine and save
p <- ggarrange(plotlist = plots, ncol = length(sig_pairs), nrow = 1,
               common.legend = TRUE, legend = "bottom")

pdf(file = file.path(PLOTS_DIR, "3_MixtureSimulations_ActivityRecovery.pdf"),
    height = 80 / 25.4, width = length(sig_pairs) * 80 / 25.4)
print(p)
dev.off()

cat("Done. Plots saved to:", PLOTS_DIR, "\n")

################################################################################
#### PLOT SIGNATURE DISTRIBUTIONS ##############################################

feature_groups <- list(
  Segment   = seg_cols,
  Changepoint = chp_cols,
  Breakpoint  = bpw_cols
)

# Build long-format normalised distribution table
sig_dists <- lapply(signatures, function(sig) {
  do.call(rbind, lapply(names(feature_groups), function(grp) {
    cols  <- feature_groups[[grp]]
    probs <- safe_normalize(matrix(as.numeric(DEFS[sig, cols]), nrow = 1))[1, ]
    data.frame(
      signature = sig,
      group     = grp,
      feature   = cols,
      feature_i = seq_along(cols),
      prob      = probs,
      stringsAsFactors = FALSE
    )
  }))
})
sig_dists_df <- do.call(rbind, sig_dists)
sig_dists_df$signature <- factor(sig_dists_df$signature, levels = signatures)
sig_dists_df$group     <- factor(sig_dists_df$group, levels = names(feature_groups))

p_dists <- ggplot(sig_dists_df, aes(x = feature_i, y = prob)) +
  geom_col(fill = "steelblue", width = 0.8) +
  facet_grid(signature ~ group, scales = "free_x", space = "free_x") +
  labs(x = "Feature bin", y = "Probability", title = "Signature feature distributions") +
  theme(strip.text = element_text(size = 6),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# pdf(file = file.path(PLOTS_DIR, "3_MixtureSimulations_SignatureDistributions.pdf"),
#     height = length(signatures) * 40 / 25.4,
#     width  = 160 / 25.4)
# print(p_dists)
# dev.off()

cat("Signature distributions saved to:", PLOTS_DIR, "\n")
