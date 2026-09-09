# =============================================================================
# Reanalyse the simulation results using a HORIZON-ANCHORED generation-based
# ground-truth threshold:
#
#   t_thresh(n, k) = (log(n) - k * log(2)) / log(N_final)
#
# Derivation. A sample of size n has a temporal horizon at the relative time
# when the population was the size of the sample:
#
#   t_horizon(n) = log(n) / log(N_final)
#
# An event is labelled "ongoing" if it arose within the last k cell-division
# generations BEFORE the horizon. One generation in relative-time units is
# log(2) / log(N_final). Subtracting gives the formula above.
#
# The horizon-anchored threshold reduces to the original generation-based
# threshold (1 - k*log(2)/log(N)) at full sampling (n = N), so the two
# previous criteria are unified into one expression that depends naturally
# on sample size and increases monotonically with n.
#
# Inputs
#   ongoing_cin_simulation_results.rds
#
# Outputs
#   threshold_diagram.pdf            visualises the horizon-anchored threshold
#   scatter_freq_vs_trel.pdf         raw mutant frequency vs t_mut
#   singleton_prob_curves.pdf  P(singleton | t_mut), per-panel threshold
#   classification_metrics.pdf per-s sens / spec / PPV / F1 / accuracy
#   pr_tradeoff.pdf            precision-recall trajectory through s
#   dfe_scenarios.pdf                the five DFE shapes
#   marginal_metrics.pdf       DFE-marginalised metrics vs sample frac
#   extinction_rate.pdf              mutation extinction rate vs t_mut, per s
#   per_s_metrics.csv          per-(s, sample fraction) numeric table
#   marginal_summary.csv       marginal numeric table
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(scales)
  library(viridis)
  library(this.path)
})

# PATHS
BASE=dirname(dirname(this.path()))
SIMS_DIR=dirname(this.path())
OUTPUT_DIR=file.path(SIMS_DIR,"outputs/selection_simulations")
PLOTS_DIR=file.path(SIMS_DIR,"plots/selection_simulations")
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)


results <- readRDS(file.path(OUTPUT_DIR,"ongoing_cin_simulation_results.rds"))

# -----------------------------------------------------------------------------
# Horizon-anchored threshold
# -----------------------------------------------------------------------------
k_generations <- 2

results <- results %>%
  mutate(
    t_horizon    = log(n_sample) / log(N_final),
    t_thresh     = (log(n_sample) - k_generations * log(2)) / log(N_final),
    ground_truth = as.integer(t_mut > t_thresh)
  )

thresh_summary <- results %>%
  group_by(sample_fraction) %>%
  summarize(
    n_sample = round(mean(n_sample)),
    horizon  = mean(t_horizon),
    t_thresh = mean(t_thresh),
    .groups  = "drop"
  )
cat("\nHorizon-anchored threshold per sample fraction (k =", k_generations, "):\n")
print(as.data.frame(thresh_summary), row.names = FALSE)

# -----------------------------------------------------------------------------
# Presentation helpers
# -----------------------------------------------------------------------------
s_order      <- sort(unique(results$s))
s_order_char <- as.character(s_order)
s_colors <- setNames(
  c("#2166AC", "#4393C3", "#92C5DE", "#999999",
    "#F4A582", "#D6604D", "#B2182B", "#67001F"),
  s_order_char
)
samp_label  <- function(sf) paste0(sf * 100, "% sampled")
samp_levels <- samp_label(sort(unique(results$sample_fraction)))
samp_colors <- setNames(
  viridis::viridis(length(samp_levels), option = "D", end = 0.85),
  samp_levels
)
prep <- function(df) {
  df %>%
    mutate(
      s_f    = factor(s, levels = s_order),
      samp_f = factor(samp_label(sample_fraction), levels = samp_levels)
    )
}

# -----------------------------------------------------------------------------
# Supp Note 1 --> Figure 1: threshold diagram
#   Visualises the horizon and the threshold (last k generations before
#   the horizon) as functions of sample fraction.
# -----------------------------------------------------------------------------
sf_grid <- 10^seq(log10(0.005), log10(1), length.out = 200)
N <- mean(results$N_final)
diag_df <- tibble(
  sample_fraction = sf_grid,
  horizon  = log(pmax(1, sf_grid * N)) / log(N),
  t_thresh = (log(pmax(1, sf_grid * N)) - k_generations * log(2)) / log(N)
) %>%
  pivot_longer(c(horizon, t_thresh),
               names_to = "boundary", values_to = "t")

p_thresh <- ggplot(diag_df, aes(y = sample_fraction, x = t,
                                color = boundary, linewidth = boundary)) +
  geom_ribbon(data = diag_df %>%
    pivot_wider(names_from = boundary, values_from = t),
              aes(y = sample_fraction, xmin = t_thresh, xmax = horizon),
              fill = "#999933", alpha = 0.15, inherit.aes = FALSE) +
  geom_line() +
  geom_point(data = thresh_summary,
             aes(y = sample_fraction, x = t_thresh),
             color = "black", size = 2.5, inherit.aes = FALSE) +
  geom_point(data = thresh_summary,
             aes(y = sample_fraction, x = horizon),
             color = "#2166AC", size = 2.5, inherit.aes = FALSE) +
  scale_color_manual(values = c(horizon = "#2166AC", t_thresh = "black"),
                     labels = c(horizon  = "Horizon  log(n)/log(N)",
                                t_thresh = sprintf("t_thresh  (last %d generations before horizon)", k_generations)),
                     name   = "Boundary") +
  scale_linewidth_manual(values = c(horizon = 0.8, t_thresh = 1.1),
                         guide = "none") +
  scale_y_log10(breaks = c(0.01, 0.05, 0.1, 0.25, 1.0),
                labels = c("1%", "5%", "10%", "25%", "100%")) +
  xlim(0, 1) +
  labs(y = "Sample fraction (log)",
       x = "t_mut") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "right")

ggsave(file.path(PLOTS_DIR,"threshold_diagram.pdf"), p_thresh, width = 7, height = 4.5)

# -----------------------------------------------------------------------------
# Supp Note 1 --> Figure 2: raw frequency scatter
# -----------------------------------------------------------------------------
df_scatter <- results %>% filter(is_detected) %>% prep()

summary_scatter <- df_scatter %>%
  group_by(t_mut, s, s_f, sample_fraction, samp_f) %>%
  summarize(median = median(freq_in_sample), .groups = "drop")

singleton_lines <- df_scatter %>%
  group_by(samp_f, n_sample) %>%
  summarize(y = 1 / n_sample[1], .groups = "drop")

p_scatter <- ggplot(df_scatter, aes(x = t_mut, y = freq_in_sample, color = s_f)) +
  geom_jitter(width = 0.008, alpha = 0.10, size = 0.5, stroke = 0) +
  geom_line(data = summary_scatter,
            aes(x = t_mut, y = median, color = s_f, group = s_f),
            linewidth = 0.75) +
  geom_hline(data = singleton_lines, aes(yintercept = y),
             linetype = "dashed", color = "grey30", linewidth = 0.35) +
  scale_y_log10(labels = label_number(drop0trailing = TRUE)) +
  scale_color_manual(values = s_colors, name = "s") +
  facet_wrap(~ samp_f, nrow = 1) +
  labs(x = "Relative time of mutation origin (t_mut)",
       y = "Mutant frequency in sample (log scale)",
       caption = "Points: individual simulations (jittered). Lines: per-s median. Dashed: singleton frequency 1/n_sample.") +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "right") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)))

ggsave(file.path(PLOTS_DIR,"scatter_freq_vs_trel.pdf"), p_scatter, width = 14, height = 3.5)

# -----------------------------------------------------------------------------
# Supp Note 1 --> Figure 4: per-(s, sample_fraction) classification metrics
# -----------------------------------------------------------------------------
metrics <- results %>%
  filter(is_detected) %>%
  group_by(s, sample_fraction) %>%
  summarize(
    n  = n(),
    TP = sum( is_singleton & ground_truth == 1),
    FP = sum( is_singleton & ground_truth == 0),
    TN = sum(!is_singleton & ground_truth == 0),
    FN = sum(!is_singleton & ground_truth == 1),
    sensitivity   = ifelse((TP + FN) == 0, NA, TP / (TP + FN)),
    specificity   = ifelse((TN + FP) == 0, NA, TN / (TN + FP)),
    PPV           = ifelse((TP + FP) == 0, NA, TP / (TP + FP)),
    F1            = ifelse(is.na(PPV) | is.na(sensitivity) | (PPV + sensitivity) == 0,
                           NA, 2 * PPV * sensitivity / (PPV + sensitivity)),
    bal_accuracy  = (sensitivity + specificity) / 2,
    base_rate     = (TP + FN) / n,
    .groups       = "drop"
  ) %>%
  mutate(samp_f = factor(samp_label(sample_fraction), levels = samp_levels))

write.csv(metrics, "per_s_metrics.csv", row.names = FALSE)

metrics_long <- metrics %>%
  select(s, samp_f, sensitivity, specificity, PPV, F1, bal_accuracy) %>%
  pivot_longer(c(PPV, sensitivity, specificity, F1, bal_accuracy),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric,
                         levels = c("PPV", "sensitivity", "specificity", "F1", "bal_accuracy"),
                         labels = c("PPV  = P(ongoing | singleton)",
                                    "Sensitivity  = P(singleton | ongoing)",
                                    "Specificity  = P(not singleton | historical)",
                                    "F1  = harmonic mean of PPV and sensitivity",
                                    "Balanced accuracy  = (sens + spec) / 2")))

p_metrics <- ggplot(metrics_long,
                    aes(x = s, y = value, color = samp_f, group = samp_f)) +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey70") +
  geom_line(linewidth = 0.85) +
  geom_point(size = 2.2) +
  facet_wrap(~ metric, nrow = 1) +
  scale_color_manual(values = samp_colors, name = "Sample fraction") +
  scale_x_continuous(breaks = s_order) +
  ylim(0, 1) +
  labs(x = "Selection coefficient s", y = NULL,
       caption = "Horizon-anchored threshold: t_thresh(n,k) = (log(n) - k*log(2))/log(N). Dotted line: random-classifier baseline at 0.5.") +
  theme_minimal(base_size = 11) +
  theme(strip.text      = element_text(face = "bold", size = 9),
        legend.position = "right",
        axis.text.x     = element_text(size = 8))

ggsave(file.path(PLOTS_DIR,"classification_metrics.pdf"), p_metrics, width = 16, height = 3.8)

# -----------------------------------------------------------------------------
# Supp Note 1 --> Figure 6: mutation extinction rate
# -----------------------------------------------------------------------------
ext_df <- results %>%
  filter(sample_fraction == 1.0) %>%
  group_by(s, t_mut) %>%
  summarize(extinction_rate = mean(mutation_extinct),
            n = n(), .groups = "drop") %>%
  mutate(s_f = factor(s, levels = s_order))

p_ext <- ggplot(ext_df, aes(x = t_mut, y = extinction_rate,
                            color = s_f, group = s_f)) +
  geom_line(linewidth = 0.7, alpha = 0.95) +
  geom_point(size = 1.6) +
  scale_color_manual(values = s_colors, name = "s") +
  ylim(0, 1) +
  labs(x = "Relative time of mutation origin (t_mut)",
       y = "Mutation extinction rate",
       caption = "Fraction of simulations where the mutant lineage went extinct before sampling.") +
  theme_minimal(base_size = 11)

ggsave(file.path(PLOTS_DIR,"extinction_rate.pdf"), p_ext, width = 7, height = 4)

# -----------------------------------------------------------------------------
# DFE marginalisation
# -----------------------------------------------------------------------------
per_s <- results %>%
  filter(is_detected) %>%
  group_by(s, sample_fraction) %>%
  summarize(TP = sum( is_singleton & ground_truth == 1),
            FP = sum( is_singleton & ground_truth == 0),
            TN = sum(!is_singleton & ground_truth == 0),
            FN = sum(!is_singleton & ground_truth == 1),
            .groups = "drop")

dfe_scenarios <- list(
  "Neutral only"                = c("0" = 1),
  "Symmetric (mild)"            = c("-0.2" = 0.10, "-0.1" = 0.20, "0" = 0.40,
                                    "0.1"  = 0.20, "0.2"  = 0.10),
  "Deleterious-skewed"          = c("-0.5" = 0.10, "-0.2" = 0.25, "-0.1" = 0.25,
                                    "0"    = 0.30, "0.1"  = 0.10),
  "Beneficial-skewed"           = c("-0.1" = 0.10, "0" = 0.30, "0.1" = 0.25,
                                    "0.2"  = 0.20, "0.5" = 0.10, "1"  = 0.05),
  "Uniform (broad uncertainty)" = c("-0.5" = 1, "-0.2" = 1, "-0.1" = 1, "0" = 1,
                                    "0.1"  = 1, "0.2"  = 1, "0.5"  = 1, "1" = 1)
)
normalise_dfe <- function(dfe) data.frame(
  s = as.numeric(names(dfe)),
  w = as.numeric(dfe) / sum(as.numeric(dfe)))
dfe_scenarios <- lapply(dfe_scenarios, normalise_dfe)

marginalise <- function(per_s, dfe) {
  per_s %>%
    inner_join(dfe, by = "s") %>%
    group_by(sample_fraction) %>%
    summarize(
      wTP = sum(w * TP), wFP = sum(w * FP),
      wTN = sum(w * TN), wFN = sum(w * FN),
      sensitivity  = wTP / (wTP + wFN),
      specificity  = wTN / (wTN + wFP),
      PPV          = ifelse(wTP + wFP == 0, NA, wTP / (wTP + wFP)),
      F1           = 2 * PPV * sensitivity / (PPV + sensitivity),
      bal_accuracy = (sensitivity + specificity) / 2,
      .groups      = "drop")
}

marg_df <- bind_rows(lapply(names(dfe_scenarios), function(name) {
  marginalise(per_s, dfe_scenarios[[name]]) %>% mutate(scenario = name)
})) %>% mutate(scenario = factor(scenario, levels = names(dfe_scenarios)))

write.csv(marg_df, file.path(OUTPUT_DIR,"marginal_summary.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# Supp Note 1 --> Figure 7: the DFE shapes
# -----------------------------------------------------------------------------
dfe_long <- bind_rows(lapply(names(dfe_scenarios), function(name) {
  dfe_scenarios[[name]] %>% mutate(scenario = name)
})) %>%
  mutate(scenario = factor(scenario, levels = names(dfe_scenarios)))

scenario_colors <- c(
  "Neutral only"                = "#404040",
  "Symmetric (mild)"            = "#4393C3",
  "Deleterious-skewed"          = "#2166AC",
  "Beneficial-skewed"           = "#B2182B",
  "Uniform (broad uncertainty)" = "#999933")

p_dfe <- ggplot(dfe_long, aes(x = s, y = w, fill = scenario)) +
  geom_col(width = 0.08) +
  geom_hline(yintercept = 0, color = "grey30", linewidth = 0.3) +
  facet_wrap(~ scenario, nrow = 1) +
  scale_fill_manual(values = scenario_colors, guide = "none") +
  scale_x_continuous(breaks = c(-0.5, -0.2, 0, 0.2, 0.5, 1)) +
  labs(x = "Selection coefficient s", y = "DFE weight p(s)",
       caption = "Discrete DFE scenarios evaluated on the simulation grid.") +
  theme_minimal(base_size = 10) +
  theme(strip.text  = element_text(face = "bold", size = 9),
        axis.text.x = element_text(size = 7))

ggsave(file.path(PLOTS_DIR,"dfe_scenarios.pdf"), p_dfe, width = 13, height = 2.5)

# -----------------------------------------------------------------------------
# Supp Note 1 --> Figure 8: marginal metrics under each DFE
# -----------------------------------------------------------------------------
marg_long <- marg_df %>%
  select(scenario, sample_fraction, PPV, sensitivity, specificity, F1, bal_accuracy) %>%
  pivot_longer(c(PPV, sensitivity, specificity, F1, bal_accuracy),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric,
                         levels = c("PPV", "sensitivity", "specificity", "F1", "bal_accuracy"),
                         labels = c("PPV  = P(ongoing | singleton)",
                                    "Sensitivity  = P(singleton | ongoing)",
                                    "Specificity  = P(not singleton | historical)",
                                    "F1  = harmonic mean of PPV and sensitivity",
                                    "Balanced accuracy  = (sens + spec) / 2")))

p_marg <- ggplot(marg_long, aes(x = sample_fraction, y = value,
                                color = scenario, group = scenario)) +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey70") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.4) +
  facet_wrap(~ metric, nrow = 1) +
  scale_color_manual(values = scenario_colors, name = "DFE scenario") +
  scale_x_log10(breaks = c(0.01, 0.05, 0.1, 0.25, 1.0),
                labels = c("1%", "5%", "10%", "25%", "100%")) +
  ylim(0, 1) +
  labs(x = "Sample fraction (log)", y = NULL,
       caption = "DFE-marginalised metrics under the horizon-anchored threshold.") +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold", size = 9))

ggsave(file.path(PLOTS_DIR,"marginal_metrics.pdf"), p_marg, width = 16, height = 3.8)

# -----------------------------------------------------------------------------
# Concise summary
# -----------------------------------------------------------------------------
cat("\nMarginal classifier metrics under the horizon-anchored threshold:\n\n")
print(as.data.frame(marg_df %>%
                      mutate(across(c(PPV, sensitivity, specificity, F1, bal_accuracy),
                                    ~ round(.x, 2))) %>%
                      arrange(scenario, sample_fraction)),
      row.names = FALSE)

message("\nSaved:",
        "\n  threshold_diagram.pdf",
        "\n  scatter_freq_vs_trel.pdf",
        "\n  classification_metrics.pdf",
        "\n  extinction_rate.pdf",
        "\n  dfe_scenarios.pdf",
        "\n  marginal_metrics.pdf",
        "\n  per_s_metrics.csv",
        "\n  marginal_summary.csv")