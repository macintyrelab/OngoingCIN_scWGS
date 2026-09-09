# =============================================================================
# Ongoing vs historical CIN: identifiability under selection
# -----------------------------------------------------------------------------
# Single-mutation birth-death Gillespie simulation to quantify how selection
# distorts the singleton-as-ongoing-CIN interpretation across sample depths.
#
# Binary classifier version (no continuous score, no ROC):
#   prediction    = mutation is a singleton in the sample (yes / no)
#   ground truth  = mutation arose after t_thresh, where
#                   t_thresh = log(n_sample) / log(N_final)
#                   (log-corrected for exponential growth; under neutrality
#                    this is the arrival time below which a mutation's
#                    expected frequency exceeds the sampling detection limit)
#
# Outputs
#   scatter_freq_vs_trel.pdf         mutation frequency vs t_mut, coloured by s
#   singleton_prob_curves.pdf        P(singleton | t_mut) per (s, sample frac)
#   classification_metrics.pdf       sens / spec / PPV per (s, sample frac)
#   extinction_rate.pdf              mutation extinction rate vs t_mut, per s
#   ongoing_cin_simulation_results.rds
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(scales)
  library(viridis)
  library(this.path)
})

set.seed(42)

# PATHS
BASE=dirname(dirname(this.path()))
SIMS_DIR=dirname(this.path())
# Create directories for outputs
OUTPUT_DIR=file.path(SIMS_DIR,"outputs/selection_simulations")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)


# -----------------------------------------------------------------------------
# 1. Core simulator: one-mutation birth-death Gillespie
# -----------------------------------------------------------------------------
simulate_bd_single_mutation <- function(
  b        = 1.0,       # wild-type birth rate
  d        = 0.2,       # death rate (same for WT and mutant)
  s        = 0,         # selection coefficient for the mutation
  t_mut    = 0.5,       # fraction of expected total time at which mutation arises
  N_target = 10000,      # stop when total cell number reaches this
  max_iter = 5e6
) {
  T_expected  <- log(N_target) / (b - d)
  t_introduce <- t_mut * T_expected
  b_mut       <- b * (1 + s)

  n_wt <- 1L; n_mut <- 0L; t <- 0
  mutation_introduced <- FALSE

  for (i in seq_len(max_iter)) {
    n_tot <- n_wt + n_mut
    if (n_tot >= N_target) break
    if (n_tot == 0)        break

    r_bw <- b     * n_wt;  r_dw <- d * n_wt
    r_bm <- b_mut * n_mut; r_dm <- d * n_mut
    R    <- r_bw + r_dw + r_bm + r_dm
    if (R == 0) break

    dt     <- rexp(1, rate = R)
    t_next <- t + dt

    if (!mutation_introduced && t_next >= t_introduce && n_wt >= 1) {
      t     <- t_introduce
      n_wt  <- n_wt  - 1L
      n_mut <- n_mut + 1L
      mutation_introduced <- TRUE
      next
    }

    t <- t_next
    u <- runif(1) * R
    if      (u < r_bw)                 n_wt  <- n_wt  + 1L
    else if (u < r_bw + r_dw)          n_wt  <- n_wt  - 1L
    else if (u < r_bw + r_dw + r_bm)   n_mut <- n_mut + 1L
    else                               n_mut <- n_mut - 1L
  }

  list(
    n_wt                = n_wt,
    n_mut               = n_mut,
    n_total             = n_wt + n_mut,
    t_final             = t,
    mutation_introduced = mutation_introduced,
    mutation_extinct    = mutation_introduced && n_mut == 0,
    pop_extinct         = (n_wt + n_mut) == 0
  )
}

# -----------------------------------------------------------------------------
# 2. Sample the final population at several depths (hypergeometric)
# -----------------------------------------------------------------------------
sample_mutation_freq <- function(
  sim,
  sample_fractions = c(0.01, 0.05, 0.10, 0.25, 1.00)
) {
  if (sim$pop_extinct || !sim$mutation_introduced) return(NULL)
  N <- sim$n_total
  do.call(rbind, lapply(sample_fractions, function(f) {
    n_sample   <- max(1L, round(f * N))
    n_mut_samp <- rhyper(1, m = sim$n_mut, n = sim$n_wt, k = n_sample)
    data.frame(
      sample_fraction = f,
      n_sample        = n_sample,
      n_mut_sampled   = n_mut_samp,
      freq_in_sample  = n_mut_samp / n_sample,
      is_singleton    = n_mut_samp == 1L,
      is_detected     = n_mut_samp >= 1L
    )
  }))
}

# -----------------------------------------------------------------------------
# 3. Grid driver
# -----------------------------------------------------------------------------
run_grid <- function(
  t_mut_grid       = seq(0.05, 0.95, by = 0.05),
  s_grid           = c(-0.5, -0.2, -0.1, 0, 0.1, 0.2, 0.5, 1.0),
  sample_fractions = c(0.01, 0.05, 0.10, 0.25, 1.00),
  n_reps           = 100,
  b = 1.0, d = 0.2, N_target = 10000,
  verbose = TRUE
) {
  grid <- expand.grid(t_mut = t_mut_grid, s = s_grid, rep = seq_len(n_reps),
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  if (verbose) message("Running ", nrow(grid), " simulations...")

  out <- vector("list", nrow(grid))
  pb  <- txtProgressBar(min = 0, max = nrow(grid), style = 3)
  for (i in seq_len(nrow(grid))) {
    setTxtProgressBar(pb, i)
    sim <- simulate_bd_single_mutation(
      b = b, d = d, s = grid$s[i], t_mut = grid$t_mut[i], N_target = N_target
    )
    samp <- sample_mutation_freq(sim, sample_fractions)
    if (is.null(samp)) next
    samp$t_mut            <- grid$t_mut[i]
    samp$s                <- grid$s[i]
    samp$rep              <- grid$rep[i]
    samp$N_final          <- sim$n_total
    samp$pop_freq         <- sim$n_mut / sim$n_total
    samp$mutation_extinct <- sim$mutation_extinct
    out[[i]] <- samp
  }
  close(pb)
  bind_rows(out)
}

# -----------------------------------------------------------------------------
# 4. Log-corrected ground truth
# -----------------------------------------------------------------------------
annotate_ground_truth <- function(df) {
  df %>%
    mutate(
      t_thresh_log = log(n_sample) / log(N_final),
      ground_truth = as.integer(t_mut > t_thresh_log)
    )
}

# -----------------------------------------------------------------------------
# 9. Run it
# -----------------------------------------------------------------------------
if (sys.nframe() == 0) {
  results <- run_grid(
    t_mut_grid       = seq(0.05, 1, by = 0.05),
    s_grid           = c(-0.5, -0.2, -0.1, 0, 0.1, 0.2, 0.5, 1.0),
    sample_fractions = c(0.01, 0.05, 0.10, 0.25, 1.00),
    n_reps           = 100,
    b = 1.0, d = 0.2, N_target = 10000
  )
  results <- annotate_ground_truth(results)
  saveRDS(results, file.path(OUTPUT_DIR,"ongoing_cin_simulation_results.rds"))

  message("Done.")
}