# =============================================================================
# RVF Climate Modelling — Bayesian ZINB and Hurdle Negative Binomial
# Fits both model families in brms and compares via LOO-CV
# =============================================================================
# Required packages:
#   install.packages(c("tidyverse", "brms", "loo", "bayesplot",
#                      "tidybayes", "patchwork", "viridis"))
# =============================================================================

library(tidyverse)
library(brms)
library(loo)
library(bayesplot)
library(tidybayes)
library(patchwork)
library(viridis)

# Resolve namespace conflicts
select <- dplyr::select
filter <- dplyr::filter

# =============================================================================
# 0. LOAD DATA
# =============================================================================
# Replace with your actual data load:
# dat <- read_csv("your_data.csv")
#
# This script assumes the same column structure as the exploration script:
#   cases, zone, year, precip, temp, pdsi, pet, soil_moisture,
#   livestock, year_period

# --- SYNTHETIC DATA (remove when using real data) ---
set.seed(42)
dat <- expand.grid(
  zone = paste0("AEZ_", 1:7),
  year = 1980:2021
) %>%
  mutate(
    precip        = rnorm(n(), 600, 150),
    temp          = rnorm(n(), 22, 3),
    pdsi          = rnorm(n(), 0, 2),
    pet           = rnorm(n(), 1200, 200) + 0.3 * temp * 10,
    soil_moisture = rnorm(n(), 30, 8) + 0.4 * (precip - 600) / 150,
    livestock     = round(runif(n(), 5000, 50000)),
    cases         = rnbinom(n(), mu = exp(0.003 * precip - 0.1 * pdsi), size = 0.5),
    cases         = ifelse(runif(n()) < 0.6, 0L, cases),
    year_period   = ifelse(year < 2000, "pre2000", "post2000"),
    zone          = as.factor(zone)
  )
# --- END SYNTHETIC DATA ---

# =============================================================================
# 1. SCALE COVARIATES
# =============================================================================
# All continuous covariates are scaled to mean = 0, SD = 1.
# This serves three purposes:
#   (a) Aids MCMC sampling and chain convergence
#   (b) Makes regression coefficients directly comparable in magnitude
#   (c) Allows meaningful weakly-informative priors on a consistent scale
#
# Scaling parameters are stored so results can be back-transformed for
# reporting if needed.

covariates <- c("precip", "temp", "pdsi", "pet", "soil_moisture")

scaling_params <- dat %>%
  summarise(across(all_of(covariates),
                   list(mean = mean, sd = sd),
                   .names = "{.col}__{.fn}"))

cat("===== SCALING PARAMETERS =====\n")
print(scaling_params)

dat_scaled <- dat %>%
  mutate(across(all_of(covariates),
                ~ (. - mean(.)) / sd(.),
                .names = "{.col}_z")) %>%
  mutate(
    log_livestock = log(livestock),   # pre-compute offset
    zone          = as.factor(zone)
  )

# Quick check — scaled covariates should have mean ~0 and SD ~1
cat("\nScaled covariate means (should be ~0):\n")
dat_scaled %>%
  summarise(across(ends_with("_z"), mean)) %>%
  print()

# =============================================================================
# 2. PRIOR SPECIFICATION
# =============================================================================
# Weakly informative priors on the standardised scale:
#   Normal(0, 1) for fixed effects — allows IRR roughly between 0.05 and 20,
#   which is broad enough to be non-constraining but prevents implausible values.
#   Exponential(1) for random effect SDs — places most prior mass below 2,
#   appropriate for log-scale variance components.
#   Beta(1, 1) / logistic-Normal(0, 1.5) for zero-inflation / hurdle
#   probability — weakly informative on the probability scale.
#
# These are the same priors for both models to ensure fair LOO comparison.

priors_shared <- c(
  prior(normal(0, 1),   class = b),           # fixed effects — count component
  prior(exponential(1), class = sd),           # random effect SDs — count component
  prior(normal(0, 1),   class = b,    dpar = zi),  # fixed effects — ZI component
  prior(exponential(1), class = sd,   dpar = zi),  # random effect SDs — ZI component
  prior(normal(0, 1),   class = b,    dpar = hu),  # fixed effects — hurdle component
  prior(exponential(1), class = sd,   dpar = hu)   # random effect SDs — hurdle component
)
# Note: brms ignores priors for dpars not present in a given model family,
# so it is safe to specify all three sets here and use the same object for both.

# =============================================================================
# 3. FIT ZINB MODEL
# =============================================================================
# Zero-inflated negative binomial.
# Count component: climate covariates + AEZ random intercept + AR(1)
# ZI component:    year_period + AEZ random intercept
#
# The ZI component models the probability of being a structural zero
# (non-detection / surveillance failure). year_period captures secular
# improvements in reporting over time.

cat("\n===== FITTING ZINB MODEL =====\n")

formula_zinb <- bf(
  cases ~ precip_z + temp_z + pdsi_z + pet_z + soil_moisture_z +
          offset(log_livestock) +
          (1 | zone) +
          ar(time = year, gr = zone, p = 1),
  zi    ~ year_period + (1 | zone)
)

fit_zinb <- brm(
  formula = formula_zinb,
  data    = dat_scaled,
  family  = zero_inflated_negbinomial(),
  prior   = priors_shared,
  chains  = 4,
  cores   = 4,           # set to number of available CPU cores
  iter    = 4000,
  warmup  = 1000,
  seed    = 42,
  control = list(
    adapt_delta   = 0.95,  # increase if you get divergent transitions
    max_treedepth = 12
  ),
  file    = "fit_zinb"    # saves fitted object — comment out to refit each time
)

cat("\n--- ZINB Model Summary ---\n")
print(summary(fit_zinb))

# =============================================================================
# 4. FIT HURDLE NEGATIVE BINOMIAL MODEL
# =============================================================================
# Hurdle negative binomial.
# Count component: climate covariates + AEZ random intercept + AR(1)
#                  (truncated at zero — only observations > 0 contribute)
# Hurdle component: year_period + AEZ random intercept
#                  (models probability of a zero — positive coefficient
#                   on year_period means MORE zeros in that period)
#
# Note on direction: in brms, hu models P(zero), so if reporting improved
# post-2000 you would expect FEWER zeros, i.e. a negative coefficient on
# the post2000 indicator. Check this is consistent with your data.

cat("\n===== FITTING HURDLE NEGATIVE BINOMIAL MODEL =====\n")

formula_hurdle <- bf(
  cases ~ precip_z + temp_z + pdsi_z + pet_z + soil_moisture_z +
          offset(log_livestock) +
          (1 | zone) +
          ar(time = year, gr = zone, p = 1),
  hu    ~ year_period + (1 | zone)
)

fit_hurdle <- brm(
  formula = formula_hurdle,
  data    = dat_scaled,
  family  = hurdle_negbinomial(),
  prior   = priors_shared,
  chains  = 4,
  cores   = 4,
  iter    = 4000,
  warmup  = 1000,
  seed    = 42,
  control = list(
    adapt_delta   = 0.95,
    max_treedepth = 12
  ),
  file    = "fit_hurdle"
)

cat("\n--- Hurdle NB Model Summary ---\n")
print(summary(fit_hurdle))

# =============================================================================
# 5. MCMC DIAGNOSTICS
# =============================================================================
# Check before trusting any results. Key diagnostics:
#   Rhat: should be < 1.01 for all parameters
#   Bulk and tail ESS: should be > 400 (ideally > 1000)
#   Divergent transitions: should be 0 (increase adapt_delta if not)
#   Trace plots: chains should mix well with no trends or sticking

check_diagnostics <- function(fit, model_name) {
  cat(sprintf("\n===== MCMC DIAGNOSTICS — %s =====\n", model_name))

  rhat_vals <- rhat(fit)
  ess_bulk  <- ess_bulk(fit)
  ess_tail  <- ess_tail(fit)

  cat(sprintf("Max Rhat          : %.4f (should be < 1.01)\n", max(rhat_vals, na.rm = TRUE)))
  cat(sprintf("Min bulk ESS      : %.0f  (should be > 400)\n",  min(ess_bulk, na.rm = TRUE)))
  cat(sprintf("Min tail ESS      : %.0f  (should be > 400)\n",  min(ess_tail, na.rm = TRUE)))
  cat(sprintf("Divergent transitions: %d  (should be 0)\n",
              sum(nuts_params(fit)$Value[nuts_params(fit)$Parameter == "divergent__"])))

  # Flag any problematic parameters
  prob_rhat <- names(rhat_vals)[rhat_vals > 1.01 & !is.na(rhat_vals)]
  if (length(prob_rhat) > 0) {
    cat("Parameters with Rhat > 1.01:\n")
    print(prob_rhat)
  } else {
    cat("All Rhat values < 1.01 -- good convergence.\n")
  }
}

check_diagnostics(fit_zinb,   "ZINB")
check_diagnostics(fit_hurdle, "Hurdle NB")

# Trace plots — visually inspect chain mixing for key parameters
# (saves a plot for each model; open and inspect manually)
p_trace_zinb <- mcmc_trace(
  fit_zinb,
  pars = c("b_precip_z", "b_temp_z", "b_pdsi_z", "b_pet_z",
           "b_soil_moisture_z", "shape"),
  facet_args = list(ncol = 2)
) +
  labs(title = "Trace plots — ZINB (count component)") +
  theme_bw(base_size = 10)

ggsave("09_trace_zinb.png", p_trace_zinb, width = 10, height = 8, dpi = 300)

p_trace_hurdle <- mcmc_trace(
  fit_hurdle,
  pars = c("b_precip_z", "b_temp_z", "b_pdsi_z", "b_pet_z",
           "b_soil_moisture_z", "shape"),
  facet_args = list(ncol = 2)
) +
  labs(title = "Trace plots — Hurdle NB (count component)") +
  theme_bw(base_size = 10)

ggsave("10_trace_hurdle.png", p_trace_hurdle, width = 10, height = 8, dpi = 300)

# =============================================================================
# 6. POSTERIOR PREDICTIVE CHECKS
# =============================================================================
# Compares observed data against data simulated from the posterior predictive
# distribution. Key checks for your data:
#   (a) Overall distribution — does the model reproduce the general shape?
#   (b) Zero proportion — does the model reproduce your ~75% zeros?
#   (c) Outbreak years — does the model capture the heavy right tail?

ppc_check <- function(fit, model_name, file_prefix) {

  y     <- dat_scaled$cases
  y_rep <- posterior_predict(fit, ndraws = 200)

  # Overall distribution
  p1 <- ppc_dens_overlay(y, y_rep[1:50, ]) +
    coord_cartesian(xlim = c(0, quantile(y, 0.99))) +
    labs(title = sprintf("%s — posterior predictive density", model_name)) +
    theme_bw(base_size = 11)

  # Zero proportion
  prop_zero <- function(x) mean(x == 0)
  p2 <- ppc_stat(y, y_rep, stat = "prop_zero") +
    labs(title = sprintf("%s — zero proportion", model_name),
         x = "Proportion of zeros") +
    theme_bw(base_size = 11)

  # Mean
  p3 <- ppc_stat(y, y_rep, stat = "mean") +
    labs(title = sprintf("%s — mean", model_name)) +
    theme_bw(base_size = 11)

  # Standard deviation (captures spread / heavy tail)
  p4 <- ppc_stat(y, y_rep, stat = "sd") +
    labs(title = sprintf("%s — SD", model_name)) +
    theme_bw(base_size = 11)

  combined <- (p1 + p2) / (p3 + p4)
  ggsave(sprintf("%s_ppc.png", file_prefix), combined,
         width = 12, height = 8, dpi = 300)

  cat(sprintf("\nPPC saved: %s_ppc.png\n", file_prefix))
}

cat("\n===== POSTERIOR PREDICTIVE CHECKS =====\n")
ppc_check(fit_zinb,   "ZINB",      "11_zinb")
ppc_check(fit_hurdle, "Hurdle NB", "12_hurdle")

# =============================================================================
# 7. LOO-CV MODEL COMPARISON
# =============================================================================
# Leave-one-out cross-validation via PSIS-LOO.
# The model with the higher ELPD (less negative) has better predictive fit.
# elpd_diff / se_diff > 2 is a commonly used threshold for a meaningful
# difference, though this is a guideline rather than a formal test.
#
# Note: if you get warnings about high Pareto k values (k > 0.7), this
# indicates some observations are highly influential. brms will suggest
# using reloo() for exact LOO on those observations — follow that advice
# if more than a handful of k values are problematic.

cat("\n===== LOO-CV MODEL COMPARISON =====\n")

loo_zinb   <- loo(fit_zinb,   moment_match = TRUE, reloo = FALSE)
loo_hurdle <- loo(fit_hurdle, moment_match = TRUE, reloo = FALSE)

cat("\n--- ZINB LOO ---\n")
print(loo_zinb)

cat("\n--- Hurdle NB LOO ---\n")
print(loo_hurdle)

cat("\n--- Model Comparison ---\n")
loo_comp <- loo_compare(loo_zinb, loo_hurdle)
print(loo_comp)

# Interpret the comparison
winner <- rownames(loo_comp)[1]
elpd_diff <- loo_comp[2, "elpd_diff"]
se_diff   <- loo_comp[2, "se"]
cat(sprintf(
  "\nPreferred model: %s\nELPD difference: %.2f (SE = %.2f)\n%s\n",
  winner, abs(elpd_diff), se_diff,
  ifelse(abs(elpd_diff) > 2 * se_diff,
         "--> Meaningful difference in predictive performance.",
         "--> Difference is within 2 SE — models perform similarly.")
))

# Pareto k diagnostic plot
png("13_loo_pareto_k.png", width = 1200, height = 500, res = 150)
par(mfrow = c(1, 2))
plot(loo_zinb,   main = "Pareto k — ZINB",      label_points = TRUE)
plot(loo_hurdle, main = "Pareto k — Hurdle NB", label_points = TRUE)
par(mfrow = c(1, 1))
dev.off()

# =============================================================================
# 8. COEFFICIENT PLOT — PREFERRED MODEL
# =============================================================================
# Plots posterior distributions of fixed effects from the count component.
# Exponentiated to incidence rate ratios (IRR) for interpretability.
# Reminder: coefficients are on the scaled covariate scale — IRR represents
# the change in case rate per 1 SD increase in each covariate.

plot_coefficients <- function(fit, model_name, file_name) {

  # Count component fixed effects only (excludes intercept and ZI/hu params)
  coef_vars <- paste0("b_", c("precip_z", "temp_z", "pdsi_z",
                               "pet_z", "soil_moisture_z"))

  draws <- as_draws_df(fit) %>%
    dplyr::select(any_of(coef_vars)) %>%
    pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
    mutate(
      IRR       = exp(value),
      parameter = str_remove(parameter, "^b_") %>% str_remove("_z$")
    )

  p <- ggplot(draws, aes(x = IRR, y = reorder(parameter, IRR))) +
    stat_halfeye(
      aes(fill = after_stat(x > 1)),
      .width       = c(0.80, 0.95),
      point_interval = median_qi,
      normalize    = "panels"
    ) +
    geom_vline(xintercept = 1, linetype = "dashed",
               colour = "grey30", linewidth = 0.7) +
    scale_fill_manual(values = c("TRUE"  = viridis(1, begin = 0.6),
                                  "FALSE" = viridis(1, begin = 0.2)),
                       guide = "none") +
    scale_x_log10() +
    labs(
      title    = sprintf("%s — posterior IRRs (count component)", model_name),
      subtitle = "IRR per 1 SD increase in covariate | inner bar = 80% CI, outer = 95% CI",
      x        = "Incidence Rate Ratio (log scale)",
      y        = NULL
    ) +
    theme_bw(base_size = 12)

  ggsave(file_name, p, width = 8, height = 5, dpi = 300)
  cat(sprintf("Coefficient plot saved: %s\n", file_name))
}

cat("\n===== COEFFICIENT PLOTS =====\n")
plot_coefficients(fit_zinb,   "ZINB",      "14_coef_zinb.png")
plot_coefficients(fit_hurdle, "Hurdle NB", "15_coef_hurdle.png")

# =============================================================================
# 9. SAVE SESSION AND MODEL SUMMARIES
# =============================================================================

sink("16_model_summaries.txt")
cat("===== ZINB MODEL SUMMARY =====\n")
print(summary(fit_zinb))
cat("\n\n===== HURDLE NB MODEL SUMMARY =====\n")
print(summary(fit_hurdle))
cat("\n\n===== LOO COMPARISON =====\n")
print(loo_comp)
sink()

cat("\nAll outputs saved.\n")
cat("Model objects cached as fit_zinb.rds and fit_hurdle.rds\n")
cat("(delete these files if you want to refit from scratch)\n")

# =============================================================================
# NOTES
# =============================================================================
# 1. COMPUTATION TIME:
#    With 4 chains x 4000 iterations on real data, expect 30-90 minutes
#    depending on your machine and data size. The file = argument caches
#    fitted models so you do not need to refit every session.
#
# 2. DIVERGENT TRANSITIONS:
#    If you get divergent transitions after fitting, increase adapt_delta
#    toward 0.99 and refit. Delete the cached .rds file first.
#
# 3. HIGH PARETO K VALUES:
#    If loo() warns about many k > 0.7, run:
#    loo_zinb <- loo(fit_zinb, moment_match = TRUE, reloo = TRUE)
#    This refits the model for influential observations exactly — slower
#    but produces reliable LOO estimates.
#
# 4. BACK-TRANSFORMING COEFFICIENTS:
#    To report effects on the original (unscaled) covariate scale, divide
#    the posterior draws by the SD of that covariate stored in scaling_params.
#    The IRR interpretation does not change, only the "per unit" reference.
#
# 5. ZI vs HU COMPONENT DIRECTION:
#    In ZINB, zi models P(structural zero). Higher zi = more structural zeros.
#    In hurdle, hu models P(zero). Higher hu = more zeros.
#    A negative coefficient on year_period (post2000) in either component
#    would indicate fewer zeros post-2000, consistent with improved reporting.
#    Check this direction against your expectations.
