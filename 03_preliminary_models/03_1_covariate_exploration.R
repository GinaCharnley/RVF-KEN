### FITTING THE PRELIMINARY MODELS 

# Load pacakges 
library(tidyverse)
library(readr)
library(corrplot)
library(car)          
library(lme4)         
library(performance)  
library(MASS)       
library(patchwork)
library(viridis)
library(pscl)         
library(AER)  

# Load data
dat <- read_csv("01_data_processing/fitting_data.csv")

# Set columns to correct names and format 
dat <- dat %>% group_by(AEZ_Name) %>%
  mutate(livestock = sum(cattle_mean, goat_mean, sheep_mean, na.rm = T)) %>%
  ungroup %>%
  rename("cases" = Number_confirmed_events,
         "zone" = AEZ_Name,
         "year" = Year,
         "precip" = pre,
         "temp" = tmp,
         "soil_moisture" = sm) %>%
  mutate(zone = as.factor(zone))

covariates <- c("precip", "temp", "pdsi", "pet", "soil_moisture")

## DESCRIPTIVE SUMMARIES

dat %>% dplyr::select(cases, all_of(covariates)) %>% summary() %>% print()

zero_prop  <- mean(dat$cases == 0)
disp_ratio <- var(dat$cases) / mean(dat$cases)

cat(sprintf("\nZero proportion in cases : %.3f\n", zero_prop))
cat(sprintf("Variance / Mean ratio    : %.2f  (>1 suggests overdispersion)\n",
            disp_ratio))

# Distribution of cases
p_dist <- ggplot(dat, aes(x = cases)) +
  geom_histogram(bins = 40, fill = viridis(1, begin = 0.4),
                 colour = "white", alpha = 0.85) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "Distribution of annual RVF cases",
       subtitle = sprintf("Zero proportion = %.2f | Variance/Mean = %.1f",
                          zero_prop, disp_ratio),
       x = "Cases", y = "Count") +
  theme_bw(base_size = 12)

p_dist_log <- ggplot(dat %>% filter(cases > 0), aes(x = log1p(cases))) +
  geom_histogram(bins = 30, fill = viridis(1, begin = 0.65),
                 colour = "white", alpha = 0.85) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "log(cases + 1) — non-zero observations only",
       x = "log(cases + 1)", y = "Count") +
  theme_bw(base_size = 12)

ggsave("03_preliminary_models/case_distribution.png", p_dist + p_dist_log,
       width = 10, height = 4, dpi = 300)


## OVERDISPERSION TEST

# Formally tests whether a Poisson model is adequate or whether the negative
# binomial is needed. A significant result justifies using NB family in brms.
# H0: equidispersion (Poisson is adequate)

fit_pois <- glm(cases ~ precip + temp + pdsi + pet + soil_moisture +
                  offset(log(livestock)),
                data = dat, family = poisson(link = "log"))
disp_test <- dispersiontest(fit_pois, trafo = 1)
print(disp_test)
cat(ifelse(disp_test$p.value < 0.05,
           "--> Significant overdispersion: negative binomial family justified.\n",
           "--> No significant overdispersion detected.\n"))

## ZERO-INFLATION ASSESSMENT

# Compares observed zero count against zeros expected under Poisson and NB.
# A large excess over NB expectation justifies the hurdle or ZI component.

mu_pois              <- fitted(fit_pois)
expected_zeros_pois  <- sum(dpois(0, lambda = mu_pois))

fit_nb               <- glm.nb(cases ~ precip + temp + pdsi + pet +
                                  soil_moisture + offset(log(livestock)),
                                data = dat)
mu_nb                <- fitted(fit_nb)
theta_nb             <- fit_nb$theta
expected_zeros_nb    <- sum(dnbinom(0, mu = mu_nb, size = theta_nb))

observed_zeros <- sum(dat$cases == 0)
n_obs          <- nrow(dat)

zero_compare <- tibble(
  model   = c("Observed", "Expected (Poisson)", "Expected (Neg. Binomial)"),
  n_zeros = round(c(observed_zeros, expected_zeros_pois, expected_zeros_nb)),
  prop    = round(c(observed_zeros, expected_zeros_pois,
                    expected_zeros_nb) / n_obs, 3)
)
print(zero_compare)

excess_over_nb <- observed_zeros - round(expected_zeros_nb)
cat(sprintf(
  "\nExcess zeros over NB expectation: %d (%.1f%% of observations)\n",
  excess_over_nb, 100 * excess_over_nb / n_obs
))
cat(ifelse(excess_over_nb > 0,
           "--> Zero inflation beyond NB: hurdle or ZI component justified.\n",
           "--> NB accounts for zeros adequately; ZI component may not be needed.\n"))

p_zeros <- ggplot(zero_compare, aes(x = model, y = prop, fill = model)) +
  geom_col(width = 0.55) +
  scale_fill_viridis_d(option = "C", begin = 0.2, end = 0.85) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0),
                     limits = c(0, max(zero_compare$prop) * 1.15)) +
  labs(title = "Observed vs expected zero proportion",
       subtitle = "Excess over NB expectation motivates hurdle / ZI component",
       x = NULL, y = "Proportion of zeros") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

ggsave("03_preliminary_models/zero_inflation_check.png", p_zeros, width = 6, height = 4, dpi = 300)

## TEMPORAL AUTOCORRELATION (ACF) WITHIN ZONES

# Justifies the AR(1) term in brms. ACF of NB model residuals within each
# zone — significant lag-1 autocorrelation supports ar(time, gr = zone, p = 1).

dat <- dat %>%
  mutate(nb_resid = residuals(fit_nb, type = "pearson"))

acf_results <- dat %>%
  group_by(zone) %>%
  arrange(year, .by_group = TRUE) %>%
  summarise(
    acf_lag1 = acf(nb_resid, lag.max = 1, plot = FALSE)$acf[2],
    .groups = "drop"
  )

cat("Lag-1 autocorrelation of NB residuals by zone:\n")
print(acf_results)
cat(sprintf("Mean lag-1 autocorrelation across zones: %.3f\n",
            mean(acf_results$acf_lag1)))
cat(ifelse(mean(abs(acf_results$acf_lag1)) > 0.1,
           "--> Meaningful autocorrelation: AR(1) term justified in brms.\n",
           "--> Low autocorrelation: AR(1) term may not be necessary.\n"))

# ACF plots per zone
png("03_preliminary_models/acf_by_zone.png", width = 1400, height = 1000, res = 150)
par(mfrow = c(ceiling(n_distinct(dat$zone) / 3), 3),
    mar = c(3, 3, 2, 1))
for (z in levels(dat$zone)) {
  resids_z <- dat %>% filter(zone == z) %>% arrange(year) %>% pull(nb_resid)
  acf(resids_z, main = paste("ACF residuals —", z), lag.max = 10,
      col = viridis(1, begin = 0.5), lwd = 2)
}
par(mfrow = c(1, 1))
dev.off()

## INTRACLASS CORRELATION (ICC)

# Proportion of total variance attributable to between-zone differences.
# ICC > 0.05 is generally considered meaningful and justifies a random effect.

fit_icc  <- lmer(log1p(cases) ~ 1 + (1 | zone), data = dat, REML = TRUE)
icc_vals <- icc(fit_icc)
print(icc_vals)
cat(sprintf(
  "\nICC = %.3f: %.1f%% of variance is between zones.\n",
  icc_vals$ICC_adjusted,
  icc_vals$ICC_adjusted * 100
))
cat(ifelse(icc_vals$ICC_adjusted > 0.05,
           "--> Non-trivial between-zone variance: random intercept justified.\n",
           "--> Low ICC: random effects may contribute little.\n"))

## CORRELATION ANALYSIS

# Spearman preferred for climate variables (non-normal distributions).
# Pearson included for comparison.

covar_mat <- dat %>% dplyr::select(all_of(covariates)) %>% as.matrix()

cor_pearson  <- cor(covar_mat, method = "pearson",  use = "complete.obs")
cor_spearman <- cor(covar_mat, method = "spearman", use = "complete.obs")

cor_out_spearman <- cor(covar_mat, dat$cases,
                        method = "spearman", use = "complete.obs")
cor_out_pearson  <- cor(covar_mat, dat$cases,
                        method = "pearson",  use = "complete.obs")


print(round(cor_out_spearman, 3))

cor_pvals <- sapply(covariates, function(v) {
  c(
    pearson_p  = cor.test(dat[[v]], dat$cases,
                          method = "pearson")$p.value,
    spearman_p = cor.test(dat[[v]], dat$cases,
                          method = "spearman", exact = FALSE)$p.value
  )
})

print(round(t(cor_pvals), 4))

png("03_preliminary_models/correlation_spearman.png", width = 900, height = 800, res = 150)
corrplot(cor_spearman,
         method = "color", type = "upper", order = "hclust",
         addCoef.col = "black", number.cex = 0.8,
         tl.col = "black", tl.srt = 45,
         col = colorRampPalette(c("#2C7BB6", "white", "#D7191C"))(200),
         title = "Spearman correlation — covariates", mar = c(0, 0, 2, 0))
dev.off()

png("03_preliminary_models/correlation_pearson.png", width = 900, height = 800, res = 150)
corrplot(cor_pearson,
         method = "color", type = "upper", order = "hclust",
         addCoef.col = "black", number.cex = 0.8,
         tl.col = "black", tl.srt = 45,
         col = colorRampPalette(c("#2C7BB6", "white", "#D7191C"))(200),
         title = "Pearson correlation — covariates", mar = c(0, 0, 2, 0))
dev.off()

cat("\n===== HIGHLY CORRELATED COVARIATE PAIRS (Spearman |r| > 0.7) =====\n")
high_cor <- which(abs(cor_spearman) > 0.7 & abs(cor_spearman) < 1,
                  arr.ind = TRUE)
if (nrow(high_cor) > 0) {
  for (i in seq(1, nrow(high_cor), by = 2)) {
    r <- round(cor_spearman[high_cor[i, 1], high_cor[i, 2]], 3)
    cat(sprintf("  %s -- %s : r = %s\n",
                covariates[high_cor[i, 1]],
                covariates[high_cor[i, 2]], r))
  }
} else {
  cat("  No pairs above threshold.\n")
}

## VARIANCE INFLATION FACTORS (VIF)

vif_model <- glm(cases ~ precip + temp + pdsi + pet + soil_moisture,
                 data = dat, family = poisson(link = "log"))
vif_vals  <- vif(vif_model)
print(round(vif_vals, 3))

vif_df <- tibble(
  variable = names(vif_vals),
  VIF      = as.numeric(vif_vals),
  concern  = case_when(
    VIF >= 10 ~ "High (>10)",
    VIF >= 5  ~ "Moderate (5-10)",
    TRUE      ~ "Acceptable (<5)"
  )
)

p_vif <- ggplot(vif_df, aes(x = reorder(variable, VIF), y = VIF,
                              fill = concern)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 5,  linetype = "dashed",
             colour = "orange", linewidth = 0.8) +
  geom_hline(yintercept = 10, linetype = "dashed",
             colour = "red",   linewidth = 0.8) +
  annotate("text", x = 0.6, y = 5.4,  label = "VIF = 5",
           size = 3, colour = "orange", hjust = 0) +
  annotate("text", x = 0.6, y = 10.4, label = "VIF = 10",
           size = 3, colour = "red",   hjust = 0) +
  scale_fill_manual(values = c("Acceptable (<5)"  = "#2ca25f",
                                "Moderate (5-10)"  = "#fc8d59",
                                "High (>10)"       = "#d7301f")) +
  coord_flip() +
  labs(title = "Variance Inflation Factors",
       subtitle = "Dashed lines at VIF = 5 and VIF = 10",
       x = NULL, y = "VIF", fill = "Concern level") +
  theme_bw(base_size = 12)

ggsave("03_preliminary_models/vif.png", p_vif, width = 7, height = 4, dpi = 300)

## UNIVARIATE NEGATIVE BINOMIAL GLMs

# Fits a separate NB GLM per covariate with a livestock offset.
# Incidence rate ratios (IRR), 95% CIs, and delta AIC vs a null model provide
# a distribution-appropriate importance ranking — directly analogous to
# univariate screening before a multivariable model, and on the correct
# scale for count data.

null_nb  <- glm.nb(cases ~ offset(log(livestock)), data = dat)
null_aic <- AIC(null_nb)
cat("Null model AIC:", round(null_aic, 2), "\n")

dat <- dat %>%
  mutate(across(all_of(covariates), scale, .names = "{.col}_z"))

covariates_z <- paste0(covariates, "_z")

univariate_nb <- lapply(covariates_z, function(v) {
  form <- as.formula(paste("cases ~", v, "+ offset(log(livestock))"))
  fit  <- tryCatch(
    withCallingHandlers(
      glm.nb(form, data = dat, control = glm.control(maxit = 100)),
      warning = function(w) {
        if (grepl("did not converge", conditionMessage(w))) {
          message("Convergence warning for: ", v)
        }
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  coef_row <- summary(fit)$coefficients[v, ]
  ci <- tryCatch(
    exp(confint.default(fit)[v, ]),
    error = function(e) c(NA_real_, NA_real_)
  )
  
  tibble(
    variable  = v,
    IRR       = exp(coef_row["Estimate"]),
    CI_lower  = ci[1],
    CI_upper  = ci[2],
    z_value   = coef_row["z value"],
    p_value   = coef_row["Pr(>|z|)"],
    AIC       = AIC(fit),
    delta_AIC = AIC(fit) - null_aic
  )
}) %>% bind_rows() %>% arrange(delta_AIC)

print(univariate_nb %>% mutate(across(where(is.numeric), ~ round(., 4))))

# Forest plot of IRRs
p_irr <- ggplot(univariate_nb,
                aes(x = reorder(variable, IRR), y = IRR,
                    ymin = CI_lower, ymax = CI_upper,
                    colour = p_value < 0.05)) +
  geom_pointrange(size = 0.8, linewidth = 0.9) +
  geom_hline(yintercept = 1, linetype = "dashed",
             colour = "grey40", linewidth = 0.7) +
  scale_colour_manual(values = c("TRUE"  = viridis(1, begin = 0.2),
                                  "FALSE" = "grey60"),
                       labels = c("TRUE"  = "p < 0.05",
                                  "FALSE" = "p \u2265 0.05")) +
  coord_flip() +
  labs(title = "Univariate NB GLM — Incidence Rate Ratios",
       subtitle = "Each model includes livestock offset; dashed line at IRR = 1",
       x = NULL, y = "IRR (95% CI)", colour = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

p_daic <- ggplot(univariate_nb,
                 aes(x = reorder(variable, delta_AIC),
                     y = delta_AIC, fill = delta_AIC < 0)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 0, linewidth = 0.6, colour = "grey40") +
  scale_fill_manual(values = c("TRUE"  = viridis(1, begin = 0.5),
                                "FALSE" = "grey70"),
                    labels = c("TRUE"  = "Better than null",
                               "FALSE" = "Worse than null")) +
  coord_flip() +
  labs(title = "\u0394AIC vs null model",
       subtitle = "More negative = stronger univariate association",
       x = NULL, y = "\u0394AIC", fill = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("03_preliminary_models/univariate_nb_glm.png", p_irr + p_daic,
       width = 12, height = 5, dpi = 300)

## COMBINED SUMMARY TABLE

summary_table <- tibble(variable = covariates) %>%
  mutate(
    spearman_r = as.numeric(cor_out_spearman)[match(variable,
                                                    rownames(cor_out_spearman))]
  ) %>%
  left_join(
    univariate_nb %>%
      mutate(variable = str_remove(variable, "_z")) %>%
      dplyr::select(variable, IRR, p_value, delta_AIC),
    by = "variable"
  ) %>%
  left_join(vif_df %>% dplyr::select(variable, VIF), by = "variable") %>%
  arrange(delta_AIC) %>%
  mutate(across(where(is.numeric), ~ round(., 4)))

print(summary_table)
write_csv(summary_table, "03_preliminary_models/covariate_summary.csv")

## INTERPRETATION NOTES

cat("
=============================================================================
INTERPRETATION NOTES
=============================================================================

1. OVERDISPERSION TEST:
   Significant result confirms NB family over Poisson. If not significant,
   a Poisson hurdle is simpler and defensible.

2. ZERO-INFLATION ASSESSMENT:
   Excess zeros beyond NB expectation directly motivates the hurdle or ZI
   component in brms. If NB already accounts for zeros adequately, a plain
   NB mixed model may suffice and is worth testing with LOO.

3. ICC:
   Values above ~0.05 indicate meaningful between-zone clustering and justify
   the random intercept. High ICC (>0.3) suggests zones differ substantially
   — consider whether random slopes on key covariates are warranted.

4. ACF PLOTS:
   Significant lag-1 autocorrelation supports ar(time = year, gr = zone, p = 1)
   in brms. Negligible ACF suggests the AR(1) term may not be needed and adds
   complexity without benefit — worth testing with LOO.

5. SPEARMAN vs PEARSON:
   Prefer Spearman for climate variables. Disagreements between the two
   suggest non-normality or outlier influence worth examining.

6. VIF:
   Collinearity between PDSI, soil moisture, and precipitation is common in
   East African climate data. If VIF > 5, consider dropping the weaker
   variable or combining indices. High VIF inflates posterior uncertainty
   in brms without biasing point estimates.

7. UNIVARIATE NB GLMs:
   IRR > 1 = more cases with increasing covariate; IRR < 1 = fewer.
   Delta AIC provides a relative importance ranking on the correct
   distributional scale. Use as a screening guide — retain covariates with
   strong a priori ecological justification regardless of univariate result.

8. BIVARIATE SCATTERPLOTS:
   LOESS curves showing threshold or non-linear behaviour (especially in
   precipitation or PDSI) suggest whether spline terms s() via the mgcv
   interface in brms would be more appropriate than linear terms.
=============================================================================
")
