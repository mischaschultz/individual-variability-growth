# ==============================================================================
# Von Bertalanffy Growth Function (VBGF) Analysis Using JAGS
# ==============================================================================
# Purpose: Fit hierarchical Bayesian VBGF models using JAGS to back-calculated 
#          size-at-age data to estimate individual and population level growth,
#          as well as conducting pairwise comparison of growth parameters
#
# Input:  Back-calculated length data from 02_BackCalculation.R
# Output: Sex-specific growth models and pairwise parameter comparisons
# ==============================================================================

library(tidyverse)
library(janitor)
library(R2jags)
library(ggdist)
library(patchwork)
library(ggridges)

# Set theme --------------------------------------------------------------------
theme_clean <- function() {
  theme_minimal() +
    theme(
      panel.grid = element_blank(),  # Remove all grid lines
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Black border around plot
      plot.margin = margin(10, 10, 10, 10),  # External margins outside the black box
      axis.title = element_text(face = "bold", margin = margin(t = 5, r = 5, b = 5, l = 5)),  # Small margin for clarity
      axis.text = element_text(margin = margin(0, 0, 0, 0)),  # Remove internal padding from axis text
      axis.ticks.length = unit(0, "pt"),  # Remove extra space from axis ticks
      panel.spacing = unit(1, "lines"),  # No extra spacing inside the panel
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold", size = rel(1), hjust = .5),
      strip.background = element_rect(fill = "grey90", color = "black"),
      legend.title = element_text(face = "bold"),
      legend.position = "right"  # Place legend on the side
    )
}

# Load and Prepare Data --------------------------------------------------------
# Use otolith_df from 02_BackCalculation.R or load csv
otolith_df <- otolith_df %>%
  mutate(
    id = sub("\\.tif$", "", id),
    id = factor(trimws(id)),
    Li_cm = Li / 10
  ) %>%
  filter(agecap > 1) %>%
  droplevels()

# ==============================================================================
# POOLED MODEL (Both Sexes)
# ==============================================================================

# Define JAGS Model ------------------------------------------------------------
sink("vonBmodel_hierarchical.txt")
cat("
model {
  # Priors for population-level means
  Linf_mu ~ dnorm(linf_prior, 0.001)
  K_mu ~ dunif(0, 0.5)
  t0_mu ~ dnorm(-0.81, 1)

  # Priors for standard deviations (hierarchical)
  sigma.y ~ dt(0, pow(10, -2), 3) T(0,)
  sigma.Linf ~ dt(0, pow(10, -2), 3) T(0,)
  sigma.K ~ dt(0, pow(10, -2), 3) T(0,)
  sigma.t0 ~ dt(0, pow(10, -2), 3) T(0,)

  # Random effects per individual
  for (j in 1:N_1) {
    z_Linf[j] ~ dnorm(0, 1)
    z_K[j] ~ dnorm(0, 1)
    z_t0[j] ~ dnorm(0, 1)

    Linf[j] <- Linf_mu + z_Linf[j] * sigma.Linf
    K[j] <- K_mu + z_K[j] * sigma.K
    t0[j] <- t0_mu + z_t0[j] * sigma.t0
  }

  # Likelihood
  for (i in 1:N) {
    mu[i] <- Linf[J[i]] * (1 - exp(-K[J[i]] * (x[i] - t0[J[i]])))
    y[i] ~ dnorm(mu[i], pow(sigma.y, -2))
  }
}
", fill = TRUE)
sink()

# Bundle Data for JAGS ---------------------------------------------------------
data_jags <- list(
  N = length(otolith_df$Li),
  y = otolith_df$Li,
  x = otolith_df$agei,
  N_1 = length(unique(otolith_df$id)),
  J = as.numeric(factor(otolith_df$id)),
  linf_prior = 1100
)

# Initial Values ---------------------------------------------------------------
inits <- function() {
  list(
    Linf_mu = 1100,
    K_mu = 0.2,
    t0_mu = -1,
    sigma.y = 10,
    sigma.Linf = 50,
    sigma.K = 0.1,
    sigma.t0 = 0.5,
    z_Linf = rnorm(data_jags$N_1, 0, 1),
    z_K = rnorm(data_jags$N_1, 0, 1),
    z_t0 = rnorm(data_jags$N_1, 0, 1)
  )
}

# MCMC Settings ----------------------------------------------------------------
ni <- 500000   # Iterations
nt <- 2      # Thinning
nb <- 100000   # Burn-in
nc <- 4      # Chains

# Parameters to Save -----------------------------------------------------------
parameters <- c(
  "Linf_mu", "K_mu", "t0_mu",
  "sigma.y", "sigma.Linf", "sigma.K", "sigma.t0",
  "Linf", "K", "t0"
)

# Run Model --------------------------------------------------------------------
out_vb <- jags(
  data = data_jags,
  inits = inits,
  parameters.to.save = parameters,
  model.file = "vonBmodel_hierarchical.txt",
  n.chains = nc,
  n.thin = nt,
  n.iter = ni,
  n.burnin = nb
)

print(out_vb, digits = 2)

# ==============================================================================
# MALE-SPECIFIC MODEL
# ==============================================================================

# Define JAGS Model for Males --------------------------------------------------
sink("vonBmodel_hierarchicalMales.txt")
cat("
model {
  # Priors for population-level means
  Linf_mu ~ dnorm(linf_prior, 0.0001)
  K_mu ~ dunif(0, 0.6)

  # Priors for standard deviations
  sigma.y ~ dt(0, pow(10, -2), 3) T(0,)
  sigma.Linf ~ dt(0, pow(10, -2), 3) T(0,)
  sigma.K ~ dt(0, pow(10, -2), 3) T(0,)

  # Random effects per individual
  for (j in 1:N_1) {
    z_Linf[j] ~ dnorm(0, 1)
    z_K[j] ~ dnorm(0, 1)

    Linf[j] <- Linf_mu + z_Linf[j] * sigma.Linf
    K[j] <- K_mu + z_K[j] * sigma.K
  }

  # Likelihood (fixed t0)
  for (i in 1:N) {
    mu[i] <- Linf[J[i]] * (1 - exp(-K[J[i]] * (x[i] - t0)))
    y[i] ~ dnorm(mu[i], pow(sigma.y, -2))
  }
}
", fill = TRUE)
sink()

# Prepare Male Data ------------------------------------------------------------
otolith_male <- otolith_df %>% filter(sex == "M")

data_jags_male <- list(
  N = length(otolith_male$Li),
  y = otolith_male$Li,
  x = otolith_male$agei,
  N_1 = length(unique(otolith_male$id)),
  J = as.numeric(factor(otolith_male$id)),
  linf_prior = 1100,
  t0 = -0.9
)

# Initial Values for Males -----------------------------------------------------
inits_male <- function() {
  list(
    Linf_mu = 1100,
    K_mu = 0.2,
    sigma.y = 10,
    sigma.Linf = 50,
    sigma.K = 0.1,
    z_Linf = rnorm(data_jags_male$N_1, 0, 1),
    z_K = rnorm(data_jags_male$N_1, 0, 1)
  )
}

# Run Male Model ---------------------------------------------------------------
out_vb_male <- jags(
  data = data_jags_male,
  inits = inits_male,
  parameters.to.save = c("Linf_mu", "K_mu", "sigma.y", "sigma.Linf", "sigma.K", "Linf", "K"),
  model.file = "vonBmodel_hierarchicalMales.txt",
  n.chains = nc,
  n.thin = nt,
  n.iter = 500000,
  n.burnin = 100000
)

print(out_vb_male, digits = 2)

# Visualize Male Growth Trajectories -------------------------------------------
mcmc_samples <- as.mcmc(out_vb_male)
posterior_mat <- as.matrix(mcmc_samples)
age_grid <- seq(0, max(data_jags_male$x), length.out = 100)
n_ind <- data_jags_male$N_1

# Extract parameters
linf_mu_samples <- posterior_mat[, "Linf_mu"]
k_mu_samples <- posterior_mat[, "K_mu"]
linf_ind_samples <- posterior_mat[, grep("^Linf\\[", colnames(posterior_mat))]
k_ind_samples <- posterior_mat[, grep("^K\\[", colnames(posterior_mat))]

# VBGF function
vbgf <- function(Linf, K, age, t0 = -0.9) {
  Linf * (1 - exp(-K * (age - t0)))
}

# Individual curves
ind_curves_male <- map_dfr(1:n_ind, function(i) {
  tibble(
    age = age_grid,
    length = vbgf(mean(linf_ind_samples[, i]), mean(k_ind_samples[, i]), age_grid),
    id = paste0("ind_", i),
    type = "individual"
  )
})

# Global curve
global_curve_male <- tibble(
  age = age_grid,
  length = vbgf(mean(linf_mu_samples), mean(k_mu_samples), age_grid),
  id = "global",
  type = "global"
)

growth_df_male <- bind_rows(ind_curves_male, global_curve_male)

# Plot
ind_male <- ggplot(growth_df_male, aes(x = age, y = length, group = id, color = type)) +
  geom_line(
    data = filter(growth_df_male, type == "individual"),
    size = 0.5, alpha = 0.6
  ) +
  geom_line(
    data = filter(growth_df_male, type == "global"),
    size = 1.5, alpha = 1
  ) +
  theme_minimal() +
  labs(title = "VBGF Male", x = "Age", y = "Length (mm)") +
  scale_color_manual(values = c("individual" = "gray70", "global" = "#5165a4"), guide = "none") +
  ylim(0, 1500) +
  scale_x_continuous(limits = c(0, 10), breaks = c(0:5, 7, 9)) +
  theme_clean()

# ==============================================================================
# FEMALE-SPECIFIC MODEL
# ==============================================================================

# Define JAGS Model for Females ------------------------------------------------
sink("vonBmodel_hierarchicalFemales.txt")
cat("
model {
  # Priors for population-level means
  Linf_mu ~ dnorm(linf_prior, 0.0001)
  K_mu ~ dunif(0, 0.6)

  # Priors for standard deviations
  sigma.y ~ dt(0, pow(10, -2), 3) T(0,)
  sigma.Linf ~ dt(0, pow(10, -2), 3) T(0,)
  sigma.K ~ dt(0, pow(10, -2), 3) T(0,)

  # Random effects per individual
  for (j in 1:N_1) {
    z_Linf[j] ~ dnorm(0, 1)
    z_K[j] ~ dnorm(0, 1)

    Linf[j] <- Linf_mu + z_Linf[j] * sigma.Linf
    K[j] <- K_mu + z_K[j] * sigma.K
  }

  # Likelihood (fixed t0)
  for (i in 1:N) {
    mu[i] <- Linf[J[i]] * (1 - exp(-K[J[i]] * (x[i] - t0)))
    y[i] ~ dnorm(mu[i], pow(sigma.y, -2))
  }
}
", fill = TRUE)
sink()

# Prepare Female Data ----------------------------------------------------------
otolith_female <- otolith_df %>% filter(sex == "F")

data_jags_female <- list(
  N = length(otolith_female$Li),
  y = otolith_female$Li,
  x = otolith_female$agei,
  N_1 = length(unique(otolith_female$id)),
  J = as.numeric(factor(otolith_female$id)),
  linf_prior = 1300,
  t0 = -0.9
)

# Initial Values for Females ---------------------------------------------------
inits_female <- function() {
  list(
    Linf_mu = 1100,
    K_mu = 0.2,
    sigma.y = 10,
    sigma.Linf = 50,
    sigma.K = 0.1,
    z_Linf = rnorm(data_jags_female$N_1, 0, 1),
    z_K = rnorm(data_jags_female$N_1, 0, 1)
  )
}

# Run Female Model -------------------------------------------------------------
out_vb_female <- jags(
  data = data_jags_female,
  inits = inits_female,
  parameters.to.save = c("Linf_mu", "K_mu", "sigma.y", "sigma.Linf", "sigma.K", "Linf", "K"),
  model.file = "vonBmodel_hierarchicalFemales.txt",
  n.chains = nc,
  n.thin = 3,
  n.iter = 600000,
  n.burnin = 100000
)

print(out_vb_female, digits = 2)

# Visualize Female Growth Trajectories -----------------------------------------
mcmc_samples_f <- as.mcmc(out_vb_female)
posterior_mat_f <- as.matrix(mcmc_samples_f)
n_ind_f <- data_jags_female$N_1

# Extract parameters
linf_mu_samples_f <- posterior_mat_f[, "Linf_mu"]
k_mu_samples_f <- posterior_mat_f[, "K_mu"]
linf_ind_samples_f <- posterior_mat_f[, grep("^Linf\\[", colnames(posterior_mat_f))]
k_ind_samples_f <- posterior_mat_f[, grep("^K\\[", colnames(posterior_mat_f))]

# Individual curves
ind_curves_female <- map_dfr(1:n_ind_f, function(i) {
  tibble(
    age = age_grid,
    length = vbgf(mean(linf_ind_samples_f[, i]), mean(k_ind_samples_f[, i]), age_grid),
    id = paste0("ind_", i),
    type = "individual"
  )
})

# Global curve
global_curve_female <- tibble(
  age = age_grid,
  length = vbgf(median(linf_mu_samples_f), median(k_mu_samples_f), age_grid),
  id = "global",
  type = "global"
)

growth_df_female <- bind_rows(ind_curves_female, global_curve_female)

# Plot
ind_female <- ggplot(growth_df_female, aes(x = age, y = length, group = id, color = type)) +
  geom_line(
    data = filter(growth_df_female, type == "individual"),
    size = 0.5, alpha = 0.6
  ) +
  geom_line(
    data = filter(growth_df_female, type == "global"),
    size = 1.5, alpha = 1
  ) +
  theme_minimal() +
  labs(title = "VBGF Female", x = "Age", y = "Length (mm)") +
  scale_color_manual(values = c("individual" = "gray70", "global" = "#fe735d"), guide = "none") +
  ylim(0, 1500) +
  scale_x_continuous(limits = c(0, 10), breaks = c(0:5, 8, 10)) +
  theme_clean()

# Combine Plots ----------------------------------------------------------------
ind_both <- ind_female + ind_male

ggsave(
  filename = "figures/vonB_Ind_Both.png",
  ind_both,
  width = 11, height = 6, dpi = 600, units = "in", device = "png", bg = "white"
)

# ==============================================================================
# CREDIBLE INTERVALS FOR PARAMETERS
# ==============================================================================

# Female Parameter Summaries ---------------------------------------------------
post_female <- out_vb_female$BUGSoutput$sims.matrix

# Extract and summarize Linf
linf_cols <- grep("^Linf\\[", colnames(post_female), value = TRUE)
linf_summary_f <- post_female[, linf_cols] %>%
  as_tibble() %>%
  pivot_longer(everything(), names_to = "param", values_to = "linf") %>%
  group_by(param) %>%
  summarise(
    linf_mean = mean(linf),
    linf_lower = quantile(linf, 0.025),
    linf_upper = quantile(linf, 0.975),
    .groups = "drop"
  ) %>%
  mutate(Individual = as_factor(1:n())) %>%
  mutate(Individual = fct_reorder(Individual, linf_mean))

# Extract and summarize K
k_cols <- grep("^K\\[", colnames(post_female), value = TRUE)
k_summary_f <- post_female[, k_cols] %>%
  as_tibble() %>%
  pivot_longer(everything(), names_to = "param", values_to = "k") %>%
  group_by(param) %>%
  summarise(
    k_mean = mean(k),
    k_lower = quantile(k, 0.025),
    k_upper = quantile(k, 0.975),
    .groups = "drop"
  ) %>%
  mutate(Individual = as_factor(1:n())) %>%
  mutate(Individual = fct_reorder(Individual, k_mean))

# Plot credible intervals
Linf_plot_f <- ggplot(linf_summary_f, aes(x = linf_mean, y = Individual)) +
  geom_pointinterval(aes(xmin = linf_lower, xmax = linf_upper)) +
  labs(x = "Linf (mm)", y = "Individual (Female)") +
  theme_minimal() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  theme_clean()

K_plot_f <- ggplot(k_summary_f, aes(x = k_mean, y = Individual)) +
  geom_pointinterval(aes(xmin = k_lower, xmax = k_upper)) +
  labs(x = "K", y = "Individual (Female)") +
  theme_minimal() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  theme_clean()

Linf_plot_f / K_plot_f

# Male Parameter Summaries -----------------------------------------------------
post_male <- out_vb_male$BUGSoutput$sims.matrix

# Extract and summarize Linf
linf_cols_m <- grep("^Linf\\[", colnames(post_male), value = TRUE)
linf_summary_m <- post_male[, linf_cols_m] %>%
  as_tibble() %>%
  pivot_longer(everything(), names_to = "param", values_to = "linf") %>%
  group_by(param) %>%
  summarise(
    linf_mean = mean(linf),
    linf_lower = quantile(linf, 0.025),
    linf_upper = quantile(linf, 0.975),
    .groups = "drop"
  ) %>%
  mutate(Individual = as_factor(1:n())) %>%
  mutate(Individual = fct_reorder(Individual, linf_mean))

# Extract and summarize K
k_cols_m <- grep("^K\\[", colnames(post_male), value = TRUE)
k_summary_m <- post_male[, k_cols_m] %>%
  as_tibble() %>%
  pivot_longer(everything(), names_to = "param", values_to = "k") %>%
  group_by(param) %>%
  summarise(
    k_mean = mean(k),
    k_lower = quantile(k, 0.025),
    k_upper = quantile(k, 0.975),
    .groups = "drop"
  ) %>%
  mutate(Individual = as_factor(1:n())) %>%
  mutate(Individual = fct_reorder(Individual, k_mean))

# Plot credible intervals
Linf_plot_m <- ggplot(linf_summary_m, aes(x = linf_mean, y = Individual)) +
  geom_pointinterval(aes(xmin = linf_lower, xmax = linf_upper)) +
  labs(x = "Linf (mm)", y = "Individual (Male)") +
  theme_minimal() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  theme_clean()

K_plot_m <- ggplot(k_summary_m, aes(x = k_mean, y = Individual)) +
  geom_pointinterval(aes(xmin = k_lower, xmax = k_upper)) +
  labs(x = "K", y = "Individual (Male)") +
  theme_minimal() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  theme_clean()

Linf_plot_m / K_plot_m

# ==============================================================================
# PAIRWISE COMPARISONS
# ==============================================================================

# Function to compute pairwise posterior probabilities -------------------------
compute_pairwise_probs <- function(param_matrix) {
  n_indiv <- ncol(param_matrix)
  diff_prob_mat <- matrix(NA, nrow = n_indiv, ncol = n_indiv)
  colnames(diff_prob_mat) <- rownames(diff_prob_mat) <- colnames(param_matrix)
  
  for (i in 1:n_indiv) {
    for (j in 1:n_indiv) {
      if (i != j) {
        diff <- param_matrix[, i] - param_matrix[, j]
        diff_prob_mat[i, j] <- mean(diff > 0)
      }
    }
  }
  
  # Identify significant pairs (> 97.5% or < 2.5%)
  signif_pairs <- which(diff_prob_mat > 0.975 | diff_prob_mat < 0.025, arr.ind = TRUE)
  
  signif_pairs_named <- data.frame(
    Indiv_A = rownames(diff_prob_mat)[signif_pairs[, 1]],
    Indiv_B = colnames(diff_prob_mat)[signif_pairs[, 2]],
    Prob_A_gt_B = diff_prob_mat[signif_pairs]
  )
  
  # Calculate percentage of significant pairs
  total_pairs <- n_indiv * (n_indiv - 1)
  n_signif <- nrow(signif_pairs_named)
  percent_signif <- 100 * n_signif / total_pairs
  
  list(
    prob_matrix = diff_prob_mat,
    signif_pairs = signif_pairs_named,
    percent_signif = percent_signif
  )
}

# Females: K comparisons -------------------------------------------------------
k_cols_f <- grep("^K\\[", colnames(post_female), value = TRUE)
K_mat_f <- post_female[, k_cols_f]
female_k_comparison <- compute_pairwise_probs(K_mat_f)

cat("\nFemales - K Parameter:\n")
cat("Percentage of significantly different pairs:", 
    round(female_k_comparison$percent_signif, 2), "%\n")

# Females: Linf comparisons ----------------------------------------------------
linf_cols_f <- grep("^Linf\\[", colnames(post_female), value = TRUE)
Linf_mat_f <- post_female[, linf_cols_f]
female_linf_comparison <- compute_pairwise_probs(Linf_mat_f)

cat("\nFemales - Linf Parameter:\n")
cat("Percentage of significantly different pairs:", 
    round(female_linf_comparison$percent_signif, 2), "%\n")

# Males: K comparisons ---------------------------------------------------------
k_cols_m <- grep("^K\\[", colnames(post_male), value = TRUE)
K_mat_m <- post_male[, k_cols_m]
male_k_comparison <- compute_pairwise_probs(K_mat_m)

cat("\nMales - K Parameter:\n")
cat("Percentage of significantly different pairs:", 
    round(male_k_comparison$percent_signif, 2), "%\n")

# Males: Linf comparisons ------------------------------------------------------
linf_cols_m <- grep("^Linf\\[", colnames(post_male), value = TRUE)
Linf_mat_m <- post_male[, linf_cols_m]
male_linf_comparison <- compute_pairwise_probs(Linf_mat_m)

cat("\nMales - Linf Parameter:\n")
cat("Percentage of significantly different pairs:", 
    round(male_linf_comparison$percent_signif, 2), "%\n")

# Save workspace ---------------------------------------------------------------

save.image("03_vgbf_jags.RData")

