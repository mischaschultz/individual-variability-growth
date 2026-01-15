# ==============================================================================
# Back-Calculation of Size-at-Age Using Modified Fry Method
# ==============================================================================
# Purpose: Estimate individual fish length at each age using otolith increments
#          and a non-linear body-otolith relationship (Modified Fry method)
#
# Input:  Tidy otolith dataset from 01_TidyData.R
# Output: Dataset with back-calculated lengths (Li) at each age (agei)
#
# Key Parameters:
#   L0p: Biological intercept length (2.6 mm)
#   R0p: Otolith radius at hatching (0.013 mm = 13 microns)
# ==============================================================================

library(tidyverse)
library(nlme)

# Load Tidy Data ---------------------------------------------------------------
# Use growth_tidy from 01_TidyData.R or load csv
# growth_tidy <- read_csv("data/greater_amberjack_growth_data_tidy.csv")


# Define Biological Constants --------------------------------------------------
L0p <- 2.6   # Biological intercept length (mm)
R0p <- 0.013 # Otolith radius at hatching (mm)

# Estimate Body-Otolith Relationship Parameters --------------------------------
# Fit non-linear model: lencap = L0p - b * R0p^c + b * radcap^c
fit1 <- nls(
  lencap ~ L0p - b * R0p^c + b * radcap^c,
  data = growth_tidy,
  start = list(b = 1016.8, c = 0.81),
  control = list(maxiter = 500, warnOnly = TRUE)
)

summary(fit1)

# Compute Non-Linear Intercept (aNL) -------------------------------------------
coef_vals <- log(coef(fit1))
b <- coef_vals["b"]
c <- coef_vals["c"]
aNL <- L0p - b * R0p^c

# Back-Calculate Lengths at Each Age ------------------------------------------
# Apply Modified Fry equation to compute length at each increment
otolith_df <- growth_tidy %>%
  mutate(
    Li = aNL + exp(
      log(L0p - aNL) +
        ((log(lencap - aNL) - log(L0p - aNL)) / (log(radcap) - log(R0p))) *
        (log(radi) - log(R0p))
    )
  )

# Add Age-0 Hatchling Lengths from Known Distribution -------------------------
# Sample from empirical hatchling length distribution
set.seed(123)
n_draws <- 1000
hatchling_lengths <- rnorm(n_draws, mean = 3.085, sd = 0.132)

# Create age-0 rows for each individual
age0_rows <- otolith_df %>%
  group_by(id) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    agei = 0,
    Li = sample(hatchling_lengths, size = n(), replace = TRUE)
  )

# Combine with back-calculated data
otolith_data_with_age0 <- bind_rows(otolith_df, age0_rows) %>%
  arrange(id, agei)

# Export Back-Calculated Dataset -----------------------------------------------
write.csv(
  otolith_data_with_age0,
  "data/back_calculated_lengths.csv",
  row.names = FALSE
)