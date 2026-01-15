# ==============================================================================
# Data Preparation for Individual Growth Analysis
# ==============================================================================
# Purpose: Clean and transform otolith measurement data for back-calculation
#          of individual size-at-age in greater amberjack (Seriola dumerili)
#
# Input:  Raw CSV with wide-format age readings and measurements
# Output: Long-format tidy dataset ready for back-calculation analysis
#
# Methods based on: Morat et al. (2020) - Individual back-calculated size-at-age
#                   from otoliths in Pacific coral reef fish species
# ==============================================================================

# Load Required Packages -------------------------------------------------------
library(tidyverse)
library(janitor)
library(stringr)

# Load and Clean Data ----------------------------------------------------------
growth <- read_csv("individual_growth20250731.csv") %>%
  clean_names() %>%
  rename(lencap = fork_length_mm) %>%
  # Remove duplicate records based on unique ID
  distinct(id_merge, .keep_all = TRUE)

# Transform Data from Wide to Long Format -------------------------------------
# Convert age columns (age_1, age_2, etc.) into individual rows
growth_tidy <- growth %>%
  pivot_longer(
    cols = starts_with("age_"),
    names_to = "agei",
    names_prefix = "age_",
    values_to = "radi"
  ) %>%
  mutate(
    agei = as.integer(agei),
    agecap = visible_annuli,
    radcap = otolith_length,
    l0p = 2.9  # Biological intercept length (mm)
  ) %>%
  select(
    id = id_merge,
    agei,
    radi,
    agecap,
    radcap,
    lencap,
    l0p,
    sex
  ) %>%
  # Convert otolith measurements from microns to mm
  mutate(
    radcap = radcap / 1000,
    radi = radi / 1000
  ) %>%
  # Remove rows with zero or negative radius measurements
  filter(radi > 0)

# Remove Incomplete Records ----------------------------------------------------
growth_tidy <- growth_tidy %>%
  filter(
    !is.na(lencap),
    !is.na(radcap),
    !is.na(radi),
    !is.na(agei),
    !is.na(id),
    !is.na(l0p)
  )

# Export Cleaned Data ----------------------------------------------------------
write.csv(
  growth_tidy,
  "data/greater_amberjack_growth_data_tidy.csv",
  row.names = FALSE
)
