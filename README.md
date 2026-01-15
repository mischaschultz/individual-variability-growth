# individual-variability-growth
Bayesian hierarchical analysis of individual growth variation in greater amberjack using otolith back-calculation

## Workflow

This repository contains some code examples for estimating individual-level growth parameters from otolith increment data. The analysis uses back-calculation methods to reconstruct size-at-age histories and applies hierarchical Bayesian models to quantify individual variation in growth trajectories for greater amberjack in the Gulf of Mexico. The analysis follows a sequential three-step process:

### 1. Data Preparation (`01_TidyData.R`)
- Transforms raw otolith measurements from wide to long format
- Cleans and validates measurement data
- Prepares dataset for back-calculation analysis

**Input:** Raw CSV with otolith increment measurements  
**Output:** Tidy long-format dataset

### 2. Back-Calculation (`02_BackCalculation.R`)
- Estimates body-otolith relationship using non-linear regression
- Applies Modified Fry method to back-calculate length-at-age
- Incorporates hatchling length distribution for age-0 fish

**Input:** Tidy dataset from Step 1  
**Output:** Dataset with back-calculated lengths (Li) at each age

### 3. Growth Modeling with brms (`03_vgbf_jags.R`)
- Fits hierarchical von Bertalanffy growth functions using JAGS
- Estimates population and individual-level growth parameters (Linf, k, t0)
- Generates sex-specific growth models
- Visualizes individual growth trajectories
- Conducts pairwise Bayesian comparisons of individual growth parameters
- Quantifies extent of individual variation within sexes

**Input:** Back-calculated lengths from Step 2  
**Output:** JAGS models, comparison statistics, publication-quality figures

## Key Methods

**Back-Calculation Method:** Modified Fry (non-linear body-otolith relationship)

**Growth Model:** Von Bertalanffy Growth Function (VBGF)
```
L(t) = Linf × (1 - exp(-k × (t - t0)))
```

Where:
- **Linf** = Asymptotic length
- **k** = Growth rate coefficient  
- **t0** = Theoretical age at length zero

**Statistical Framework:** Hierarchical Bayesian modeling with:
- Population-level (fixed) effects
- Individual-level (random) effects
- Sex-specific parameterization

## Requirements

### R Packages

```r
# Data manipulation
tidyverse
janitor
stringr

# Statistical modeling
R2jags
nlme

# Visualization
ggplot2
ggdist
patchwork
ggridges
```

### JAGS Installation

- **JAGS:** Requires JAGS installation (see [JAGS homepage](https://mcmc-jags.sourceforge.io/))

## Usage

Run scripts sequentially in numerical order:

```r
# Step 1: Prepare data
source("01_TidyData.R")

# Step 2: Back-calculate lengths
source("02_BackCalculation.R")

# Step 3: Fit JAGS models and conduct comparisons
source("03_vgbf_jags.R")
```

## Data Format

### Input Data Structure

The raw data file should contain:
- `id_merge`: Unique fish identifier
- `fork_length_mm`: Measured fork length at capture
- `visible_annuli`: Age determined from otolith
- `otolith_length`: Maximum otolith radius (microns)
- `age_1`, `age_2`, ..., `age_n`: Otolith radius at each annulus (microns)
- `sex`: Fish sex (M/F)

## Key Results

The analysis produces:

1. **Individual growth curves** showing variation in growth trajectories
2. **Sex-specific growth parameters** with credible intervals
3. **Pairwise comparisons** quantifying the proportion of individuals with significantly different growth rates
4. **Publication-quality figures** of individual and population growth curves

## References

**Methodological Framework:**
- Morat, F., Letourneur, Y., Nervousness, M., Banaru, D., & Batjakas, I. E. (2020). Individual back-calculated size-at-age based on otoliths from Pacific coral reef fish species. *Scientific Data*, 7(1), 370.

**Growth Modeling:**
- von Bertalanffy, L. (1938). A quantitative theory of organic growth. *Human Biology*, 10(2), 181-213.

