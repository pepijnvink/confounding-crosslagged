###############################################################################-
#### Supplement: Methods to Account for Unobserved Baseline Confounders with 
####             Time-Invariant and Time-Varying Effects in Cross-Lagged Panel 
####             Research
####
#### Script:Create R2 tables from run population level analyses

# functions --------------------------------------------------------------------
source("00_helper_functions_PopAnalysis.R")

### generate R2 tables for main results (no Y->X, large confounding) -----------
load("00_Data/01_PopulationAnalysis_NoYX_Large_C.rda")
create_R2_table(dataR2,
                condition_name = "Large-Confounding",
                label_suffix = "large")

### generate R2 tables for  additional results (no Y->X, reduced confounding) ---
load("00_Data/02_PopulationAnalysis_NoYX_Reduced.rda")
create_R2_table(dataR2,
                condition_name = "Reduced-Confounding",
                label_suffix = "reduced")

### generate R2 tables for  additional results (reciprocal, reduced confounding) ---
load("00_Data/03_PopulationAnalysis_Reciprocal_Reduced.rda")
create_R2_table(dataR2,
                condition_name = "Reduced-Confounding (Reciprocal Cross-Lag)",
                label_suffix = "reduced_bothCL")
