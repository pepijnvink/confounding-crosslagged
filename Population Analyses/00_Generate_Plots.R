###############################################################################-
#### Supplement: Methods to Account for Unobserved Baseline Confounders with 
####             Time-Invariant and Time-Varying Effects in Cross-Lagged Panel 
####             Research
####
#### Script: Population analysis: Plot results
####

### plot main results (no Y->X, large confounding) -----------------------------
load("00_Data/01_PopulationAnalysis_NoYX_Large_C.rda")
Prefix <- "1_NoYX_Large_C_"
Add_Selection_Plots <- TRUE # add on to differentiate between models
removeAuto <- FALSE; removeZeroCross <- FALSE
source("00_helper_Plotting.R") |> suppressMessages() |> suppressWarnings()
load("00_Data/01_PopulationAnalysis_NoYX_Large_C.rda"); Prefix <- "1a_NoYX_Large_C_"; removeAuto <- TRUE; removeZeroCross <- TRUE
source("00_helper_Plotting.R") |> suppressMessages() |> suppressWarnings()
load("00_Data/01_PopulationAnalysis_NoYX_Large_C.rda"); Prefix <- "1b_NoYX_Large_C_"; removeAuto <- TRUE; removeZeroCross <- FALSE
source("00_helper_Plotting.R") |> suppressMessages() |> suppressWarnings()

### plot additional results (no Y->X, reduced confounding) ---------------------
load("00_Data/02_PopulationAnalysis_NoYX_Reduced.rda"); Prefix <- "2_NoYX_reduced_"; Add_Selection_Plots <- TRUE # add on to differentiate between models
removeAuto <- FALSE; removeZeroCross <- FALSE
source("00_helper_Plotting.R") |> suppressMessages() |> suppressWarnings()
load("00_Data/02_PopulationAnalysis_NoYX_Reduced.rda"); Prefix <- "2a_NoYX_reduced_"; removeAuto <- TRUE; removeZeroCross <- TRUE
source("00_helper_Plotting.R") |> suppressMessages() |> suppressWarnings()
load("00_Data/02_PopulationAnalysis_NoYX_Reduced.rda"); Prefix <- "2b_NoYX_reduced_"; removeAuto <- TRUE; removeZeroCross <- FALSE
source("00_helper_Plotting.R") |> suppressMessages() |> suppressWarnings()



### plot additional results (reciprocal, reduced confounding) ------------------
removeZeroCross <- FALSE # nothing to remove
load("00_Data/03_PopulationAnalysis_Reciprocal_Reduced.rda"); Prefix <- "3_reciprocal_reduced"; Add_Selection_Plots <- TRUE # add on to differentiate between models
removeAuto <- FALSE
source("00_helper_Plotting.R") |> suppressMessages() |> suppressWarnings()
load("00_Data/03_PopulationAnalysis_Reciprocal_Reduced.rda"); Prefix <- "3a_reciprocal_reduced"; removeAuto <- TRUE
source("00_helper_Plotting.R") |> suppressMessages() |> suppressWarnings()


