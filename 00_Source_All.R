###############################################################################-
#### Supplement: Methods to Account for Unobserved Baseline Confounders with 
####             Time-Invariant and Time-Varying Effects in Cross-Lagged Panel 
####             Research
####
#### Script: Source all population level analyses and description plot files. 
####         Made for complete reproducability

# set wd -----------------------------------------------------------------------
rstudioapi::documentPath() |> dirname() |> setwd()

# check whether folders are present and if not create them ---------------------
if (!dir.exists("00_Figures")) dir.create("00_Figures")
if (!dir.exists("00_Data")) dir.create("00_Data")
if (!dir.exists("00_Data/PlotData")) dir.create("00_Data/PlotData")
if (!dir.exists("00_Tables")) dir.create("00_Tables")

# descriptive plots ------------------------------------------------------------
source("00_Describe_Pop_Values.R")

# Run population level analysis ------------------------------------------------
# main article condition: large confounding with no Y -> X paths
source("01_Population Analysis_NoYX_Large_C.R") 

# additional condition: reduced confounding with no Y -> X paths
source("02_Population Analysis_NoYX_Reduced.R")

# additional condition: reciprocal effects with reduced confounding
source("03_Population Analysis_Reciprocal_Reduced.R")

# Plot population level analyses -----------------------------------------------
source("00_Generate_Plots.R") 

# Generate R2-tables -----------------------------------------------------------
library(kableExtra)
source("00_Generate_Tables.R") 
