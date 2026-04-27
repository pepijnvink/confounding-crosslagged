# README

Code repository for *Methods to Account for Unobserved Baseline Confounders with Time-Invariant and Time-Varying Effects in Cross-Lagged Panel Research*.

Run `00_Source_All.R` to reproduce the full workflow for the population level analyses and simulation descriptives.
The file order below follows `00_Source_All.R`.

## Population Level Analyses, Plotting, and Descriptives

| File | Brief description |
|---|---|
| `00_Source_All.R` | Main wrapper script; creates output folders and runs all analyses, plots, and tables. |
| `00_Describe_Pop_Values.R` | Describes the inverse-Wishart correlation distributions used in the simulation setup. |
| `00_helper_functions_PopAnalysis.R` | Core helper functions for simulation setup, model fitting, and table generation. |
| `01_Population Analysis_NoYX_Large_C.R` | Main population analysis for the no-`Y -> X` condition with large confounding. |
| `02_Population Analysis_NoYX_Reduced.R` | Additional population analysis for the no-`Y -> X` condition with reduced confounding. |
| `03_Population Analysis_Reciprocal_Reduced.R` | Additional population analysis with reciprocal cross-lagged effects and reduced confounding. |
| `00_Generate_Plots.R` | Loads saved results and creates manuscript and supplemental figures. |
| `00_helper_Plotting.R` | Helper functions and plotting code for bias, relative bias, R² figures and more. |
| `00_Generate_Tables.R` | Generates the R² summary tables from the saved population analysis results. |

### Output folders

The scripts save results to:

- `00_Data/`
- `00_Data/PlotData/`
- `00_Figures/`
- `00_Tables/`

These folders are created automatically by `00_Source_All.R` if they do not already exist.

`00_Figures` contains the main manuscript and supplemental figures, while `00_Tables` contains the R² summary tables for the population analyses. 
The `00_Data` folder contains intermediate data files used for plotting and table generation.