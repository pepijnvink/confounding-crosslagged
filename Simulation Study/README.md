# README

Code repository for *Methods to Account for Unobserved Baseline Confounders with Time-Invariant and Time-Varying Effects in Cross-Lagged Panel Research*.

Open a terminal in the current directory. Then run `bash run_all.sh` in the terminal to run all script files and obtain results and visualizations for the simulation study. The folder contains code for three studies:
- `main`: the main analyses in the manuscript.
- `bothdirec`: supplemental analyses with cross-lagged effects in both directions.
- `highconf`: supplemental analyses with a higher amount of confounding

## Folders after running

| Folder | Brief description |
|---|---|
| `figures` | After running the script, contains figures for each of the three studies.
| `Functions` | Contains functions necessary for simulation and analysis
| `output` | After running the script, contains simulated data, analyzed models, and results.
| `params` | Contains parameter settings for each study.
| `scripts` | Contains the scripts to run.

These folders are created automatically by `run_scripts.R` if they do not already exist.

## Simulation Scripts

For each study, the scripts to run are contained in `scripts` in a seperate folder (e.g. scripts/main/...). Each folder contains the following R files
| File | Brief description |
|---|---|
| `00_*_multiplication_tables.R` | Create parameter values for the study.
| `01_*_simulate.R` | Simulate data
| `02_*_analyses.R` | Analyze simulated data and obtain results
| `03_*_plot_results.R` | Create visualizations of results
| `04_*_plot_power.R` | Create visualizations, specifically for power (or Type II error).