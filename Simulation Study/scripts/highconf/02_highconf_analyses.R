argv = commandArgs(trailingOnly = TRUE)
N = as.numeric(argv[1])
oldw <- getOption("warn")
options(warn = -1)
# Load packages
## Package names
packages <- c("ggplot2",
              "dplyr",
              "tidyr",
              "purrr",
              "cli",
              "progressr",
              "lavaan",
              "MASS",
              "furrr",
              "broom",
              "matrixcalc",
              "rlang",
              "readr",
              "jtools")

## Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages], repos = "http://cran.us.r-project.org")
}

## Load Packages Not Yet Installed
invisible(lapply(packages, library, character.only = TRUE))
rm(list = c("installed_packages", "packages")) # remove objects

##options(future.globals.maxSize = 3000 * 1024^2) # Increase memory allocation

# Load function
source("Functions/R/analyze_data.R") # function to analyze data
source("Functions/R/check_convergence.R") # function to check convergence
source("Functions/R/remove_invalid.R") # function to check convergence
source("Functions/R/output_measures.R") # function to output measures (bias, mse etc.)

# Load models
riclpm <- read_lines("Models/riclpm.txt") # RI-CLPM
clpm <- read_lines("Models/clpm.txt") # CLPM
clpm_lag2 <- read_lines("Models/clpm2.txt") # CLPM with lag 2 effects
riclpm_free <- read_lines("Models/riclpm_free.txt") # RI-CLPM with residuals
dpm <- read_lines("Models/dpm.txt") # DPM
dpm_free <-read_lines("Models/dpm_free.txt") # DPM with residuals

# Object with variable names (used by `output_measures()`)
output_vars <- list(riclpm = list(xvar = "wx",
                                  yvar = "wy"),
                    clpm = list(xvar = "x",
                                yvar = "y"),
                    clpm_lag2 = list(xvar = "x",
                                     yvar = "y"),
                    riclpm_free = list(xvar = "wx",
                                      yvar = "wy"),
                    dpm = list(xvar = "x",
                               yvar = "y"),
                    dpm_free = list(xvar = "x",
                                   yvar = "y"))


# Scenario 1
## Choose confounders for propensity scores

##progressr::handlers("cli") # progress bar type
data_s1 <- readRDS(paste0("output/highconf/simdat/N", N, "/scenario1_N", N, ".rds"))
data_s1 <- analyze_data(data = data_s1,#analyze data
                        models = list(clpm = clpm,
                                      clpm_lag2 = clpm_lag2,
                                      riclpm = riclpm,
                                      riclpm_free = riclpm_free,
                                      dpm = dpm,
                                      dpm_free = dpm_free))
#plan(sequential) # back to sequential processing

data_s1 <- check_convergence(data_s1, # check convergence of models
                             latent = c(F, F, T, T, T, T),
                             filter = T)
data_s1 <- remove_invalid(data_s1) # remove invalid models
saveRDS(data_s1, file = paste0("output/highconf/analyzed_dat/N", N, "/analyzed_scenario1_N", N, ".rds")) # analyzed data
cat("Convergence and valid models of scenario 1:\n")
map(data_s1$convergence, function(x) summary(dplyr::select(x, -dataset_id))) # summary of convergence (i.e. how many of all models converged)
#data_s1$valid_mods_index$dpm_free <- data_s1$valid_mods_index$dpm_free[-c(253, 290)] # remove outliers
res_s1 <- output_measures(data_s1, xyvar = output_vars, filtered = T) # output measures (bias, mse etc.)
saveRDS(res_s1, file = paste0("output/highconf/results/N", N, "/results_scenario1_N", N, ".rds")) # save results
rm(data_s1, res_s1)

# Scenario 2
data_s2 <- readRDS(paste0("output/highconf/simdat/N", N, "/scenario2_N", N, ".rds"))
data_s2 <- analyze_data(data_s2, models = list(clpm = clpm, # analyze data
                                               clpm_lag2 = clpm_lag2,
                                               riclpm = riclpm,
                                               riclpm_free = riclpm_free,
                                               dpm = dpm,
                                               dpm_free = dpm_free))
data_s2 <- check_convergence(data_s2, latent = c(F, F, T, T, T, T), filter = T) # check convergence of models
data_s2 <- remove_invalid(data_s2) # remove invalid models
saveRDS(data_s2, file = paste0("output/highconf/analyzed_dat/N", N, "/analyzed_scenario2_N", N, ".rds")) # analyzed data
cat("Convergence and valid models of scenario 2:\n")
map(data_s2$convergence, function(x) summary(dplyr::select(x, -dataset_id)))
res_s2 <- output_measures(data_s2, xyvar = output_vars) # output measures (bias, mse etc.)
saveRDS(res_s2, file = paste0("output/highconf/results/N", N, "/results_scenario2_N", N, ".rds")) # save results
rm(data_s2)

# Scenario 3
data_s3 <- readRDS(paste0("output/highconf/simdat/N", N, "/scenario3_N", N, ".rds"))
data_s3 <- analyze_data(data_s3, models = list(clpm = clpm, # analyze data
                                               clpm_lag2 = clpm_lag2,
                                               riclpm = riclpm,
                                               riclpm_free = riclpm_free,
                                               dpm = dpm,
                                               dpm_free = dpm_free))
#plan(sequential) # back to sequential processing
data_s3 <- check_convergence(data_s3, latent = c(F, F, T, T, T, T)) # check convergence of models
data_s3 <- remove_invalid(data_s3) # remove invalid models
saveRDS(data_s3, file = paste0("output/highconf/analyzed_dat/N", N, "/analyzed_scenario3_N", N, ".rds")) # analyzed data
res_s3 <- output_measures(data_s3, xyvar = output_vars) # output measures (bias, mse etc.)
cat("Convergence and valid models of scenario 3:\n")
map(data_s3$convergence, function(x) summary(dplyr::select(x, -dataset_id)))
saveRDS(res_s3, file = paste0("output/highconf/results/N", N, "/results_scenario3_N", N, ".rds")) # save results
rm(data_s3)


# Scenario 4
data_s4 <- readRDS(paste0("output/highconf/simdat/N", N, "/scenario4_N", N, ".rds"))
data_s4 <- analyze_data(data_s4, models = list(clpm = clpm, # analyze data
                                               clpm_lag2 = clpm_lag2,
                                               riclpm = riclpm,
                                               riclpm_free = riclpm_free,
                                               dpm = dpm,
                                               dpm_free = dpm_free))
data_s4 <- check_convergence(data_s4, latent = c(F, F, T, T, T, T)) # check convergence of models
data_s4 <- remove_invalid(data_s4) # remove invalid models
saveRDS(data_s4, file = paste0("output/highconf/analyzed_dat/N", N, "/analyzed_scenario4_N", N, ".rds")) # analyzed data
cat("Convergence and valid models of scenario 4:\n")
map(data_s4$convergence, function(x) summary(dplyr::select(x, -dataset_id)))
res_s4 <- output_measures(data_s4, xyvar = output_vars) # output measures (bias, mse etc.)
saveRDS(res_s4, file = paste0("output/highconf/results/N", N, "/results_scenario4_N", N, ".rds")) # save results
rm(data_s4, res_s4)

# Scenario 5
data_s5 <- readRDS(paste0("output/highconf/simdat/N", N, "/scenario5_N", N, ".rds"))
data_s5 <- analyze_data(data_s5, models = list(clpm = clpm, # analyze data
                                               clpm_lag2 = clpm_lag2,
                                               riclpm = riclpm,
                                               riclpm_free = riclpm_free,
                                               dpm = dpm,
                                               dpm_free = dpm_free))
data_s5 <- check_convergence(data_s5, latent = c(F, F, T, T, T, T)) # check convergence of models
data_s5 <- remove_invalid(data_s5) # remove invalid models
saveRDS(data_s5, file = paste0("output/highconf/analyzed_dat/N", N, "/analyzed_scenario5_N", N, ".rds")) # analyzed data
cat("Convergence and valid models of scenario 5:\n")
map(data_s5$convergence, function(x) summary(dplyr::select(x, -dataset_id)))
res_s5 <- output_measures(data_s5, xyvar = output_vars) # output measures (bias, mse etc.)
saveRDS(res_s5, file = paste0("output/highconf/results/N", N, "/results_scenario5_N", N, ".rds")) # save results
rm(data_s5, res_s5)

# Scenario 6
data_s6 <- readRDS(paste0("output/highconf/simdat/N", N, "/scenario6_N", N, ".rds"))
data_s6 <- analyze_data(data_s6, models = list(clpm = clpm, # analyze data
                                               clpm_lag2 = clpm_lag2,
                                               riclpm = riclpm,
                                               riclpm_free = riclpm_free,
                                               dpm = dpm,
                                               dpm_free = dpm_free))

data_s6 <- check_convergence(data_s6, latent = c(F, F, T, T, T, T)) # check convergence of models
data_s6 <- remove_invalid(data_s6) # remove invalid models
saveRDS(data_s6, file = paste0("output/highconf/analyzed_dat/N", N, "/analyzed_scenario6_N", N, ".rds")) # analyzed data
cat("Convergence and valid models of scenario 6:\n")
map(data_s6$convergence, function(x) summary(dplyr::select(x, -dataset_id)))
res_s6 <- output_measures(data_s6, xyvar = output_vars) # output measures (bias, mse etc.)
saveRDS(res_s6, file = paste0("output/highconf/results/N", N, "/results_scenario6_N", N, ".rds")) # save results
rm(data_s6, res_s6)
options(warn = oldw)
