argv = commandArgs(trailingOnly=TRUE)
# load packages
# Package names
packages <- c("ggplot2",
              "dplyr",
              "tidyr",
              "purrr",
              "cli",
              "Rcpp",
              "RcppArmadillo",
              "RcppProgress",
              "lavaan",
              "MASS",
              "furrr",
              "broom",
              "matrixcalc",
              "rlang",
              "readr",
              "furrr",
              "future",
              "jtools",
              "progressr",
              "LaplacesDemon")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages], repos = "http://cran.us.r-project.org")
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))
rm(list = c("installed_packages", "packages")) # remove objects

# Load necessary functions
## R
source("Functions/R/sim_stable_R.R") # function to simulate data with confounders with time-invariant effects
source("Functions/R/sim_varying_R.R") # function to simulate data with confounders with time-varying effects

## Load parameters
psi_all <- readRDS("params/sigmas.rds")
multi_mats <- readRDS("params/multi_mats.rds")
beta_all <- readRDS("params/betaAll.rds")

# Specify general hyperparameters
# N <- 500 # sample size
N <- as.numeric(argv[1])
ndat <- as.numeric(argv[2]) # number of datasets to simulate
burnin <- 50 # burn-in period
t <- 5 # number of timepoints to keep

## model parameters
phix <- 0.2 # AR for x
phiy <- 0.3 # AR for y
gammax <- 0.15 # CL on x
gammay <- 0.1 # CL on y
## phi matrix (AR and CL)
phi <- matrix(c(phix, gammax,
                gammay, phiy),
              nrow = 2,
              byrow = T)
## Covariance matrix of confounders
nu = 500
mean_denom <- nu-100-1
sigma_wis <- diag(100)
sigma_wis[lower.tri(sigma_wis)] <- 0
sigma_wis[upper.tri(sigma_wis)] <- 0
set.seed(42)
corr_mat <- rinvwishart(nu = nu, S = sigma_wis) %>% cov2cor()

# Scenario 1: Stable
## simulate data for scenario 1
cat("Simulating Scenario 1\n")
data_s1 <- sim_stable_R(timepoints = t,
                      burnin = burnin,
                      N = N,
                      ndat = ndat,
                      phi = phi,
                      betac_cont = beta_all$beta_cont,
                      betac_dich = beta_all$beta_dich,
                      psi = psi_all[[1]]$Sigma_Eps_List,
                      intercepts = c(0, 0),
                      meansc = 0,
                      sigmaC = corr_mat,
                      seed = 42)
saveRDS(data_s1, paste0("output/bothdirec/simdat/N", N, "/scenario1_N", N, ".rds"))
rm(data_s1)
# Scenario 2: Proportional Step
cat("Simulating Scenario 2\n")
data_s2 <- sim_varying_R(timepoints = t,
                         burnin = burnin,
                         N = N,
                         ndat = ndat,
                         phi = phi,
                         betac_cont = beta_all$beta_cont,
                         betac_dich = beta_all$beta_dich,
                         mult_mat = multi_mats[[2]],
                         psi = psi_all[[2]]$Sigma_Eps_List,
                         intercepts = c(0, 0),
                         meansc = 0,
                         sigmaC = corr_mat,
                         seed = 42)
saveRDS(data_s2, paste0("output/bothdirec/simdat/N", N, "/scenario2_N", N, ".rds"))
rm(data_s2)
# Scenario 3: Proportional Change
cat("Simulating Scenario 3\n")
data_s3 <- sim_varying_R(timepoints = t,
                         burnin = burnin,
                         N = N,
                         ndat = ndat,
                         phi = phi,
                         betac_cont = beta_all$beta_cont,
                         betac_dich = beta_all$beta_dich,
                         mult_mat = multi_mats[[3]],
                         psi = psi_all[[3]]$Sigma_Eps_List,
                         intercepts = c(0, 0),
                         meansc = 0,
                         sigmaC = corr_mat,
                         seed = 42)
saveRDS(data_s3, paste0("output/bothdirec/simdat/N", N, "/scenario3_N", N, ".rds"))
rm(data_s3)
# Scenario 4: Partially Stable Step
cat("Simulating Scenario 4\n")
data_s4 <- sim_varying_R(timepoints = t,
                         burnin = burnin,
                         N = N,
                         ndat = ndat,
                         phi = phi,
                         betac_cont = beta_all$beta_cont,
                         betac_dich = beta_all$beta_dich,
                         mult_mat = multi_mats[[4]],
                         psi = psi_all[[4]]$Sigma_Eps_List,
                         intercepts = c(0, 0),
                         meansc = 0,
                         sigmaC = corr_mat,
                         seed = 42)
saveRDS(data_s4, paste0("output/bothdirec/simdat/N", N, "/scenario4_N", N, ".rds"))
rm(data_s4)
# Scenario 5: Partially Stable Proportional Change
cat("Simulating Scenario 5\n")
data_s5 <- sim_varying_R(timepoints = t,
                         burnin = burnin,
                         N = N,
                         ndat = ndat,
                         phi = phi,
                         betac_cont = beta_all$beta_cont,
                         betac_dich = beta_all$beta_dich,
                         mult_mat = multi_mats[[5]],
                         psi = psi_all[[5]]$Sigma_Eps_List,
                         intercepts = c(0, 0),
                         meansc = 0,
                         sigmaC = corr_mat,
                         seed = 42)
saveRDS(data_s5, paste0("output/bothdirec/simdat/N", N, "/scenario5_N", N, ".rds"))
rm(data_s5)
# Scenario 6: Random directions with general increasing trend
cat("Simulating Scenario 6\n")
data_s6 <- sim_varying_R(timepoints = t,
                         burnin = burnin,
                         N = N,
                         ndat = ndat,
                         phi = phi,
                         betac_cont = beta_all$beta_cont,
                         betac_dich = beta_all$beta_dich,
                         mult_mat = multi_mats[[6]],
                         psi = psi_all[[6]]$Sigma_Eps_List,
                         intercepts = c(0, 0),
                         meansc = 0,
                         sigmaC = corr_mat,
                         seed = 42)
saveRDS(data_s6, paste0("output/bothdirec/simdat/N", N, "/scenario6_N", N, ".rds"))
rm(data_s6)
