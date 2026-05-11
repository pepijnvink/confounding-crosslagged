library(LaplacesDemon)
source("Functions/R/00_helper_functions_PopAnalysis.R")
library(tidyverse)
nu = 500
mean_denom <- nu-100-1
sigma_wis <- diag(100)
sigma_wis[lower.tri(sigma_wis)] <- 0
sigma_wis[upper.tri(sigma_wis)] <- 0
set.seed(42)
corr_mat <- rinvwishart(nu = nu, S = sigma_wis) %>% cov2cor()

set.seed(1234)
beta_cont_s1 <- rbind(c(runif(10, min = -0.03, max = -0.01),
                        runif(40, min = 0.01, max = 0.11)),   # effects on x
                      c(runif(40, min = 0.01, max = 0.11),
                        runif(10, min = -0.03, max = -0.01))) # effects on y
# ## effects of (50) binary confounders
beta_dich_s1 <- rbind(c(runif(10, min = -0.03, max = -0.01),
                        runif(40, min = 0.01, max = 0.11)),   # effects on x
                      c(runif(40, min = 0.01, max = 0.11),
                        runif(10, min = -0.03, max = -0.01))) # effects on y

BetaAll <- cbind(beta_cont_s1, beta_dich_s1)

# Change factors Previous Scenario 2 (change in different directions). Variables are multiplied by this factor at t=3
mult.mat_s2_cont <- rbind(matrix(rep(c(2, 0.5, 1, 1, 1, 2, 0.5, 1, 1, 1), times = 5), nrow = 1, byrow = T), matrix(rep(c(0.5, 2, 1, 1, 1, 2, 0.5, 1, 1, 1), times = 5), nrow = 1, byrow = T)) # multiplication matrix for continuous confounders (increase of effects)
mult.mat_s2_dich <- rbind(matrix(rep(c(2, 0.5, 1, 1, 1, 2, 0.5, 1, 1, 1), times = 5), nrow = 1, byrow = T), matrix(rep(c(0.5, 2, 1, 1, 1, 2, 0.5, 1, 1, 1), times = 5), nrow = 1, byrow = T)) # multiplication matrix for dichotomous confounders (increase of effects)

# Change factors for Previous Scenario 3 (change in same direction)
mult.mat_s3_cont <- rbind(matrix(rep(c(2, 2, 1, 1, 1, 2, 2, 1, 1, 1), times = 5), nrow = 1, byrow = T), matrix(rep(c(2, 2, 1, 1, 1, 2, 2, 1, 1, 1), times = 5), nrow = 1, byrow = T)) # multiplication matrix for continuous confounders (increase of effects)
mult.mat_s3_dich <- rbind(matrix(rep(c(1.5, 1.5, 1, 1, 1, 1.5, 1.5, 1, 1, 1), times = 5), nrow = 1, byrow = T), matrix(rep(c(2, 2, 1, 1, 1, 2, 2, 1, 1, 1), times = 5), nrow = 1, byrow = T)) # multiplication matrix for dichotomous confounders (increase of effects)

# compute theoretical correlation
# matrix across binary and continuous confounders ------------------------------
corr_mat_theo <- transform_mixed_corr(Sigma = corr_mat, n_cont = 50, n_bin = 50)


## model parameters for cross lagged strucuture --------------------------------
phix <- 0.2 # AR for x
phiy <- 0.3 # AR for y
gammax <- 0 # CL on x
gammay <- 0.1 # CL on y
## phi matrix (AR and CL)
phi <- matrix(c(phix, gammax,
                gammay, phiy),
              nrow = 2,
              byrow = T)


## Create multiplication factors
# Scenario 1: Stable -----------------------------------------------------------
multi_scenario1_List <- lapply(1:5, FUN = function(i) matrix(1, nrow = 2, ncol = 100))


# Scenario 2: Proportional Step Function ---------------------------------------
multi_scenario2_List <- lapply(1:5, FUN = function(i) rbind(rep(ifelse(i < 3, 1, 1.3), 100),
                                                            rep(ifelse(i < 3, 1, .6), 100)))

# Scenario 3: Proportional Change ----------------------------------------------
# proportional but chaotic change

multi_scenario3_List <- lapply(1:5, FUN = function(i) rbind(rep(1 + .03*i*(-1.2)^(i-1)-.03, 100),
                                                            rep(1 - .01*i*(-2)^(i-1)+.01, 100)))

###############----
# partially stable scenarios ---------------------------------------------------
# Scenario 4: Partially Stable Step Function -----------------------------------
# partial change at t = 3 and in opposite direction
tempOther <- cbind(mult.mat_s2_cont, mult.mat_s2_dich); tempOther[tempOther == 1] <- 0;
tempOther[tempOther > 1] <- 1; tempOther[tempOther < 1 & tempOther != 0] <- -1
temp_List <- lapply(1:5, FUN = function(i) multi_scenario2_List[[i]]-1)
multi_scenario4_List <- lapply(1:5, FUN = function(i){temp1 <- temp_List[[i]]
1 + tempOther*temp1})

# Scenario 5: Partially Stable Proportional Change -----------------------------
# same changing behavior as in scenario 4: some change in opposite direction
tempOther <- multi_scenario4_List[[5]]; tempOther[tempOther == 1] <- 0;
tempOther[tempOther > 1] <- 1; tempOther[tempOther < 1 & tempOther != 0] <- -1
temp_List <- lapply(1:5, FUN = function(i) multi_scenario3_List[[i]]-1)
multi_scenario5_List <- lapply(1:5, FUN = function(i){temp1 <- temp_List[[i]]
1 + tempOther*temp1})

# Scenario 6: Random Directions with general increasing trend ------------------
set.seed(1234)
multi_scenario6_List <- lapply(1:5, FUN = function(i) rbind(runif(n = 100, min = 0, max = 2.2),
                                                             runif(n = 100, min = 0, max = 1.1)))
multi_scenario6_List[[1]][,] <- 1

## Get Sigma matrices
SigmaTarget <- matrix(.5, 2, 2); diag(SigmaTarget) <- 1 ## observed (population) covariance matrix
ParameterValues_1 <- computeSigmaWithBurnin_linear(nT = 5,
                                                   BetaC_List = get_coefs(multi_scenario1_List, BetaAll),
                                                   CrossLagged_List = phi,    # if is not list, then all are equal
                                                   Sigma_C = corr_mat_theo,   # theoretical correlation matrix among Cs (including binary)
                                                   SigmaTarget = SigmaTarget, # Observed Covariance Among Xt,Yt
                                                   burnin = 49)

ParameterValues_2 <- computeSigmaWithBurnin_linear(nT = 5,
                                                   BetaC_List = get_coefs(multi_scenario2_List, BetaAll),
                                                   CrossLagged_List = phi,    # if is not list, then all are equal
                                                   Sigma_C = corr_mat_theo,   # theoretical correlation matrix among Cs (including binary)
                                                   SigmaTarget = SigmaTarget, # Observed Covariance Among Xt,Yt
                                                   burnin = 49)

ParameterValues_3 <- computeSigmaWithBurnin_linear(nT = 5,
                                                   BetaC_List = get_coefs(multi_scenario3_List, BetaAll),
                                                   CrossLagged_List = phi,    # if is not list, then all are equal
                                                   Sigma_C = corr_mat_theo,   # theoretical correlation matrix among Cs (including binary)
                                                   SigmaTarget = SigmaTarget, # Observed Covariance Among Xt,Yt
                                                   burnin = 49)

ParameterValues_4 <- computeSigmaWithBurnin_linear(nT = 5,
                                                   BetaC_List = get_coefs(multi_scenario4_List, BetaAll),
                                                   CrossLagged_List = phi,    # if is not list, then all are equal
                                                   Sigma_C = corr_mat_theo,   # theoretical correlation matrix among Cs (including binary)
                                                   SigmaTarget = SigmaTarget, # Observed Covariance Among Xt,Yt
                                                   burnin = 49)

ParameterValues_5 <- computeSigmaWithBurnin_linear(nT = 5,
                                                   BetaC_List = get_coefs(multi_scenario5_List, BetaAll),
                                                   CrossLagged_List = phi,    # if is not list, then all are equal
                                                   Sigma_C = corr_mat_theo,   # theoretical correlation matrix among Cs (including binary)
                                                   SigmaTarget = SigmaTarget, # Observed Covariance Among Xt,Yt
                                                   burnin = 49)

ParameterValues_6 <- computeSigmaWithBurnin_linear(nT = 5,
                                                   BetaC_List = get_coefs(multi_scenario6_List, BetaAll),
                                                   CrossLagged_List = phi,    # if is not list, then all are equal
                                                   Sigma_C = corr_mat_theo,   # theoretical correlation matrix among Cs (including binary)
                                                   SigmaTarget = SigmaTarget, # Observed Covariance Among Xt,Yt
                                                   burnin = 49)

sigmas <- list(ParameterValues_1,
                   ParameterValues_2,
                   ParameterValues_3,
                   ParameterValues_4,
                   ParameterValues_5,
                   ParameterValues_6)

saveRDS(sigmas, "params/main/sigmas.rds")

multi_mats <- list(
  multi_scenario1_List,
  multi_scenario2_List,
  multi_scenario3_List,
  multi_scenario4_List,
  multi_scenario5_List,
  multi_scenario6_List
)
saveRDS(multi_mats, "params/main/multi_mats.rds")

saveRDS(list(beta_cont = beta_cont_s1, beta_dich = beta_dich_s1), "params/main/betaAll.rds")

