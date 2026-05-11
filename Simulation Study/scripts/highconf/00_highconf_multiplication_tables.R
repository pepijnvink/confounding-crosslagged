library(LaplacesDemon)
source("Functions/R/00_helper_functions_PopAnalysis.R")
library(tidyverse)
nu = 500
mean_denom <- nu-100-1
sigma_wis <- diag(100)
sigma_wis[lower.tri(sigma_wis)] <- 0
sigma_wis[upper.tri(sigma_wis)] <- 0
set.seed(42)
cor_mat <- rinvwishart(nu = nu, S = sigma_wis) %>% cov2cor()

# compute theoretical correlation
# matrix across binary and continuous confounders ------------------------------
cov_mat_theo <- transform_mixed_corr(Sigma = cor_mat, n_cont = 50, n_bin = 50)


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
# Cov[X_t,Y_t]
SigmaTarget <- matrix(.5, 2, 2); diag(SigmaTarget) <- 1

### Effects of betas on X,Y ----------------------------------------------------
## Simulate effects on x and y for continuous and binary
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

# Change factors Scenario 2 (change in different directions). Variables are multiplied by this factor at t=3
mult.mat_s2_cont <- rbind(matrix(rep(c(2, 0.5, 1, 1, 1, 2, 0.5, 1, 1, 1), times = 5), nrow = 1, byrow = T), matrix(rep(c(0.5, 2, 1, 1, 1, 2, 0.5, 1, 1, 1), times = 5), nrow = 1, byrow = T)) # multiplication matrix for continuous confounders (increase of effects)
mult.mat_s2_dich <- rbind(matrix(rep(c(2, 0.5, 1, 1, 1, 2, 0.5, 1, 1, 1), times = 5), nrow = 1, byrow = T), matrix(rep(c(0.5, 2, 1, 1, 1, 2, 0.5, 1, 1, 1), times = 5), nrow = 1, byrow = T)) # multiplication matrix for dichotomous confounders (increase of effects)

# Change factors for Scenario 3 (change in same direction)
mult.mat_s3_cont <- rbind(matrix(rep(c(2, 2, 1, 1, 1, 2, 2, 1, 1, 1), times = 5), nrow = 1, byrow = T), matrix(rep(c(2, 2, 1, 1, 1, 2, 2, 1, 1, 1), times = 5), nrow = 1, byrow = T)) # multiplication matrix for continuous confounders (increase of effects)
mult.mat_s3_dich <- rbind(matrix(rep(c(1.5, 1.5, 1, 1, 1, 1.5, 1.5, 1, 1, 1), times = 5), nrow = 1, byrow = T), matrix(rep(c(2, 2, 1, 1, 1, 2, 2, 1, 1, 1), times = 5), nrow = 1, byrow = T)) # multiplication matrix for dichotomous confounders (increase of effects)

# For comparison
VarComp <- 4 # This is the variance given to the latent factors that should
# explain the variation of the confounders
# -> equal for comparsion (so that the chisq has a meaning!)
# Beta that are scaled
BetaAll <- cbind(beta_cont_s1, beta_dich_s1)*sqrt(2)

# effects should not be too small - otherwise not visible in plots
BetaAll[abs(BetaAll) < .002]

# explained variance (if X,Y are standardized)
diag(BetaAll %*% cov_mat_theo %*% t(BetaAll))
diag(phi %*% SigmaTarget %*% t(phi))

# Scenario 1: Stable -----------------------------------------------------------
multi_scenario1_List <- lapply(1:5, FUN = function(i) matrix(1, nrow = 2, ncol = 100))
StepOne_1 <- automate_scenario_cfa(BetaBase = BetaAll,
                                   cov_mat_theo = cov_mat_theo,
                                   VarComp = VarComp,
                                   multi_List = multi_scenario1_List)


# Scenario 2: Proportional Step Function ---------------------------------------
# jump from 1 to 2 at t = 3
multi_scenario2_List <- lapply(1:5, FUN = function(i) rbind(rep(ifelse(i < 3, 1, 1.1), 100),
                                                            rep(ifelse(i < 3, 1, .8), 100))) # edited
StepOne_2 <- automate_scenario_cfa(BetaBase = BetaAll,
                                   cov_mat_theo = cov_mat_theo,
                                   VarComp = VarComp,
                                   multi_List = multi_scenario2_List)


# Scenario 3: Proportional Change ----------------------------------------------
# proportional but chaotic change

multi_scenario3_List <- lapply(1:5, FUN = function(i) rbind(rep(1 + .005*i*(-1.15)^(i-1)-.005, 100),
                                                            rep(1 + .005*i*(-1.25)^(i-1)-.005, 100)))
StepOne_3 <- automate_scenario_cfa(BetaBase = BetaAll,
                                   cov_mat_theo = cov_mat_theo,
                                   VarComp = VarComp,
                                   multi_List = multi_scenario3_List)


######## From here on, freely estimated DPM should no longer work --------------

###############----
# partially stable scenarios ---------------------------------------------------

# Scenario 4: Partially Stable Step Function -----------------------------------
# partial change at t = 3 and in opposite direction
tempOther <- cbind(mult.mat_s2_cont, mult.mat_s2_dich); tempOther[tempOther == 1] <- 0;
tempOther[tempOther > 1] <- 1; tempOther[tempOther < 1 & tempOther != 0] <- -1
temp_List <- lapply(1:5, FUN = function(i) multi_scenario2_List[[i]]-1)
multi_scenario4_List <- lapply(1:5, FUN = function(i){temp1 <- temp_List[[i]]
1 + tempOther*temp1})
StepOne_4 <- automate_scenario_cfa(BetaBase = BetaAll,
                                   cov_mat_theo = cov_mat_theo,
                                   VarComp = VarComp,
                                   multi_List = multi_scenario4_List)

# Scenario 5: Partially Stable Proportional Change -----------------------------
# same changing behavior as in scenario 8: some change in opposite direction
tempOther <- multi_scenario4_List[[5]]; tempOther[tempOther == 1] <- 0;
tempOther[tempOther > 1] <- 1; tempOther[tempOther < 1 & tempOther != 0] <- -1
temp_List <- lapply(1:5, FUN = function(i) multi_scenario3_List[[i]]-1)
multi_scenario5_List <- lapply(1:5, FUN = function(i){temp1 <- temp_List[[i]]
1 + tempOther*temp1})
StepOne_5 <- automate_scenario_cfa(BetaBase = BetaAll,
                                    cov_mat_theo = cov_mat_theo,
                                    VarComp = VarComp,
                                    multi_List = multi_scenario5_List)



# Scenario 11: Random Directions with general increasing trend in X and decreasing in Y ------------
set.seed(1234)
d <- .32
multi_scenario6_List <- lapply(1:5, FUN = function(i) rbind(runif(n = 100, min = 1-d, max = 1+d),
                                                             runif(n = 100, min = 1-d, max = 1+d)))
multi_scenario6_List[[1]][,] <- 1
StepOne_6 <- automate_scenario_cfa(BetaBase = BetaAll,
                                    cov_mat_theo = cov_mat_theo,
                                    VarComp = VarComp,
                                    multi_List = multi_scenario6_List)
## Get Sigma matrices
ParameterValues_1 <- computeSigmaWithBurnin_linear(nT = 5,
                                                   BetaC_List = get_coefs(multi_scenario1_List, BetaAll),
                                                   CrossLagged_List = phi,    # if is not list, then all are equal
                                                   Sigma_C = cov_mat_theo,   # theoretical correlation matrix among Cs (including binary)
                                                   SigmaTarget = SigmaTarget, # Observed Covariance Among Xt,Yt
                                                   burnin = 49)

ParameterValues_2 <- computeSigmaWithBurnin_linear(nT = 5,
                                                   BetaC_List = get_coefs(multi_scenario2_List, BetaAll),
                                                   CrossLagged_List = phi,    # if is not list, then all are equal
                                                   Sigma_C = cov_mat_theo,   # theoretical correlation matrix among Cs (including binary)
                                                   SigmaTarget = SigmaTarget, # Observed Covariance Among Xt,Yt
                                                   burnin = 49)

ParameterValues_3 <- computeSigmaWithBurnin_linear(nT = 5,
                                                   BetaC_List = get_coefs(multi_scenario3_List, BetaAll),
                                                   CrossLagged_List = phi,    # if is not list, then all are equal
                                                   Sigma_C = cov_mat_theo,   # theoretical correlation matrix among Cs (including binary)
                                                   SigmaTarget = SigmaTarget, # Observed Covariance Among Xt,Yt
                                                   burnin = 49)

ParameterValues_4 <- computeSigmaWithBurnin_linear(nT = 5,
                                                   BetaC_List = get_coefs(multi_scenario4_List, BetaAll),
                                                   CrossLagged_List = phi,    # if is not list, then all are equal
                                                   Sigma_C = cov_mat_theo,   # theoretical correlation matrix among Cs (including binary)
                                                   SigmaTarget = SigmaTarget, # Observed Covariance Among Xt,Yt
                                                   burnin = 49)

ParameterValues_5 <- computeSigmaWithBurnin_linear(nT = 5,
                                                   BetaC_List = get_coefs(multi_scenario5_List, BetaAll),
                                                   CrossLagged_List = phi,    # if is not list, then all are equal
                                                   Sigma_C = cov_mat_theo,   # theoretical correlation matrix among Cs (including binary)
                                                   SigmaTarget = SigmaTarget, # Observed Covariance Among Xt,Yt
                                                   burnin = 49)

ParameterValues_6 <- computeSigmaWithBurnin_linear(nT = 5,
                                                   BetaC_List = get_coefs(multi_scenario6_List, BetaAll),
                                                   CrossLagged_List = phi,    # if is not list, then all are equal
                                                   Sigma_C = cov_mat_theo,   # theoretical correlation matrix among Cs (including binary)
                                                   SigmaTarget = SigmaTarget, # Observed Covariance Among Xt,Yt
                                                   burnin = 49)

sigmas <- list(ParameterValues_1,
                   ParameterValues_2,
                   ParameterValues_3,
                   ParameterValues_4,
                   ParameterValues_5,
                   ParameterValues_6)

saveRDS(sigmas, "params/highconf/sigmas.rds")

multi_mats <- list(
  multi_scenario1_List,
  multi_scenario2_List,
  multi_scenario3_List,
  multi_scenario4_List,
  multi_scenario5_List,
  multi_scenario6_List
)
saveRDS(multi_mats, "params/highconf/multi_mats.rds")

saveRDS(list(beta_cont = beta_cont_s1, beta_dich = beta_dich_s1), "params/highconf/betaAll.rds")

