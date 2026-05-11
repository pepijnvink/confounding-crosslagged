###############################################################################-
#### Supplement: Methods to Account for Unobserved Baseline Confounders with 
####             Time-Invariant and Time-Varying Effects in Cross-Lagged Panel 
####             Research
####
#### Script: Population analysis: main article analysis
####         Model without Y->X path -> large effects of auto and cross

# load helper functions --------------------------------------------------------
rstudioapi::documentPath() |> dirname() |> setwd()
source("00_helper_functions_PopAnalysis.R")

# load packages ----------------------------------------------------------------
library(LaplacesDemon)
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

#############################################################################---
############### Scenarios ------------------------------------------------------
# Scenario 0: No Confounding ---------------------------------------------------
multi_scenario0_List <- lapply(1:5, FUN = function(i) matrix(0, nrow = 2, ncol = 100))
StepOne_0 <- automate_scenario_cfa(BetaBase = BetaAll, 
                                   cov_mat_theo = cov_mat_theo,
                                   VarComp = VarComp,
                                   multi_List = multi_scenario0_List)

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


# Scenario 3: Proportional Linear Increase -------------------------------------
# linear increase
multi_scenario3_List <- lapply(1:5, FUN = function(i) rbind(rep(1 + .025*(i-1), 100),
                                                            rep(1 - .1*(i-1), 100)))
StepOne_3 <- automate_scenario_cfa(BetaBase = BetaAll, 
                                   cov_mat_theo = cov_mat_theo, 
                                   VarComp = VarComp,
                                   multi_List = multi_scenario3_List)

# Scenario 4: Proportional Change ----------------------------------------------
# proportional but chaotic change

multi_scenario4_List <- lapply(1:5, FUN = function(i) rbind(rep(1 + .005*i*(-1.15)^(i-1)-.005, 100),
                                                            rep(1 + .005*i*(-1.25)^(i-1)-.005, 100)))
StepOne_4 <- automate_scenario_cfa(BetaBase = BetaAll, 
                                   cov_mat_theo = cov_mat_theo, 
                                   VarComp = VarComp,
                                   multi_List = multi_scenario4_List)


######## From here on, freely estimated DPM should no longer work --------------

###############----
# partially stable scenarios ---------------------------------------------------

# Scenario 5: Partially Stable Step Function -----------------------------------
# partial change in same direction at t = 3
multi_scenario5_List <- lapply(1:5, FUN = function(i) matrix(1, nrow = 2, 
                                                             ncol = 100))
multi_scenario5_3to5_blueprint <- cbind(mult.mat_s3_cont, mult.mat_s3_dich)
multi_scenario5_3to5 <- multi_scenario2_List[[5]]
multi_scenario5_3to5[multi_scenario5_3to5_blueprint == 1] <- 1
multi_scenario5_List[[3]] <- multi_scenario5_3to5
multi_scenario5_List[[4]] <- multi_scenario5_3to5
multi_scenario5_List[[5]] <- multi_scenario5_3to5
StepOne_5 <- automate_scenario_cfa(BetaBase = BetaAll, 
                                   cov_mat_theo = cov_mat_theo, 
                                   VarComp = VarComp,
                                   multi_List = multi_scenario5_List)

# Scenario 6: Partially Stable Linear Increase ---------------------------------
# same changing behavior as in scenario 5: some change
multi_scenario6_List <- multi_scenario3_List
multi_scenario6_List <- lapply(seq_along(multi_scenario6_List),
                               FUN = function(i){
                                 temp <- multi_scenario6_List[[i]]
                                 tempOther <- multi_scenario5_List[[5]]
                                 temp[tempOther == 1] <- 1
                                 temp})
StepOne_6 <- automate_scenario_cfa(BetaBase = BetaAll, 
                                   cov_mat_theo = cov_mat_theo, 
                                   VarComp = VarComp,
                                   multi_List = multi_scenario6_List)

# Scenario 7: Partially Stable Proportional Change -----------------------------
# same changing behavior as in scenario 5: some change
multi_scenario7_List <- multi_scenario4_List
multi_scenario7_List <- lapply(seq_along(multi_scenario7_List),
                               FUN = function(i){
                                 temp <- multi_scenario7_List[[i]]
                                 tempOther <- multi_scenario5_List[[5]]
                                 temp[tempOther == 1] <- 1
                                 temp})
StepOne_7 <- automate_scenario_cfa(BetaBase = BetaAll, 
                                   cov_mat_theo = cov_mat_theo, 
                                   VarComp = VarComp,
                                   multi_List = multi_scenario7_List)



# Scenario 8: Partially Stable Step Function -----------------------------------
# partial change at t = 3 and in opposite direction
tempOther <- cbind(mult.mat_s2_cont, mult.mat_s2_dich); tempOther[tempOther == 1] <- 0;
tempOther[tempOther > 1] <- 1; tempOther[tempOther < 1 & tempOther != 0] <- -1
temp_List <- lapply(1:5, FUN = function(i) multi_scenario2_List[[i]]-1)
multi_scenario8_List <- lapply(1:5, FUN = function(i){temp1 <- temp_List[[i]]
1 + tempOther*temp1})
StepOne_8 <- automate_scenario_cfa(BetaBase = BetaAll, 
                                   cov_mat_theo = cov_mat_theo, 
                                   VarComp = VarComp,
                                   multi_List = multi_scenario8_List)

# Scenario 9: Partially Stable Linear Increase ---------------------------------
# same changing behavior as in scenario 8: some change in opposite direction
tempOther <- multi_scenario8_List[[5]]; tempOther[tempOther == 1] <- 0;
tempOther[tempOther > 1] <- 1; tempOther[tempOther < 1 & tempOther != 0] <- -1
temp_List <- lapply(1:5, FUN = function(i) multi_scenario3_List[[i]]-1)
multi_scenario9_List <- lapply(1:5, FUN = function(i){temp1 <- temp_List[[i]]
1 + tempOther*temp1})

StepOne_9 <- automate_scenario_cfa(BetaBase = BetaAll, 
                                   cov_mat_theo = cov_mat_theo, 
                                   VarComp = VarComp,
                                   multi_List = multi_scenario9_List)

# Scenario 10: Partially Stable Proportional Change -----------------------------
# same changing behavior as in scenario 8: some change in opposite direction
tempOther <- multi_scenario8_List[[5]]; tempOther[tempOther == 1] <- 0;
tempOther[tempOther > 1] <- 1; tempOther[tempOther < 1 & tempOther != 0] <- -1
temp_List <- lapply(1:5, FUN = function(i) multi_scenario4_List[[i]]-1)
multi_scenario10_List <- lapply(1:5, FUN = function(i){temp1 <- temp_List[[i]]
1 + tempOther*temp1})
StepOne_10 <- automate_scenario_cfa(BetaBase = BetaAll, 
                                    cov_mat_theo = cov_mat_theo, 
                                    VarComp = VarComp,
                                    multi_List = multi_scenario10_List)



# Scenario 11: Random Directions with general increasing trend in X and decreasing in Y ------------
set.seed(1234)
d <- .32
multi_scenario11_List <- lapply(1:5, FUN = function(i) rbind(runif(n = 100, min = 1-d, max = 1+d),
                                                             runif(n = 100, min = 1-d, max = 1+d)))
multi_scenario11_List[[1]][,] <- 1
StepOne_11 <- automate_scenario_cfa(BetaBase = BetaAll, 
                                    cov_mat_theo = cov_mat_theo, 
                                    VarComp = VarComp,
                                    multi_List = multi_scenario11_List)



ListAll <- list(StepOne_0,
                StepOne_1, 
                StepOne_2, 
                StepOne_3, 
                StepOne_4, 
                StepOne_5, 
                StepOne_6, 
                StepOne_7,
                StepOne_8,
                StepOne_9,
                StepOne_10,
                StepOne_11)
FitMeasures <- t(sapply(ListAll, "[[", "fit_measures"))
FitMeasures

###############################################################################-
##### generate covariance matrices --------------------------------------------
# Scenario 0
ParameterValues_0 <- computeSigmaWithBurnin_linear(
  nT = 5,
  BetaC_List = StepOne_0$B_XY_List,
  CrossLagged_List = phi,
  Sigma_C = cov_mat_theo,
  SigmaTarget = SigmaTarget,
  burnin = 49
)
ParameterValues_0$Scenario <- "Scenario 0"
ParameterValues_0$ScenarioSubtitle <- "No Confounding"

# Scenario 1
ParameterValues_1 <- computeSigmaWithBurnin_linear(
  nT = 5,
  BetaC_List = StepOne_1$B_XY_List,
  CrossLagged_List = phi,
  Sigma_C = cov_mat_theo,
  SigmaTarget = SigmaTarget,
  burnin = 49
)
ParameterValues_1$Scenario <- "Scenario 1"
ParameterValues_1$ScenarioSubtitle <- "Stable Confounder Effects"

# Scenario 2
ParameterValues_2 <- computeSigmaWithBurnin_linear(
  nT = 5,
  BetaC_List = StepOne_2$B_XY_List,
  CrossLagged_List = phi,
  Sigma_C = cov_mat_theo,
  SigmaTarget = SigmaTarget,
  burnin = 49
)
ParameterValues_2$Scenario <- "Scenario 2"
ParameterValues_2$ScenarioSubtitle <- "Step"

# Scenario 3
ParameterValues_3 <- computeSigmaWithBurnin_linear(
  nT = 5,
  BetaC_List = StepOne_3$B_XY_List,
  CrossLagged_List = phi,
  Sigma_C = cov_mat_theo,
  SigmaTarget = SigmaTarget,
  burnin = 49
)
ParameterValues_3$Scenario <- "Scenario 3"
ParameterValues_3$ScenarioSubtitle <- "Linear"

# Scenario 4
ParameterValues_4 <- computeSigmaWithBurnin_linear(
  nT = 5,
  BetaC_List = StepOne_4$B_XY_List,
  CrossLagged_List = phi,
  Sigma_C = cov_mat_theo,
  SigmaTarget = SigmaTarget,
  burnin = 49
)
ParameterValues_4$Scenario <- "Scenario 4"
ParameterValues_4$ScenarioSubtitle <- "Proportional"

# Scenario 5
ParameterValues_5 <- computeSigmaWithBurnin_linear(
  nT = 5,
  BetaC_List = StepOne_5$B_XY_List,
  CrossLagged_List = phi,
  Sigma_C = cov_mat_theo,
  SigmaTarget = SigmaTarget,
  burnin = 49
)
ParameterValues_5$Scenario <- "Scenario 5"
ParameterValues_5$ScenarioSubtitle <- "Part-Step"

# Scenario 6
ParameterValues_6 <- computeSigmaWithBurnin_linear(
  nT = 5,
  BetaC_List = StepOne_6$B_XY_List,
  CrossLagged_List = phi,
  Sigma_C = cov_mat_theo,
  SigmaTarget = SigmaTarget,
  burnin = 49
)
ParameterValues_6$Scenario <- "Scenario 6"
ParameterValues_6$ScenarioSubtitle <- "Part-Linear"

# Scenario 7
ParameterValues_7 <- computeSigmaWithBurnin_linear(
  nT = 5,
  BetaC_List = StepOne_7$B_XY_List,
  CrossLagged_List = phi,
  Sigma_C = cov_mat_theo,
  SigmaTarget = SigmaTarget,
  burnin = 49
)
ParameterValues_7$Scenario <- "Scenario 7"
ParameterValues_7$ScenarioSubtitle <- "Part-Proportional"

# Scenario 8
ParameterValues_8 <- computeSigmaWithBurnin_linear(
  nT = 5,
  BetaC_List = StepOne_8$B_XY_List,
  CrossLagged_List = phi,
  Sigma_C = cov_mat_theo,
  SigmaTarget = SigmaTarget,
  burnin = 49
)
ParameterValues_8$Scenario <- "Scenario 8"
ParameterValues_8$ScenarioSubtitle <- "Part-Step Opposite"

# Scenario 9
ParameterValues_9 <- computeSigmaWithBurnin_linear(
  nT = 5,
  BetaC_List = StepOne_9$B_XY_List,
  CrossLagged_List = phi,
  Sigma_C = cov_mat_theo,
  SigmaTarget = SigmaTarget,
  burnin = 49
)
ParameterValues_9$Scenario <- "Scenario 9"
ParameterValues_9$ScenarioSubtitle <- "Part-Linear Opposite"

# Scenario 10
ParameterValues_10 <- computeSigmaWithBurnin_linear(
  nT = 5,
  BetaC_List = StepOne_10$B_XY_List,
  CrossLagged_List = phi,
  Sigma_C = cov_mat_theo,
  SigmaTarget = SigmaTarget,
  burnin = 49
)
ParameterValues_10$Scenario <- "Scenario 10"
ParameterValues_10$ScenarioSubtitle <- "Part-Proportional Opposite"

# Scenario 11
ParameterValues_11 <- computeSigmaWithBurnin_linear(
  nT = 5,
  BetaC_List = StepOne_11$B_XY_List,
  CrossLagged_List = phi,
  Sigma_C = cov_mat_theo,
  SigmaTarget = SigmaTarget,
  burnin = 49
)
ParameterValues_11$Scenario <- "Scenario 11"
ParameterValues_11$ScenarioSubtitle <- "Chaotic Effects"

ParList <- list(ParameterValues_0,
                ParameterValues_1,
                ParameterValues_2,
                ParameterValues_3,
                ParameterValues_4,
                ParameterValues_5,
                ParameterValues_6,
                ParameterValues_7,
                ParameterValues_8,
                ParameterValues_9,
                ParameterValues_10,
                ParameterValues_11)
R2_List <- lapply(ParList, 
                  FUN =  function(PL){
                    Res_R2_C <- data.frame(t(sapply(PL$VarExpC_List, diag)))
                    Res_R2 <- data.frame(t(1-sapply(PL$Sigma_Eps_List, diag)))
                    Res_R2_incr_XY <- Res_R2-Res_R2_C
                    names(Res_R2_C) <- c("R2_C_X", "R2_C_Y")
                    names(Res_R2) <- c("R2_X", "R2_Y")
                    names(Res_R2_incr_XY) <- c("R2_incr_X", "R2_incr_Y")
                    Res_out <- cbind(Res_R2_C, Res_R2, Res_R2_incr_XY)
                    Res_out$Scenario <- PL$Scenario
                    Res_out$ScenarioSubtitle <- PL$ScenarioSubtitle
                    return(Res_out)})
dataR2 <- do.call(rbind, R2_List)
dataR2 <- reshape2::melt(data = dataR2, 
                         measure.vars = c("R2_C_X", "R2_C_Y",
                                          "R2_X", "R2_Y",
                                          "R2_incr_X", "R2_incr_Y"), 
                         value.name = "R2", variable.name = "R2_name")
dataR2$XY <- "X"; dataR2$XY[grepl(pattern = "Y", dataR2$R2_name)] <- "Y"
dataR2$R2_type <- dataR2$R2_name
dataR2$R2_type <- stringr::str_remove_all(string = dataR2$R2_type, pattern = "_X")
dataR2$R2_type <- stringr::str_remove_all(string = dataR2$R2_type, pattern = "_Y")
dataR2$time <- 0:5; dataR2 <- dataR2[dataR2$time>0, ]


# eigenvalues
SmallestEigenValues_List <- lapply(ParList, 
                                   FUN =  function(PL){
                                     Res <- data.frame(min(eigen(PL$Sigma_full)$values))
                                     names(Res) <- c("EValue")
                                     Res$Scenario <- PL$Scenario; Res$ScenarioSubtitle <- PL$ScenarioSubtitle
                                     Res[, c(2,3, 1)]})
SmallestEigenValues <- do.call(rbind, SmallestEigenValues_List)
SmallestEigenValues # looks fine

###############################################################################-
#### fit models -----
ListFit <- function(L)
{
  out <- FitModels(Sigma_full = L$Sigma_full)
  out$Scenario <- L$Scenario
  out$ScenarioSubtitle <- L$ScenarioSubtitle
  return(out)
}

# parallel setup
library(future)
library(future.apply)

# Use multisession (works on Windows/macOS/Linux). Spins up 6 workers.
plan(multisession, workers = 6)
on.exit(plan(sequential), add = TRUE)  # restore after we're done

# If ListFit / FitModels need packages on workers, list them here
pkgs_needed <- c("lavaan")  # add others if used inside (e.g., "Matrix", "dplyr")

# Parallel apply
res_list <- future_lapply(
  ParList,
  FUN = ListFit,
  future.seed = TRUE,                 # reproducible RNG per worker
  future.packages = pkgs_needed       # load on each worker
)

# Combine results (base R)
Results <- do.call(rbind, res_list)
names(Results)

save(list = c("Results", "SmallestEigenValues", "dataR2", 
              "ParList", "FitMeasures", "ListAll", 
              "StepOne_0", 
              "StepOne_1",
              "StepOne_2",
              "StepOne_3",
              "StepOne_4",
              "StepOne_5",
              "StepOne_6",
              "StepOne_7",
              "StepOne_8",
              "StepOne_9",
              "StepOne_10",
              "StepOne_11"),
     file = "00_Data/01_PopulationAnalysis_NoYX_Large_C.rda")
