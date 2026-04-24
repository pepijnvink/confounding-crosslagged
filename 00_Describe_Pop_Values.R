###############################################################################-
#### Supplement: Methods to Account for Unobserved Baseline Confounders with 
####             Time-Invariant and Time-Varying Effects in Cross-Lagged Panel 
####             Research
####
#### Script: Distribution of Correlation Coefficients from Inverse Wishart 
####         Simulation 
####
# load packages ----------------------------------------------------------------
library(LaplacesDemon)
library(tidyverse)

# setwd ------------------------------------------------------------------------
rstudioapi::documentPath() |> dirname() |> setwd()

# source functions -------------------------------------------------------------
source("00_helper_functions_PopAnalysis.R")
nu = 500
mean_denom <- nu-100-1
sigma_wis <- diag(100)
sigma_wis[lower.tri(sigma_wis)] <- 0
sigma_wis[upper.tri(sigma_wis)] <- 0
set.seed(42)
cor_mat <- rinvwishart(nu = nu, S = sigma_wis) %>% cov2cor()
cov_mat_theo <- transform_mixed_corr(Sigma = cor_mat, n_cont = 50, n_bin = 50)
cor_mat_theo <- cov2cor(cov_mat_theo)

Cors <- cor_mat[lower.tri(cor_mat)]
Cors_theo <- cor_mat_theo[lower.tri(cor_mat_theo)]
rbind(psych::describe(Cors), psych::describe(Cors_theo))

library(ggplot2)
p1 <- ggplot(data.frame(Cors), aes(x = Cors)) +
  geom_histogram(bins = 30, color = "black", fill = "lightblue") +
  labs(x = "Correlation Coefficients of Continuous and Latent Confounders",
       y = "Frequency") +
  theme_minimal() + coord_cartesian(xlim = c(-0.2, 0.2))
p2 <- ggplot(data.frame(Cors_theo), aes(x = Cors_theo)) +
  geom_histogram(bins = 30, color = "black", fill = "skyblue") +
  labs(x = "Correlation Coefficients of COntinuous and Binary Confounders",
       y = "Frequency") +
  theme_minimal() + coord_cartesian(xlim = c(-0.2, 0.2))
p <- ggpubr::ggarrange(p1, p2, ncol = 2, nrow = 1, labels = c("A", "B"))
ggsave("00_Figures/Inverse_Wishart_Correlation_Distribution.pdf", 
       plot = p, width = 10, height = 5)
