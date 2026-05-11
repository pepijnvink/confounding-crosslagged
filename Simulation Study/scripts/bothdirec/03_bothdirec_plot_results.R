argv = commandArgs(trailingOnly = TRUE)
N = as.numeric(argv[1])
outcome = argv[2]
##N=250
library(dplyr)
library(ggplot2)
library(jtools)
library(ggforce)
library(patchwork)
model_colors <- c(
  "RI-CLPM"                          = "salmon",
  "RI-CLPM (Free Loadings)"          = "red",
  "DPM"                             = "blue",
  "DPM (Free Loadings)"                        = "skyblue",
  "CLPM"                            = "gold",
  "CL2PM"                      = "gold3"
)
allres <- vector("list", 6)
for(i in 1:6){
  allres[[i]] <- readRDS(paste0("output/bothdirec/results/N", N, "/results_scenario", i, "_N", N, ".rds"))
  allres[[i]] <- bind_rows(allres[[i]]$output, .id = "model")
}
names(allres) <- paste("Scenario", c(1, 2, "3b", "4b", "5d", 6))
measures_all_df <- bind_rows(allres, .id = "Scenario") %>%
  mutate(Scenario = factor(Scenario, levels = paste("Scenario", c(1, 2, "3b", "4b", "5d", 6)), ordered = TRUE))

measures_rel_bias <- measures_all_df %>%
  dplyr::select(c(model, variable, time, true, estimate, rel_bias, rel_bias_mcse, Scenario)) %>%
  mutate(model = factor(model,
                        levels = c("clpm", "clpm_lag2", "riclpm", "riclpm_free", "dpm", "dpm_free"),
                        labels = c("CLPM", "CL2PM", "RI-CLPM", "RI-CLPM (Free Loadings)", "DPM", "DPM (Free Loadings)")),
         variable = factor(variable),
         Measure = "relative bias",
         value = rel_bias,
         mcse = rel_bias_mcse) %>%
  dplyr::select(-c(rel_bias, rel_bias_mcse))

measures_rmse <-  measures_all_df %>%
  dplyr::select(c(model, variable, time, true, estimate, rmse_estimate, rmse_mcse, Scenario)) %>%
  mutate(model = factor(model,
                        levels = c("clpm", "clpm_lag2", "riclpm", "riclpm_free", "dpm", "dpm_free"),
                        labels = c("CLPM", "CL2PM", "RI-CLPM", "RI-CLPM (Free Loadings)", "DPM", "DPM (Free Loadings)")),
         variable = factor(variable),
         Measure = "rmse",
         value = rmse_estimate,
         mcse = rmse_mcse) %>%
  dplyr::select(-c(rmse_estimate, rmse_mcse))

measures_coverage <-  measures_all_df %>%
  dplyr::select(c(model, variable, time, true, estimate, coverage_estimate, coverage_mcse, Scenario)) %>%
  mutate(model = factor(model,
                        levels = c("clpm", "clpm_lag2", "riclpm", "riclpm_free", "dpm", "dpm_free"),
                        labels = c("CLPM", "CL2PM", "RI-CLPM", "RI-CLPM (Free Loadings)", "DPM", "DPM (Free Loadings)")),
         variable = factor(variable),
         Measure = "coverage",
         value = coverage_estimate,
         mcse = coverage_mcse) %>%
  dplyr::select(-c(coverage_estimate, coverage_mcse))

measures_power <-  measures_all_df %>%
  dplyr::select(c(model, variable, time, true, estimate, power05_estimate, power05_mcse, Scenario)) %>%
  mutate(model = factor(model,
                        levels = c("clpm", "clpm_lag2", "riclpm", "riclpm_free", "dpm", "dpm_free"),
                        labels = c("CLPM", "CL2PM", "RI-CLPM", "RI-CLPM (Free Loadings)", "DPM", "DPM (Free Loadings)")),
         variable = factor(variable),
         Measure = "power",
         value = power05_estimate,
         mcse = power05_mcse) %>%
  dplyr::select(-c(power05_estimate, power05_mcse))

all_measures_long <- bind_rows(measures_rel_bias,
                               measures_rmse,
                               measures_coverage,
                               measures_power) %>%
  mutate(Scenario = factor(Scenario))


  ### plot for bias
plot_bias <-  all_measures_long %>%
    dplyr::mutate(Measure = factor(Measure,
                                   levels = c("relative bias", "rmse", "coverage", "power"),
                                   labels = c("Relative Bias", "RMSE", "Coverage Rate", "Power")),
                  yintercept = case_when(
                    Measure == "Relative Bias" ~ 0,
                    Measure == "RMSE" ~ 0.000001,
                    Measure == "Coverage Rate" ~ 0.95,
                    Measure == "Power" ~ 0.8
                  )) %>%
    dplyr::filter(variable == outcome, Measure == "Relative Bias") %>%
  ggplot(aes(y = value,
             x = time,
             color = model)) +
  geom_point(size = 2) +
  geom_line(orientation = "x", linewidth = 0.75) +
  facet_grid(Measure~Scenario,
             scales = "free_y") +
  coord_cartesian(ylim = c(-0.3,0.3)) +
geom_hline(aes(yintercept = yintercept), linetype = "dashed") +
  scale_color_manual(values = model_colors) +
  geom_errorbar(aes(ymin = value-1.96*mcse, ymax = value+1.96*mcse), width = 0.25) +
  labs(x = "Timepoint",
       y = "Value") +
  theme_minimal() +
  theme(axis.ticks = element_line(linewidth = 0.3, color = "grey40"),
      axis.ticks.length = unit(2, "pt"),
      strip.text = element_text(size = 9),
      panel.spacing = unit(0.6, "lines"),
      legend.text.align = 0,
      axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
    legend.position = "none")


### plot for RMSE
plot_rmse <-  all_measures_long %>%
  dplyr::mutate(Measure = factor(Measure,
                                 levels = c("relative bias", "rmse", "coverage", "power"),
                                 labels = c("Relative Bias", "RMSE", "Coverage Rate", "Power")),
                yintercept = case_when(
                  Measure == "Relative Bias" ~ 0,
                  Measure == "RMSE" ~ 0.000001,
                  Measure == "Coverage Rate" ~ 0.95,
                  Measure == "Power" ~ 0.8
                )) %>%
  dplyr::filter(variable == outcome, Measure == "RMSE") %>%
  ggplot(aes(y = value,
             x = time,
             color = model)) +
  geom_point(size = 2) +
  geom_line(orientation = "x", linewidth = 0.75) +
  facet_grid(Measure~Scenario,
             scales = "free_y") +
  coord_cartesian(ylim = c(0, 0.5)) +
  geom_hline(aes(yintercept = yintercept), linetype = "dashed") +
  scale_color_manual(values = model_colors) +
  geom_errorbar(aes(ymin = value-1.96*mcse, ymax = value+1.96*mcse), width = 0.25) +
  coord_cartesian(ylim = c(0,0.12)) +
  labs(x = "Timepoint",
       y = "Value") +
  theme_minimal() +
  theme(axis.ticks = element_line(linewidth = 0.3, color = "grey40"),
      axis.ticks.length = unit(2, "pt"),
      strip.text = element_text(size = 9),
      strip.background.x = element_blank(),
      strip.text.x = element_blank(),
      panel.spacing = unit(0.6, "lines"),
      legend.text.align = 0,
      axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
    legend.position = "none")

### plot for coverage
plot_coverage <-  all_measures_long %>%
  dplyr::mutate(Measure = factor(Measure,
                                 levels = c("relative bias", "rmse", "coverage", "power"),
                                 labels = c("Relative Bias", "RMSE", "Coverage Rate", "Power")),
                yintercept = case_when(
                  Measure == "Relative Bias" ~ 0,
                  Measure == "RMSE" ~ 0.000001,
                  Measure == "Coverage Rate" ~ 0.95,
                  Measure == "Power" ~ 0.8
                )) %>%
  dplyr::filter(variable == outcome, Measure == "Coverage Rate") %>%
  ggplot(aes(y = value,
             x = time,
             color = model)) +
  geom_point(size = 2) +
  geom_line(orientation = "x", linewidth = 0.75) +
  facet_grid(Measure~Scenario,
             scales = "free_y") +
  geom_hline(aes(yintercept = yintercept), linetype = "dashed") +
  scale_color_manual(values = model_colors) +
  geom_errorbar(aes(ymin = value-1.96*mcse, ymax = value+1.96*mcse), width = 0.25) +
  labs(x = "Timepoint",
       y = "Value") +
  theme_minimal() +
  theme(axis.ticks = element_line(linewidth = 0.3, color = "grey40"),
      axis.ticks.length = unit(2, "pt"),
      strip.text = element_text(size = 9),
      strip.background.x = element_blank(),
      strip.text.x = element_blank(),
      panel.spacing = unit(0.6, "lines"),
      legend.text.align = 0,
      axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none")

### plot for Power
plot_power <-  all_measures_long %>%
  dplyr::mutate(Measure = factor(Measure,
                                 levels = c("relative bias", "rmse", "coverage", "power"),
                                 labels = c("Relative Bias", "RMSE", "Coverage Rate", "Power")),
                yintercept = case_when(
                  Measure == "Relative Bias" ~ 0,
                  Measure == "RMSE" ~ 0.000001,
                  Measure == "Coverage Rate" ~ 0.95,
                  Measure == "Power" ~ 0.8
                )) %>%
  dplyr::filter(variable == outcome, Measure == "Power") %>%
  ggplot(aes(y = value,
             x = time,
             color = model)) +
  geom_point(size = 2) +
  geom_line(orientation = "x", linewidth = 0.75) +
  facet_grid(Measure~Scenario,
             scales = "free_y") +
  geom_hline(aes(yintercept = yintercept), linetype = "dashed") +
  scale_color_manual(values = model_colors) +
  geom_errorbar(aes(ymin = value-1.96*mcse, ymax = value+1.96*mcse), width = 0.25) +
  labs(x = "Timepoint",
       y = "Value") +
  theme_minimal() +
  theme(axis.ticks = element_line(linewidth = 0.3, color = "grey40"),
      axis.ticks.length = unit(2, "pt"),
      strip.text = element_text(size = 9),
      strip.background.x = element_blank(),
      strip.text.x = element_blank(),
      panel.spacing = unit(0.6, "lines"),
      legend.text.align = 0,
      legend.position = "bottom") +
  guides(color = guide_legend(nrow=2),
         shape = guide_legend(nrow = 2))



plot_all <- plot_bias/plot_rmse/plot_coverage/plot_power + plot_layout(axis_titles = "collect")
plot_all
ggsave(filename = paste0("figures/bothdirec/png/", outcome, "/N", N, "_", outcome, ".png"), plot = plot_all, width = 9, height = 7, units = "in")

ggsave(filename = paste0("figures/bothdirec/svg/", outcome, "/N", N, "_", outcome, ".svg"), plot = plot_all, width = 9, height = 7, units = "in")
