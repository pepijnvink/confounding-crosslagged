#argv = commandArgs(trailingOnly = TRUE)
#outcome = argv[1]
##N=250
library(dplyr)
library(ggplot2)
library(jtools)
library(ggforce)
library(patchwork)
outcome = 'x'
model_colors <- c(
  "RI-CLPM"                          = "salmon",
  "RI-CLPM (Free Loadings)"          = "red",
  "DPM"                             = "blue",
  "DPM (Free Loadings)"                        = "skyblue",
  "CLPM"                            = "gold",
  "CL2PM"                      = "gold3"
)
allres <- vector("list", 4)
measures_power <- vector('list', 4)
plot_power <- vector('list', 4)
for(j in 1:4){
  N <- c(250, 500, 1000, 10000)[j]
  allres[[j]] <- vector("list", 4)
  for(i in 1:6){
  allres[[j]][[i]] <- readRDS(paste0("output/main/results/N", N, "/results_scenario", i, "_N", N, ".rds"))
  allres[[j]][[i]] <- bind_rows(allres[[j]][[i]]$output, .id = "model")
  }
  names(allres[[j]]) <- paste("Scenario", c(1, 2, "3b", "4b", "5d", 6))
  measures_all_df <- bind_rows(allres[[j]], .id = "Scenario") %>%
    mutate(Scenario = factor(Scenario, levels = paste("Scenario", c(1, 2, "3b", "4b", "5d", 6)), ordered = TRUE),
  N = paste("N =", N))
  measures_power[[j]] <- measures_all_df %>%
  dplyr::select(c(model, variable, time, true, estimate, power05_estimate, power05_mcse, Scenario)) %>%
  mutate(model = factor(model,
                        levels = c("clpm", "clpm_lag2", "riclpm", "riclpm_free", "dpm", "dpm_free"),
                        labels = c("CLPM", "CL2PM", "RI-CLPM", "RI-CLPM (Free Loadings)", "DPM", "DPM (Free Loadings)")),
         variable = factor(variable),
         Measure = "power",
         N = paste('N =', N),
         value = power05_estimate,
         mcse = power05_mcse) %>%
  dplyr::select(-c(power05_estimate, power05_mcse))
  gg <-  measures_power[[j]] %>%
  dplyr::mutate(yintercept = 0.05) %>%
  dplyr::filter(variable == outcome) %>%
  ggplot(aes(y = value,
             x = time,
             color = model)) +
  geom_point(size = 2) +
  geom_line(orientation = "x", linewidth = 0.75) +
  facet_grid(N~Scenario,
             scales = "free_y") +
  geom_hline(aes(yintercept = yintercept), linetype = "dashed") +
  scale_color_manual(values = model_colors) +
  geom_errorbar(aes(ymin = value-1.96*mcse, ymax = value+1.96*mcse), width = 0.25) +
  labs(x = "Timepoint",
       y = "Value")
  if(j == 1){
gg <- gg +
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
  } else if(j < 4){
      gg <- gg +
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
    } else {
      gg <- gg +
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
    }
plot_power[[j]] <- gg
}
fig <- wrap_plots(plot_power, ncol = 1, nrow = 4) + plot_layout(axis_titles = "collect")
ggsave(filename = paste0("figures/main/power/main_power_", outcome, ".png"), plot = fig, width = 9, height = 7, units = "in")
ggsave(filename = paste0("figures/main/power/main_power_", outcome, ".svg"), plot = fig, width = 9, height = 7, units = "in")
