##############################################################################--
#### Supplement: Methods to Account for Unobserved Baseline Confounders with 
####             Time-Invariant and Time-Varying Effects in Cross-Lagged Panel 
####             Research
#### 
#### Script: helper functions for plotting population level analyses
####

###############################################################################-
# Packages ---------------------------------------------------------------------
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(scales)
library(reshape2)

# RelBias-correction (remove infinite)
Results$RelBiasPercent[abs(Results$RelBiasPercent) > 10^5] <- NA

# remove auto effects ----------------------------------------------------------
if(removeAuto)
{
  Results <- subset(Results, effect != "X[t-1] -> X[t]" & 
                      effect != "Y[t-1] -> Y[t]")
}

# remove zero cross-lagged effect ----------------------------------------------
if(removeZeroCross)
{
  Results <- subset(Results, effect != "Y[t-1] -> X[t]")
}

# adapt vertical height of figures ---------------------------------------------
vertical_scaling <- sqrt(1 - removeAuto * 0.5 - removeZeroCross * 0.25)

# --- for two kinds of plot-----------------------------------------------------
# Scenario relabeling + ordering (used across all plots)
scenario_map <- tibble::tibble(
  Scenario_raw = paste0("Scenario ", 0:11),
  Scenario_new = c("Scenario 0", "Scenario 1", "Scenario 2", "Scenario 3A",
                   "Scenario 3B", "Scenario 4A", "Scenario 5A", "Scenario 5C",
                   "Scenario 4B", "Scenario 5B", "Scenario 5D", "Scenario 6")
)
 # Explicit sorting for plotting: alphabetical within each number
scenario_levels <- c(
    "Scenario 0", "Scenario 1", "Scenario 2",
    "Scenario 3A", "Scenario 3B",
    "Scenario 4A", "Scenario 4B",
    "Scenario 5A", "Scenario 5B", "Scenario 5C", "Scenario 5D",
    "Scenario 6")

apply_scenario_map <- function(df) {
  df |>
    dplyr::left_join(scenario_map, by = c("Scenario" = "Scenario_raw")) |>
    dplyr::mutate(
      Scenario = factor(Scenario_new, levels = scenario_levels)
    ) |>
    dplyr::select(-Scenario_new)
}

group_A <- c("Scenario 0", "Scenario 1", "Scenario 2", "Scenario 3A", "Scenario 3B")
group_B <- setdiff(scenario_levels, group_A)
# empty plot for better overview
pEmpty <- ggplot() + theme_void() + coord_fixed(ratio = 1)

###############################################################################-
# Visualization of Parameters of C across time ---------------------------------

# ---- Function to create one plot given BX, BY and titles ----
make_plot_and_data <- function(StepOne, scenario, subtitle, 
                               return_data = FALSE) {
  BX <- data.frame(StepOne$BX); BY <- data.frame(StepOne$BY)
  BX$time <- 1:5; BY$time <- 1:5
  
  temp1 <- melt(BX, id.vars = "time", value.name = "Effect")
  temp2 <- melt(BY, id.vars = "time", value.name = "Effect")
  
  temp1$DV <- "X"; temp2$DV <- "Y"
  df <- rbind(temp1, temp2)
  df$Scenario <- scenario
  df$ScenarioSubtitle <- subtitle
  
  if(return_data) return(df)
  
  p <- ggplot(df, aes(x = time, y = Effect, col = variable)) +
    geom_line() +
    facet_grid(DV ~ .) +
    theme_minimal() +
    theme(legend.position = "none") +
    ggtitle(scenario, subtitle)
  return(p)
}


# ---- Extract Data ----
data0  <- make_plot_and_data(StepOne_0,  ParList[[1]]$Scenario,  
                             ParList[[1]]$ScenarioSubtitle,  
                             return_data = TRUE)

data1  <- make_plot_and_data(StepOne_1,  ParList[[2]]$Scenario,  
                             ParList[[2]]$ScenarioSubtitle,  
                             return_data = TRUE)

data2  <- make_plot_and_data(StepOne_2,  ParList[[3]]$Scenario,  
                             ParList[[3]]$ScenarioSubtitle,  
                             return_data = TRUE)

data3  <- make_plot_and_data(StepOne_3,  ParList[[4]]$Scenario,  
                             ParList[[4]]$ScenarioSubtitle,  
                             return_data = TRUE)

data4  <- make_plot_and_data(StepOne_4,  ParList[[5]]$Scenario,  
                             ParList[[5]]$ScenarioSubtitle,  
                             return_data = TRUE)

data5  <- make_plot_and_data(StepOne_5,  ParList[[6]]$Scenario,  
                             ParList[[6]]$ScenarioSubtitle,  
                             return_data = TRUE)

data6  <- make_plot_and_data(StepOne_6,  ParList[[7]]$Scenario,  
                             ParList[[7]]$ScenarioSubtitle,  
                             return_data = TRUE)

data7  <- make_plot_and_data(StepOne_7,  ParList[[8]]$Scenario,  
                             ParList[[8]]$ScenarioSubtitle,  
                             return_data = TRUE)

data8  <- make_plot_and_data(StepOne_8,  ParList[[9]]$Scenario,  
                             ParList[[9]]$ScenarioSubtitle,  
                             return_data = TRUE)

data9  <- make_plot_and_data(StepOne_9,  ParList[[10]]$Scenario,  
                             ParList[[10]]$ScenarioSubtitle,  
                             return_data = TRUE)

data10 <- make_plot_and_data(StepOne_10, ParList[[11]]$Scenario, 
                             ParList[[11]]$ScenarioSubtitle, 
                             return_data = TRUE)

data11 <- make_plot_and_data(StepOne_11, ParList[[12]]$Scenario, 
                             ParList[[12]]$ScenarioSubtitle, 
                             return_data = TRUE)

data <- rbind(data0, data1, data2, data3, data4, data5, data6, data7, data8, data9,
              data10, data11)

# Apply scenario labels/order consistently
# For plots based on `data`
data <- apply_scenario_map(data)

# Helper vectors to define which scenarios are selected / not selected
all_scenarios          <- scenario_levels
selected_scenarios     <- c("Scenario 1", "Scenario 2", "Scenario 3B", "Scenario 4B", "Scenario 5D", "Scenario 6")
not_selected_scenarios <- setdiff(all_scenarios, selected_scenarios)


# Width factors for selected / not selected plots (relative to "all")
width_factor_selected     <- length(selected_scenarios) / length(all_scenarios)
width_factor_not_selected <- length(not_selected_scenarios) / length(all_scenarios)

dataOK   <- data[(data$Scenario %in% group_A), ]
dataBias <- data[!(data$Scenario %in% group_A), ]
data$Scenario     <- factor(data$Scenario,     levels = scenario_levels)
dataOK$Scenario   <- factor(dataOK$Scenario,   levels = scenario_levels)
dataBias$Scenario <- factor(dataBias$Scenario, levels = scenario_levels)


# ---- Arrange them all in one figure (ALL) ----
p_Scenarios_all <- ggplot(data, aes(x = time, y = Effect, col = variable)) +
  geom_line() +
  facet_grid(DV ~ Scenario) +
  theme_minimal() +
  theme(legend.position = "none") 

ggsave(paste0("00_Figures/", Prefix, "ScenariosAll.pdf"), 
       plot = p_Scenarios_all, device = pdf, 
       width = 22, height = 5)

# ---- Split A/B effects of C (ALL) ----
pA_Betas <- ggplot(dataOK, aes(x = time, y = Effect, col = variable)) +
  geom_line() +
  facet_grid(DV ~ Scenario) +
  theme_minimal() +
  theme(legend.position = "none") 

pB_Betas <- ggplot(dataBias, aes(x = time, y = Effect, col = variable)) +
  geom_line() +
  facet_grid(DV ~ Scenario) +
  theme_minimal() +
  theme(legend.position = "none") 

pA_Betas_combined <- ggarrange(pA_Betas, pEmpty, ncol = 2, widths = c(5, 1.8),
                               common.legend = TRUE, legend = "none")
pA_Betas_double <- ggarrange(pA_Betas_combined, pB_Betas, ncol = 1, labels = "AUTO")
ggsave(paste0("00_Figures/", Prefix, "Betas.pdf"), pA_Betas_double, 
       width = 16, height = 8, device = pdf)

###############################################################################-
## Fancy Plot with variance decomposition (split A/B) -------------------------
dataR2$R2_type <- factor(dataR2$R2_type,
                         levels = c("R2_C", "R2", "R2_incr"))
dataR2$R2_type_lab <- dplyr::recode(
  dataR2$R2_type,
  "R2_C"   = "R[C]^2",
  "R2"     = "R[full]^2",
  "R2_incr"= "R[increment]^2"
)
dataR2$R2_type_lab <- factor(
  dataR2$R2_type_lab,
  levels = c("R[C]^2", "R[full]^2", "R[increment]^2")
)

# Apply scenario relabeling to dataR2

dataR2 <- dataR2 %>%
  mutate(Scenario = factor(Scenario, levels = unique(dataR2$Scenario)))

# Apply scenario labels/order consistently

dataR2 <- apply_scenario_map(dataR2)

wide <- dataR2 %>%
  filter(R2_type %in% c("R2_C", "R2")) %>%
  mutate(R2_type = recode(R2_type, "R2" = "R2_full")) %>%
  select(Scenario, ScenarioSubtitle, time, XY, R2_type, R2) %>%
  tidyr::pivot_wider(names_from = R2_type, values_from = R2) %>%
  mutate(R2_incr = R2_full - R2_C)

lines_long <- wide %>%
  select(Scenario, ScenarioSubtitle, time, XY, R2_full, R2_C) %>%
  pivot_longer(c(R2_full, R2_C), names_to = "what", values_to = "R2") %>%
  mutate(what = recode(what,
                       "R2_full" = "R[full]^2",
                       "R2_C"    = "R[C]^2"))

make_R2_plot <- function(data_wide, data_lines) {
  ggplot() +
    geom_ribbon(
      data = data_wide,
      aes(x = time, ymin = R2_full, ymax = 1, fill = "Unexplained"),
      alpha = 0.25
    ) +
    geom_ribbon(
      data = data_wide,
      aes(x = time, ymin = 0, ymax = R2_C, fill = "R[C]^2"),
      alpha = 0.6
    ) +
    geom_ribbon(
      data = data_wide,
      aes(x = time, ymin = R2_C, ymax = R2_full, fill = "R[incr]^2"),
      alpha = 0.8
    ) +
    geom_line(
      data = data_lines,
      aes(x = time, y = R2, color = what),
      linewidth = 0.6
    ) +
    geom_hline(yintercept = c(0, 1), color = "grey70", linewidth = 0.3) +
    facet_grid(
      XY ~ Scenario,
      labeller = labeller(.rows = label_value, .cols = label_value)
    ) +
    scale_color_manual(
      values = c("R[C]^2" = "blue", "R[full]^2" = "blue4"),
      labels = label_parse()
    ) +
    scale_fill_manual(
      values = c(
        "R[C]^2"      = "skyblue1",
        "R[incr]^2"   = "deepskyblue2",
        "Unexplained" = "grey90"
      ),
      labels = c(
        expression(R[C]^2),
        expression(R[incr]^2),
        "Unexplained"
      )
    ) +
    scale_x_continuous(breaks = 1:5, expand = c(0, 0)) +
    scale_y_continuous(
      breaks = seq(0, 1, 0.1),
      limits = c(0, 1),
      expand = c(0, 0),
      labels = scales::number_format(accuracy = 0.1)
    ) +
    labs(
      x = "time",
      y = expression(R^2),
      color = " ",
      fill  = " "
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_line(linewidth = 0.3, color = "grey40"),
      axis.ticks.length = unit(2, "pt"),
      strip.text = element_text(size = 9),
      panel.spacing = unit(0.6, "lines"),
      legend.text.align = 0,
      legend.position = "bottom"
    )
}

# Plain labels for splitting into Scenario A (0–4) and B (5–11)
scenarios_A_plain <- group_A
scenarios_B_plain <- setdiff(scenario_levels, group_A)

# --- Scenario A (0–4), ALL ---
wide_A  <- subset(wide,  Scenario %in% scenarios_A_plain)
lines_A <- subset(lines_long, Scenario %in% scenarios_A_plain)
pA_R2   <- make_R2_plot(wide_A, lines_A)

# --- Scenario B (5–11), ALL ---
wide_B  <- subset(wide,  Scenario %in% scenarios_B_plain)
lines_B <- subset(lines_long, Scenario %in% scenarios_B_plain)
pB_R2   <- make_R2_plot(wide_B, lines_B)

pA_R2_combined <- ggarrange(pA_R2, pEmpty, ncol = 2, widths = c(5, 1.8),
                            common.legend = TRUE, legend = "none")
pA_R2_double <- ggarrange(pA_R2_combined, pB_R2, ncol = 1,
                          labels = "AUTO", heights = c(1, 1.35))

ggsave(
  paste0("00_Figures/", Prefix, "ScenariosR2All_incr_double.pdf"),
  plot = pA_R2_double, device = pdf,
  width = 16, height = 8
)

write.csv(data.frame(rbind(wide_A, wide_B)), paste0("00_Data/PlotData/", Prefix, 
                       "ScenariosR2All_incr_double_data.csv"), 
          row.names = FALSE)

###############################################################################-
### Prepare Results object
# 1) Pretty facet labels using plotmath
Results <- Results %>%
  mutate(
    effect_parsed = factor(
      effect,
      levels = c("X[t-1] -> X[t]", "Y[t-1] -> Y[t]", "Y[t-1] -> X[t]", "X[t-1] -> Y[t]"),
      labels = c(
        "X[t-1] %->% X[t]",
        "Y[t-1] %->% Y[t]",
        "Y[t-1] %->% X[t]",
        "X[t-1] %->% Y[t]"
      )
    )
  )

# 2) Rename model values to human-readable labels and align with colors
Results <- Results %>%
  mutate(
    model_lab = recode(
      model,
      "model_covariates"        = "True Model",
      "model_clpm"              = "CLPM",
      "model_clpm_lag2"         = "CLPM lag 2",
      "model_dpm"               = "DPM",
      "model_dpm_free"          = "DPM-free",
      "model_riclpm"            = "RICLPM",
      "model_riclpm_free"       = "RICLPM-free"
    )
  )

# 3) Extend your color map to include the two CLPM+Covariates variants
model_colors <- c(
  "True Model"                      = "black",
  "RICLPM"                          = "salmon",
  "RICLPM-free"                     = "red",
  "DPM"                             = "blue",
  "DPM-free"                        = "skyblue",
  "CLPM"                            = "gold",
  "CLPM lag 2"                      = "gold3"
)

# Optional: fix model order in legend
Results$model_lab <- factor(
  Results$model_lab,
  levels = c(
    "True Model", 
    "CLPM", "CLPM lag 2",
    "DPM", "DPM-free",
    "RICLPM", "RICLPM-free"
  )
)

# Keep original scenario labels (with spaces) for selection logic
Results <- apply_scenario_map(Results)
Results$Scenario_plain <- Results$Scenario
Results$Scenario <- stringr::str_replace_all(string = Results$Scenario, 
                                             pattern = " ", replacement = "~")
Results$ScenarioSubtitle <- stringr::str_replace_all(string = Results$ScenarioSubtitle, 
                                                     pattern = " ", replacement = "~")
Results$Scenario <- gsub("([0-9])([A-Za-z])", "\\1*italic('\\2')", Results$Scenario)
scenario_levels_parsed <- stringr::str_replace_all(scenario_levels, " ", "~")
scenario_levels_parsed <- gsub("([0-9])([A-Za-z])", "\\1*italic('\\2')", scenario_levels_parsed)
Results$Scenario <- factor(Results$Scenario, levels = scenario_levels_parsed)
group_A <- scenario_levels_parsed[1:5]
group_B <- setdiff(scenario_levels_parsed, group_A)

# Relative bias plotting data (drop NA panels)
Results_rel <- Results |>
  dplyr::filter(!is.na(RelBiasPercent)) |>
  dplyr::mutate(
    Scenario = forcats::fct_drop(Scenario),
    effect_parsed = forcats::fct_drop(effect_parsed)
  )

write.csv(Results_rel, paste0("00_Data/PlotData/", Prefix, "Results_rel.csv"), 
          row.names = FALSE)

###############################################################################-
#### Plot Bias, RelBias, Expected Z and Average squared Bias -------------------

# ============================================================================ #
# 1) Bias -------------------------------------------------------------------- #
# ============================================================================ #
pA_Bias <- ggplot(subset(Results, Scenario %in% group_A),
                  aes(x = time, y = Bias, col = model_lab, group = model_lab)) +
  geom_line() + geom_point() +
  facet_grid(effect_parsed ~ Scenario, 
             labeller = label_parsed, scales = "free_y") +
  scale_color_manual(values = model_colors, breaks = names(model_colors)) +
  labs(color = "Models", x = "Time t", y = "Bias in Coefficient") +
  theme_minimal() + theme(legend.position = "none")

pB_Bias <- ggplot(subset(Results, Scenario %in% group_B),
                  aes(x = time, y = Bias, col = model_lab, group = model_lab)) +
  geom_line() + geom_point() +
  facet_grid(effect_parsed ~ Scenario, 
             labeller = label_parsed, scales = "free_y") +
  scale_color_manual(values = model_colors, breaks = names(model_colors)) +
  labs(color = "Models", x = "Time t", y = "Bias in Coefficient") +
  theme_minimal() + theme(legend.position = "bottom")

pA_Bias_combined <- ggarrange(pA_Bias, pEmpty, ncol = 2, widths = c(5, 1.8),
                              common.legend = TRUE, legend = "none")

p_Bias_double <- ggarrange(pA_Bias_combined, pB_Bias, ncol = 1, labels = "AUTO", 
                           heights = c(1, 1.35))
ggsave(paste0("00_Figures/", Prefix, "Bias.pdf"), p_Bias_double, 
       width = 16, height = vertical_scaling*16, device = pdf)

# ============================================================================ #
# 2) Expected Z-Value -------------------------------------------------------- #
# ============================================================================ #
pA_Z <- ggplot(subset(Results, Scenario %in% group_A),
               aes(x = time, y = z, col = model_lab, group = model_lab)) +
  geom_line() + geom_point() +
  facet_grid(effect_parsed ~ Scenario,
             labeller = label_parsed, scales = "free_y") +
  scale_color_manual(values = model_colors, breaks = names(model_colors)) +
  labs(color = "Models", x = "Time t", y = "Expected Z-Value for N = 1000") +
  theme_minimal() + theme(legend.position = "none")

pB_Z <- ggplot(subset(Results, Scenario %in% group_B),
               aes(x = time, y = z, col = model_lab, group = model_lab)) +
  geom_line() + geom_point() +
  facet_grid(effect_parsed ~ Scenario,
             labeller = label_parsed, scales = "free_y") +
  scale_color_manual(values = model_colors, breaks = names(model_colors)) +
  labs(color = "Models", x = "Time t", y = "Expected Z-Value for N = 1000") +
  theme_minimal() + theme(legend.position = "bottom")

pA_Z_combined <- ggarrange(pA_Z, pEmpty, ncol = 2, widths = c(5, 1.8),
                           common.legend = TRUE, legend = "none")
p_Z_double <- ggarrange(pA_Z_combined, pB_Z, ncol = 1, labels = "AUTO",
                        heights = c(1, 1.35))
ggsave(paste0("00_Figures/", Prefix, "Z.pdf"), p_Z_double, 
       width = 16, height = vertical_scaling*16, device = pdf)

# ============================================================================ #
# 3) Relative Bias ----------------------------------------------------------- #
# ============================================================================ #
pA_RelBias <- ggplot(subset(Results_rel, Scenario %in% group_A),
                     aes(x = time, y = RelBiasPercent, col = model_lab, group = model_lab)) +
  geom_line() + geom_point() +
  facet_grid(effect_parsed ~ Scenario, labeller = label_parsed) +
  scale_color_manual(values = model_colors, breaks = names(model_colors)) +
  labs(color = "Models", x = "Time t",
       y = "Relative Bias\nin Coefficient [in %]") +
  theme_minimal() + theme(legend.position = "none")

pB_RelBias <- ggplot(subset(Results_rel, Scenario %in% group_B),
                     aes(x = time, y = RelBiasPercent, col = model_lab, group = model_lab)) +
  geom_line() + geom_point() +
  facet_grid(effect_parsed ~ Scenario, labeller = label_parsed) +
  scale_color_manual(values = model_colors, breaks = names(model_colors)) +
  labs(color = "Models", x = "Time t",
       y = "Relative Bias\nin Coefficient [in %]") +
  theme_minimal() + theme(legend.position = "bottom")

pA_RelBias_combined <- ggarrange(pA_RelBias, pEmpty, ncol = 2, widths = c(5, 1.8),
                                 common.legend = TRUE, legend = "none")
p_RelBias_double <- ggarrange(pA_RelBias_combined, pB_RelBias, ncol = 1, labels = "AUTO",
                              heights = c(1, 1.35))
ggsave(paste0("00_Figures/", Prefix, "RelBias.pdf"), p_RelBias_double, 
       width = 16, height = vertical_scaling*16, device = pdf)

# Zoom variants (only y-limits change, rest identical) for ALL
p_RelBias_double_100 <- ggarrange(
  ggarrange(
    pA_RelBias + coord_cartesian(ylim = c(-100, 100)),
    pEmpty, ncol = 2, 
    widths = c(5, 1.8),
    common.legend = TRUE, legend = "none"
  ),
  pB_RelBias + coord_cartesian(ylim = c(-100, 100)),
  ncol = 1, labels = "AUTO",
  heights = c(1, 1.35)
)
ggsave(paste0("00_Figures/", Prefix, "RelBiasZoomed100.pdf"), p_RelBias_double_100, 
       width = 16, height = vertical_scaling*16, device = pdf)

p_RelBias_double_30 <- ggarrange(
  ggarrange(
    pA_RelBias + coord_cartesian(ylim = c(-30, 30)),
    pEmpty, ncol = 2, 
    widths = c(5, 1.8),
    common.legend = TRUE, legend = "none"
  ),
  pB_RelBias + coord_cartesian(ylim = c(-30, 30)),
  ncol = 1, labels = "AUTO",
  heights = c(1, 1.35)
)
ggsave(paste0("00_Figures/", Prefix, "RelBiasZoomed30.pdf"), p_RelBias_double_30, 
       width = 16, height = vertical_scaling*16, device = pdf)

# ============================================================================ #
# 4) Average Squared Biased across coefficients ------------------------------ #
# ============================================================================ #
RMSE_df <- Results %>%
  group_by(model_lab, time, Scenario, ScenarioSubtitle, Scenario_plain) %>%
  summarise(RMSE = sqrt(mean(Bias^2)), .groups = "drop")

write.csv(RMSE_df, paste0("00_Data/PlotData/", Prefix, "RMSE_df.csv"), 
          row.names = FALSE)

pA_RMSE <- ggplot(subset(RMSE_df, Scenario %in% group_A),
                  aes(x = time, y = RMSE, col = model_lab)) +
  geom_point() + geom_line() +
  facet_grid(. ~ Scenario, labeller = label_parsed) +
  scale_color_manual(values = model_colors, breaks = names(model_colors)) +
  labs(color = "Models", x = "Time t",
       y = "Average Root Mean Square
Bias across Coefficients") +
  theme_minimal() + theme(legend.position = "none")

pB_RMSE <- ggplot(subset(RMSE_df, Scenario %in% group_B),
                  aes(x = time, y = RMSE, col = model_lab)) +
  geom_point() + geom_line() +
  facet_grid(. ~ Scenario, labeller = label_parsed) +
  scale_color_manual(values = model_colors, breaks = names(model_colors)) +
  labs(color = "Models", x = "Time t",
       y = "Average Root Mean Square
Bias across Coefficients") +
  theme_minimal() + theme(legend.position = "bottom")

pA_RMSE_combined <- ggarrange(pA_RMSE, pEmpty, ncol = 2, widths = c(5, 1.8),
                              common.legend = TRUE, legend = "none")
p_RMSE_double <- ggarrange(pA_RMSE_combined, pB_RMSE, ncol = 1, labels = "AUTO",
                           heights = c(1, 1.35))
ggsave(paste0("00_Figures/", Prefix, "SquaredBias.pdf"), p_RMSE_double, 
       width = 16, height = 8, device = pdf)

# Zoom for RMSE (ALL)
p_RMSE_double_zoom <-  ggarrange(
  ggarrange(
    pA_RMSE + coord_cartesian(ylim = c(0, .03)),
    pEmpty, 
    ncol = 2, widths = c(5, 1.8),
    common.legend = TRUE, legend = "none"
  ), 
  pB_RMSE + coord_cartesian(ylim = c(0, .03)),
  ncol = 1, labels = "AUTO",
  heights = c(1, 1.35)
)
ggsave(paste0("00_Figures/", Prefix, "SquaredBiasZoomed.pdf"), 
       p_RMSE_double_zoom, 
       width = 16, height = 8, device = pdf)

###############################################################################-
# Optional selection-specific plots (selected / not_selected) ------------------
###############################################################################-

if (Add_Selection_Plots) {
  
  # --- Data subsets for selected / not selected scenarios ---
  data_selected         <- data[data$Scenario %in% selected_scenarios, ]
  dataBias_selected     <- dataBias[dataBias$Scenario %in% selected_scenarios, ]
  dataOK_selected       <- dataOK[dataOK$Scenario %in% selected_scenarios, ]
  
  data_not_selected     <- data[!(data$Scenario %in% selected_scenarios), ]
  dataBias_not_selected <- dataBias[!(dataBias$Scenario %in% selected_scenarios), ]
  dataOK_not_selected   <- dataOK[!(dataOK$Scenario %in% selected_scenarios), ]
  
  dataR2_selected       <- dataR2[dataR2$Scenario %in% selected_scenarios, ]
  dataR2_not_selected   <- dataR2[!(dataR2$Scenario %in% selected_scenarios), ]
  
  # ---- Arrange only selected scenarios ----
  p_Scenarios_selected <- ggplot(
    data_selected,
    aes(x = time, y = Effect, col = variable)
  ) +
    geom_line() +
    facet_grid(DV ~ Scenario) +
    theme_minimal() +
    theme(legend.position = "none") 
  
  ggsave(
    paste0("00_Figures/", Prefix, "Scenarios_selected.pdf"), 
    plot   = p_Scenarios_selected, device = pdf, 
    width  = 22 * width_factor_selected, height = 5
  )
  
  # ---- Arrange only non-selected scenarios ----
  p_Scenarios_not_selected <- ggplot(
    data_not_selected,
    aes(x = time, y = Effect, col = variable)
  ) +
    geom_line() +
    facet_grid(DV ~ Scenario) +
    theme_minimal() +
    theme(legend.position = "none") 
  
  ggsave(
    paste0("00_Figures/", Prefix, "Scenarios_not_selected.pdf"), 
    plot   = p_Scenarios_not_selected, device = pdf, 
    width  = 22 * width_factor_not_selected, height = 5
  )
  
  # ---- Split A/B effects of C (SELECTED) ----
  pA_Betas_selected <- ggplot(
    dataOK_selected,
    aes(x = time, y = Effect, col = variable)
  ) +
    geom_line() +
    facet_grid(DV ~ Scenario) +
    theme_minimal() +
    theme(legend.position = "none") 
  
  pB_Betas_selected <- ggplot(
    dataBias_selected,
    aes(x = time, y = Effect, col = variable)
  ) +
    geom_line() +
    facet_grid(DV ~ Scenario) +
    theme_minimal() +
    theme(legend.position = "none") 
  
  # pA_Betas_combined_selected <- ggarrange(
  #   pA_Betas_selected, pEmpty, ncol = 2,
  #   widths = c(5, 1.8),
  #   common.legend = TRUE, legend = "none"
  # )
  pA_Betas_double_selected <- ggarrange(
    pA_Betas_selected, pB_Betas_selected,
    ncol = 1, labels = "AUTO",
    heights = c(1, 1)
  )
  ggsave(
    paste0("00_Figures/", Prefix, "Betas_selected.pdf"), 
    pA_Betas_double_selected, 
    width = 16 * width_factor_selected, height = 8, device = pdf
  )
  
  # ---- Split A/B effects of C (NOT SELECTED) ----
  pA_Betas_not_selected <- ggplot(
    dataOK_not_selected,
    aes(x = time, y = Effect, col = variable)
  ) +
    geom_line() +
    facet_grid(DV ~ Scenario) +
    theme_minimal() +
    theme(legend.position = "none") 
  
  pB_Betas_not_selected <- ggplot(
    dataBias_not_selected,
    aes(x = time, y = Effect, col = variable)
  ) +
    geom_line() +
    facet_grid(DV ~ Scenario) +
    theme_minimal() +
    theme(legend.position = "none") 
  
  pA_Betas_combined_not_selected <- ggarrange(
    pA_Betas_not_selected, pEmpty, ncol = 2,
    widths = c(5, 4.1),
    common.legend = TRUE, legend = "none"
  )
  pA_Betas_double_not_selected <- ggarrange(
    pA_Betas_combined_not_selected,
    pB_Betas_not_selected,
    ncol = 1, labels = "AUTO",
    heights = c(1, 1)
  )
  ggsave(
    paste0("00_Figures/", Prefix, "Betas_not_selected.pdf"), 
    pA_Betas_double_not_selected, 
    width = 16 * width_factor_not_selected, height = 8, device = pdf
  )
  
  # --- R2: SELECTED scenarios ---
  wide_selected  <- subset(wide,  Scenario %in% selected_scenarios)
  lines_selected <- subset(lines_long, Scenario %in% selected_scenarios)
  
  wide_A_selected  <- subset(wide_selected,  Scenario %in% scenarios_A_plain)
  lines_A_selected <- subset(lines_selected, Scenario %in% scenarios_A_plain)
  pA_R2_selected   <- make_R2_plot(wide_A_selected, lines_A_selected)
  
  wide_B_selected  <- subset(wide_selected,  Scenario %in% scenarios_B_plain)
  lines_B_selected <- subset(lines_selected, Scenario %in% scenarios_B_plain)
  pB_R2_selected   <- make_R2_plot(wide_B_selected, lines_B_selected)
  
   pA_R2_combined_selected <- pA_R2_selected + theme(legend.position = "none")
     
   #  ggarrange(
   #  pA_R2_selected, pEmpty, ncol = 2,
   #  common.legend = TRUE, legend = "none"
   # )
  
  pA_R2_double_selected <- ggarrange(
    pA_R2_combined_selected, pB_R2_selected,
    ncol = 1, labels = "AUTO",
    heights = c(1, 1.2)
  )
  
  ggsave(
    paste0("00_Figures/", Prefix, "ScenariosR2_incr_double_selected.pdf"),
    plot = pA_R2_double_selected, device = pdf,
    width = 16 * width_factor_selected, height = 8
  )
  
  # --- R2: NOT SELECTED scenarios ---
  wide_not_selected  <- subset(wide,  Scenario %in% not_selected_scenarios)
  lines_not_selected <- subset(lines_long, Scenario %in% not_selected_scenarios)
  
  wide_A_not_selected  <- subset(wide_not_selected,  Scenario %in% scenarios_A_plain)
  lines_A_not_selected <- subset(lines_not_selected, Scenario %in% scenarios_A_plain)
  pA_R2_not_selected   <- make_R2_plot(wide_A_not_selected, lines_A_not_selected)
  
  wide_B_not_selected  <- subset(wide_not_selected,  Scenario %in% scenarios_B_plain)
  lines_B_not_selected <- subset(lines_not_selected, Scenario %in% scenarios_B_plain)
  pB_R2_not_selected   <- make_R2_plot(wide_B_not_selected, lines_B_not_selected)
  
  pA_R2_combined_not_selected <- ggarrange(
    pA_R2_not_selected, pEmpty, ncol = 2,
    widths = c(5, 4.1),
    common.legend = TRUE, legend = "none"
  )
  pA_R2_double_not_selected <- ggarrange(
    pA_R2_combined_not_selected, pB_R2_not_selected,
    ncol = 1, labels = "AUTO",
    heights = c(1, 1)
  )
  
  ggsave(
    paste0("00_Figures/", Prefix, "ScenariosR2_incr_double_not_selected.pdf"),
    plot = pA_R2_double_not_selected, device = pdf,
    width = 16 * width_factor_not_selected, height = 8
  )
  
  # --- Bias (SELECTED scenarios) ---
  pA_Bias_selected <- ggplot(
    subset(Results, Scenario %in% group_A &
             Scenario_plain %in% selected_scenarios),
    aes(x = time, y = Bias, col = model_lab, group = model_lab)
  ) +
    geom_line() + geom_point() +
    facet_grid(effect_parsed ~ Scenario, 
               labeller = label_parsed, scales = "free_y") +
    scale_color_manual(values = model_colors, breaks = names(model_colors)) +
    labs(color = "Models", x = "Time t", y = "Bias in Coefficient") +
    theme_minimal() + theme(legend.position = "none")
  
  pB_Bias_selected <- ggplot(
    subset(Results, Scenario %in% group_B &
             Scenario_plain %in% selected_scenarios),
    aes(x = time, y = Bias, col = model_lab, group = model_lab)
  ) +
    geom_line() + geom_point() +
    facet_grid(effect_parsed ~ Scenario, 
               labeller = label_parsed, scales = "free_y") +
    scale_color_manual(values = model_colors, breaks = names(model_colors)) +
    labs(color = "Models", x = "Time t", y = "Bias in Coefficient") +
    theme_minimal() + theme(legend.position = "bottom")
  
  # pA_Bias_combined_selected <- ggarrange(
  #   pA_Bias_selected, pEmpty, ncol = 2,
  #   widths = c(5, 1.8),
  #   common.legend = TRUE, legend = "none"
  # )
  p_Bias_double_selected <- ggarrange(
    pA_Bias_selected, pB_Bias_selected,
    ncol = 1, labels = "AUTO", 
    heights = c(1, 1)
  )
  ggsave(
    paste0("00_Figures/", Prefix, "Bias_selected.pdf"), 
    p_Bias_double_selected, 
    width = 16 * width_factor_selected, height = vertical_scaling*12, device = pdf
  )
  
  # --- Bias (NOT SELECTED scenarios) ---
  pA_Bias_not_selected <- ggplot(
    subset(Results, Scenario %in% group_A &
             Scenario_plain %in% not_selected_scenarios),
    aes(x = time, y = Bias, col = model_lab, group = model_lab)
  ) +
    geom_line() + geom_point() +
    facet_grid(effect_parsed ~ Scenario, 
               labeller = label_parsed, scales = "free_y") +
    scale_color_manual(values = model_colors, breaks = names(model_colors)) +
    labs(color = "Models", x = "Time t", y = "Bias in Coefficient") +
    theme_minimal() + theme(legend.position = "none")
  
  pB_Bias_not_selected <- ggplot(
    subset(Results, Scenario %in% group_B &
             Scenario_plain %in% not_selected_scenarios),
    aes(x = time, y = Bias, col = model_lab, group = model_lab)
  ) +
    geom_line() + geom_point() +
    facet_grid(effect_parsed ~ Scenario, 
               labeller = label_parsed, scales = "free_y") +
    scale_color_manual(values = model_colors, breaks = names(model_colors)) +
    labs(color = "Models", x = "Time t", y = "Bias in Coefficient") +
    theme_minimal() + theme(legend.position = "bottom")
  
  pA_Bias_combined_not_selected <- ggarrange(
    pA_Bias_not_selected, pEmpty, ncol = 2,
    widths = c(5, 3.9),
    common.legend = TRUE, legend = "none"
  )
  p_Bias_double_not_selected <- ggarrange(
    pA_Bias_combined_not_selected,
    pB_Bias_not_selected,
    ncol = 1, labels = "AUTO", 
    heights = c(1, 1.35)
  )
  ggsave(
    paste0("00_Figures/", Prefix, "Bias_not_selected.pdf"), 
    p_Bias_double_not_selected, 
    width = 16 * width_factor_not_selected, height = vertical_scaling*12, device = pdf
  )
  
  # --- Expected Z (SELECTED scenarios) ---
  pA_Z_selected <- ggplot(
    subset(Results, Scenario %in% group_A &
             Scenario_plain %in% selected_scenarios),
    aes(x = time, y = z, col = model_lab, group = model_lab)
  ) +
    geom_line() + geom_point() +
    facet_grid(effect_parsed ~ Scenario,
               labeller = label_parsed, scales = "free_y") +
    scale_color_manual(values = model_colors, breaks = names(model_colors)) +
    labs(color = "Models", x = "Time t",
         y = "Expected Z-Value for N = 1000") +
    theme_minimal() + theme(legend.position = "none")
  
  pB_Z_selected <- ggplot(
    subset(Results, Scenario %in% group_B &
             Scenario_plain %in% selected_scenarios),
    aes(x = time, y = z, col = model_lab, group = model_lab)
  ) +
    geom_line() + geom_point() +
    facet_grid(effect_parsed ~ Scenario,
               labeller = label_parsed, scales = "free_y") +
    scale_color_manual(values = model_colors, breaks = names(model_colors)) +
    labs(color = "Models", x = "Time t",
         y = "Expected Z-Value for N = 1000") +
    theme_minimal() + theme(legend.position = "bottom")
  
  # pA_Z_combined_selected <- ggarrange(
  #   pA_Z_selected, pEmpty, ncol = 2,
  #   widths = c(5, 1.8),
  #   common.legend = TRUE, legend = "none"
  # )
  p_Z_double_selected <- ggarrange(
    pA_Z_selected, pB_Z_selected,
    ncol = 1, labels = "AUTO",
    heights = c(1, 1.35)
  )
  ggsave(
    paste0("00_Figures/", Prefix, "Z_selected.pdf"), 
    p_Z_double_selected, 
    width = 16 * width_factor_selected, height = vertical_scaling*12, device = pdf
  )
  
  # --- Expected Z (NOT SELECTED scenarios) ---
  pA_Z_not_selected <- ggplot(
    subset(Results, Scenario %in% group_A &
             Scenario_plain %in% not_selected_scenarios),
    aes(x = time, y = z, col = model_lab, group = model_lab)
  ) +
    geom_line() + geom_point() +
    facet_grid(effect_parsed ~ Scenario,
               labeller = label_parsed, scales = "free_y") +
    scale_color_manual(values = model_colors, breaks = names(model_colors)) +
    labs(color = "Models", x = "Time t",
         y = "Expected Z-Value for N = 1000") +
    theme_minimal() + theme(legend.position = "none")
  
  pB_Z_not_selected <- ggplot(
    subset(Results, Scenario %in% group_B &
             Scenario_plain %in% not_selected_scenarios),
    aes(x = time, y = z, col = model_lab, group = model_lab)
  ) +
    geom_line() + geom_point() +
    facet_grid(effect_parsed ~ Scenario,
               labeller = label_parsed, scales = "free_y") +
    scale_color_manual(values = model_colors, breaks = names(model_colors)) +
    labs(color = "Models", x = "Time t",
         y = "Expected Z-Value for N = 1000") +
    theme_minimal() + theme(legend.position = "bottom")
  
  pA_Z_combined_not_selected <- ggarrange(
    pA_Z_not_selected, pEmpty, ncol = 2,
    widths = c(5, 4.1),
    common.legend = TRUE, legend = "none"
  )
  p_Z_double_not_selected <- ggarrange(
    pA_Z_combined_not_selected, pB_Z_not_selected,
    ncol = 1, labels = "AUTO",
    heights = c(1, 1.35)
  )
  ggsave(
    paste0("00_Figures/", Prefix, "Z_not_selected.pdf"), 
    p_Z_double_not_selected, 
    width = 16 * width_factor_not_selected, height = vertical_scaling*12, device = pdf
  )
  
  # --- Relative Bias (SELECTED scenarios) ---
  pA_RelBias_selected <- ggplot(
    subset(Results_rel, Scenario %in% group_A &
             Scenario_plain %in% selected_scenarios),
    aes(x = time, y = RelBiasPercent, col = model_lab, group = model_lab)
  ) +
    geom_line() + geom_point() +
    facet_grid(effect_parsed ~ Scenario,
               labeller = label_parsed) +
    scale_color_manual(values = model_colors, breaks = names(model_colors)) +
    labs(color = "Models", x = "Time t",
         y = "Relative Bias\nin Coefficient [in %]") +
    theme_minimal() + theme(legend.position = "none")
  
  pB_RelBias_selected <- ggplot(
    subset(Results_rel, Scenario %in% group_B &
             Scenario_plain %in% selected_scenarios),
    aes(x = time, y = RelBiasPercent, col = model_lab, group = model_lab)
  ) +
    geom_line() + geom_point() +
    facet_grid(effect_parsed ~ Scenario,
               labeller = label_parsed) +
    scale_color_manual(values = model_colors, breaks = names(model_colors)) +
    labs(color = "Models", x = "Time t",
         y = "Relative Bias\nin Coefficient [in %]") +
    theme_minimal() + theme(legend.position = "bottom")
  
  # pA_RelBias_combined_selected <- ggarrange(
  #   pA_RelBias_selected, pEmpty, ncol = 2,
  #   widths = c(5, 1.8),
  #   common.legend = TRUE, legend = "none"
  # )
  p_RelBias_double_selected <- ggarrange(
    pA_RelBias_selected, pB_RelBias_selected,
    ncol = 1, labels = "AUTO",
    heights = c(1, 1.35)
  )
  ggsave(
    paste0("00_Figures/", Prefix, "RelBias_selected.pdf"), 
    p_RelBias_double_selected, 
    width = 16 * width_factor_selected, height = vertical_scaling*12, device = pdf
  )
  
  # --- Relative Bias (NOT SELECTED scenarios) ---
  pA_RelBias_not_selected <- ggplot(
    subset(Results_rel, Scenario %in% group_A &
             Scenario_plain %in% not_selected_scenarios),
    aes(x = time, y = RelBiasPercent, col = model_lab, group = model_lab)
  ) +
    geom_line() + geom_point() +
    facet_grid(effect_parsed ~ Scenario,
               labeller = label_parsed) +
    scale_color_manual(values = model_colors, breaks = names(model_colors)) +
    labs(color = "Models", x = "Time t",
         y = "Relative Bias\nin Coefficient [in %]") +
    theme_minimal() + theme(legend.position = "none")
  
  pB_RelBias_not_selected <- ggplot(
    subset(Results_rel, Scenario %in% group_B &
             Scenario_plain %in% not_selected_scenarios),
    aes(x = time, y = RelBiasPercent, col = model_lab, group = model_lab)
  ) +
    geom_line() + geom_point() +
    facet_grid(effect_parsed ~ Scenario,
               labeller = label_parsed) +
    scale_color_manual(values = model_colors, breaks = names(model_colors)) +
    labs(color = "Models", x = "Time t",
         y = "Relative Bias\nin Coefficient [in %]") +
    theme_minimal() + theme(legend.position = "bottom")
  
  pA_RelBias_combined_not_selected <- ggarrange(
    pA_RelBias_not_selected, pEmpty, ncol = 2,
    widths = c(5, 4.1),
    common.legend = TRUE, legend = "none"
  )
  p_RelBias_double_not_selected <- ggarrange(
    pA_RelBias_combined_not_selected,
    pB_RelBias_not_selected,
    ncol = 1, labels = "AUTO",
    heights = c(1, 1.35)
  )
  ggsave(
    paste0("00_Figures/", Prefix, "RelBias_not_selected.pdf"), 
    p_RelBias_double_not_selected, 
    width = 16 * width_factor_not_selected, height = vertical_scaling*12, device = pdf
  )
  
  # -------------------------------------------------------------------------- #
  # Relative Bias (SELECTED) – Zoomed to [-100, 100] and [-30, 30]
  # -------------------------------------------------------------------------- #
  
  p_RelBias_double_selected_100 <- ggarrange(
    #ggarrange(
      pA_RelBias_selected + coord_cartesian(ylim = c(-100, 100)),
      #pEmpty,
      #ncol = 2,
      #widths = c(5, 1.8),
      #common.legend = TRUE,
      #legend = "none"
    #),
    pB_RelBias_selected + coord_cartesian(ylim = c(-100, 100)),
    ncol = 1,
    labels = "AUTO",
    heights = c(1, 1.35)
  )
  
  ggsave(
    paste0("00_Figures/", Prefix, "RelBiasZoomed100_selected.pdf"),
    p_RelBias_double_selected_100,
    width  = 16 * width_factor_selected, height = vertical_scaling*12,
    device = pdf
  )
  
  p_RelBias_double_selected_30 <- ggarrange(
    #ggarrange(
      pA_RelBias_selected + coord_cartesian(ylim = c(-30, 30)),
     # pEmpty,
    #  ncol = 2,
    #  widths = c(5, 1.8),
    #  common.legend = TRUE,
    #  legend = "none"
    #),
    pB_RelBias_selected + coord_cartesian(ylim = c(-30, 30)),
    ncol = 1,
    labels = "AUTO",
    heights = c(1, 1.35)
  )
  
  ggsave(
    paste0("00_Figures/", Prefix, "RelBiasZoomed30_selected.pdf"),
    p_RelBias_double_selected_30,
    width  = 16 * width_factor_selected, height = vertical_scaling*12,
    device = pdf
  )
  
  # -------------------------------------------------------------------------- #
  # Relative Bias (NOT SELECTED) – Zoomed to [-100, 100] and [-30, 30]
  # -------------------------------------------------------------------------- #
  
  p_RelBias_double_not_selected_100 <- ggarrange(
    ggarrange(
      pA_RelBias_not_selected + coord_cartesian(ylim = c(-100, 100)),
      pEmpty,
      ncol = 2,
      widths = c(5, 4.1),
      common.legend = TRUE,
      legend = "none"
    ),
    pB_RelBias_not_selected + coord_cartesian(ylim = c(-100, 100)),
    ncol = 1,
    labels = "AUTO",
    heights = c(1, 1.35)
  )
  
  ggsave(
    paste0("00_Figures/", Prefix, "RelBiasZoomed100_not_selected.pdf"),
    p_RelBias_double_not_selected_100,
    width  = 16 * width_factor_not_selected, height = vertical_scaling*12,
    device = pdf
  )
  
  p_RelBias_double_not_selected_30 <- ggarrange(
    ggarrange(
      pA_RelBias_not_selected + coord_cartesian(ylim = c(-30, 30)),
      pEmpty,
      ncol = 2,
      widths = c(5, 4.1),
      common.legend = TRUE,
      legend = "none"
    ),
    pB_RelBias_not_selected + coord_cartesian(ylim = c(-30, 30)),
    ncol = 1,
    labels = "AUTO",
    heights = c(1, 1.35)
  )
  
  ggsave(
    paste0("00_Figures/", Prefix, "RelBiasZoomed30_not_selected.pdf"),
    p_RelBias_double_not_selected_30,
    width  = 16 * width_factor_not_selected, height = vertical_scaling*12,
    device = pdf
  )
  
  # --- RMSE (SELECTED scenarios) ---
  pA_RMSE_selected <- ggplot(
    subset(RMSE_df, Scenario %in% group_A &
             Scenario_plain %in% selected_scenarios),
    aes(x = time, y = RMSE, col = model_lab)
  ) +
    geom_point() + geom_line() +
    facet_grid(. ~ Scenario, labeller = label_parsed) +
    scale_color_manual(values = model_colors, breaks = names(model_colors)) +
    labs(color = "Models", x = "Time t", 
         y = "Average Root Mean Square
Bias across Coefficients") +
    theme_minimal() + theme(legend.position = "none")
  
  pB_RMSE_selected <- ggplot(
    subset(RMSE_df, Scenario %in% group_B &
             Scenario_plain %in% selected_scenarios),
    aes(x = time, y = RMSE, col = model_lab)
  ) +
    geom_point() + geom_line() +
    facet_grid(. ~ Scenario, labeller = label_parsed) +
    scale_color_manual(values = model_colors, breaks = names(model_colors)) +
    labs(color = "Models", x = "Time t", 
         y = "Average Root Mean Square
Bias across Coefficients") +
    theme_minimal() + theme(legend.position = "bottom")
  
  # pA_RMSE_combined_selected <- ggarrange(
  #   pA_RMSE_selected, pEmpty, ncol = 2,
  #   widths = c(5, 1.8),
  #   common.legend = TRUE, legend = "none"
  # )
  p_RMSE_double_selected <- ggarrange(
    pA_RMSE_selected, pB_RMSE_selected,
    ncol = 1, labels = "AUTO",
    heights = c(1, 1.35)
  )
  ggsave(
    paste0("00_Figures/", Prefix, "SquaredBias_selected.pdf"), 
    p_RMSE_double_selected, 
    width = 16 * width_factor_selected, height = 8, device = pdf
  )
  
  # --- RMSE (NOT SELECTED scenarios) ---
  pA_RMSE_not_selected <- ggplot(
    subset(RMSE_df, Scenario %in% group_A &
             Scenario_plain %in% not_selected_scenarios),
    aes(x = time, y = RMSE, col = model_lab)
  ) +
    geom_point() + geom_line() +
    facet_grid(. ~ Scenario, labeller = label_parsed) +
    scale_color_manual(values = model_colors, breaks = names(model_colors)) +
    labs(color = "Models", x = "Time t", 
         y = "Average Root Mean Square
Bias across Coefficients") +
    theme_minimal() + theme(legend.position = "none")
  
  pB_RMSE_not_selected <- ggplot(
    subset(RMSE_df, Scenario %in% group_B &
             Scenario_plain %in% not_selected_scenarios),
    aes(x = time, y = RMSE, col = model_lab)
  ) +
    geom_point() + geom_line() +
    facet_grid(. ~ Scenario, labeller = label_parsed) +
    scale_color_manual(values = model_colors, breaks = names(model_colors)) +
    labs(color = "Models", x = "Time t", 
         y = "Average Root Mean Square
Bias across Coefficients") +
    theme_minimal() + theme(legend.position = "bottom")
  
  pA_RMSE_combined_not_selected <- ggarrange(
    pA_RMSE_not_selected, pEmpty, ncol = 2,
    widths = c(5, 4.1),
    common.legend = TRUE, legend = "none"
  )
  p_RMSE_double_not_selected <- ggarrange(
    pA_RMSE_combined_not_selected, pB_RMSE_not_selected,
    ncol = 1, labels = "AUTO",
    heights = c(1, 1.35)
  )
  ggsave(
    paste0("00_Figures/", Prefix, "SquaredBias_not_selected.pdf"), 
    p_RMSE_double_not_selected, 
    width = 16 * width_factor_not_selected, height = 8, device = pdf
  )
}

