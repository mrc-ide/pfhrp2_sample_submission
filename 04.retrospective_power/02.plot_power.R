# 02.plot_power.R
#
# Author: Bob Verity
# Date: 2025-07-11
#
# Inputs:
# - 04.retrospective_power/outputs/retro_power.rds
#
# Outputs:
# - 04.retrospective_power/outputs/power_hist.png
#
# Purpose:
# Plots power from retrospective analysis.
#
# ------------------------------------------------------------------
library(tidyverse)
library(forcats)

# read in data
dat <- readRDS("04.retrospective_power/outputs/retro_power.rds")

# make the plot
dat |>
  rename(z_test = pow_z_Dplugin,
         DRpower = pow_Bayes) |>
  select(z_test, DRpower) |>
  pivot_longer(cols = 1:2, names_to = "Method") |>
  mutate(Method = factor(Method, levels = c("z_test", "DRpower")),
         Method = fct_recode(Method,
                             "Z-test" = "z_test",
                             "DRpower" = "DRpower")) |>
  ggplot() + theme_bw() +
  geom_histogram(aes(x = value, fill = Method), breaks = seq(0, 100, 5), colour = "black") +
  geom_vline(xintercept = 80, linetype = "dashed") +
  facet_wrap(~Method) +
  scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 45), expand = c(0, 0)) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(size = 14, face = "bold")
  ) + theme(
    plot.margin = margin(t = 10, r = 15, b = 10, l = 10),  # bottom margin increased
    axis.text.x = element_text(margin = margin(t = 5))    # space between x-axis text and axis
  ) + 
  theme(panel.spacing = unit(1, "lines")) +
  xlab("Power (%)") + ylab("Studies")

# save to file
ggsave(filename = "04.retrospective_power/outputs/power_hist.png", height = 5, width = 9)
