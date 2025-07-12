# 01.make_figure.R
#
# Author: Bob Verity
# Date: 2024-06-10
#
# Inputs: (none)
#
# Outputs: outputs/CI_Z_plot.png
#
# Purpose:
# Creates a plot that illustrates the connection between the confidence interval
# (CI) based approach and the z-test.
#
# ------------------------------------------------------------------
library(tidyverse)

data.frame(y = as.factor(1:3), x = c(5.5, 7, 3.5)) |>
  ggplot() + theme_classic() +
  geom_errorbar(aes(x = x, y = y, xmin = x - 1, xmax = x + 1), width = 0.1) +
  geom_point(aes(x = x, y = y)) +
  geom_vline(xintercept = 5, linetype = "dashed") +
  xlim(c(0, 10)) +
  xlab("Prevalence of deletions (%)") + ylab("") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

ggsave("01.figure_CI_Ztest/outputs/CI_Z_plot.png", scale = 0.6)
