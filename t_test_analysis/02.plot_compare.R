# 02.plot_compare.R
#
# Author: Bob Verity
# Date: 2023-03-03
#
# Purpose:
# Read in results of 01.simulate_t_test.R, produce plot and save to file.
#
# ------------------------------------------------------------------

# load packages
library(tidyverse)

# load DRpower package from local version
#devtools::load_all("/Users/rverity/Dropbox/Bob/Work/My Programs/Sample design/DRpower")

# define plotting parameters
line_size <- 0.5

# read data from file
df_pow <- readRDS("t_test_analysis/data/df_pow.rds")

# get into long format
df_long <- df_pow %>%
  dplyr::select(n_clust, p, t_power, t_analytical, DR_power) %>%
  pivot_longer(-c(n_clust, p), names_to = "method", values_to = "power") %>%
  mutate(method = replace(method, method == "t_power", "t-test, empirical"),
         method = replace(method, method == "t_analytical", "t-test, analytical"),
         method = replace(method, method == "DR_power", "DRpower"),
         method = factor(method, levels = c("DRpower", "t-test, empirical", "t-test, analytical")),
         n_clust_title = sprintf("clusters = %s", n_clust),
         n_clust_title = factor(n_clust_title, levels = sprintf("clusters = %s", c(2, 10))))

# plot
plot1 <- df_long %>%
  #dplyr::filter(method != "t-test, analytical") %>%
  ggplot() + theme_bw() +
  geom_rect(data = mtcars[1,], aes(xmin = 0.05, xmax = 0.2, ymin = 0.8, ymax = 1), fill = "dodgerblue2", alpha = 0.2) +
  geom_line(aes(x = p, y = power, col = method), size = line_size) +
  ylim(c(0, 1)) +
  facet_wrap(~n_clust_title) +
  scale_color_discrete(name = "Analysis method") +
  xlab("Prevalence of pfhrp2 deletions") + ylab("Power")
plot1

# save plot to file
ggsave(plot = plot1, filename = "t_test_analysis/outputs/power_compare.pdf", width = 8, height = 4)
