# 04.plot_ICC.R
#
# Author: Bob Verity
# Date: 2025-07-11
#
# Inputs:
# - 03.ICC_posteriors/outputs/dat_ICCfilter.rds
#
# Outputs:
# - 03.ICC_posteriors/outputs/historical_ICC.png
#
# Purpose:
# Reads in filtered historical data, estimates the full posterior distribution
# of ICC for each distinct region, and plots these results. Figures are saved to
# file as an image.
#
# ------------------------------------------------------------------
library(cowplot)
library(tidyverse)
#devtools::install_github("mrc-ide/DRpower@v1.0.3")
library(DRpower)

# read in data
dat <- readRDS("03.ICC_posteriors/outputs/dat_ICCfilter.rds")

# create unique group based on ADMIN1 and year
dat <- dat |>
  mutate(group = sprintf("%s: %s, %s", COUNTRY_NAME, ADMIN1_NAME, YEAR_START))

# split by group
dat_list <- split(dat, f = dat$group)

# loop through studies and get full posterior in each case
x <- seq(0, 1, l = 1001)
l <- list()
for (i in 1:length(dat_list)) {
  
  # subset to this study and get raw posterior ICC
  post_df <- get_ICC(n = dat_list[[i]]$HRP2_NUM_DELETION,
                     N = dat_list[[i]]$HRP2_TESTED,
                     prior_ICC_shape2 = 1,
                     n_intervals = 100,
                     post_full_on = TRUE,
                     post_full_breaks = x)
  
  # save results
  y <- post_df$post_full[[1]]
  l[[i]] <- data.frame(col = "data",
                       group = dat_list[[i]]$group[1],
                       x = x,
                       y = y / sum(y))
}
df1 <- l %>%
  bind_rows()

# add combined posterior
post_comb <- df1 |>
  pivot_wider(names_from = group, values_from = y) |>
  select(-col, -x) |>
  t() |>
  log() |>
  colSums() |>
  exp()

df1 <- rbind(df1, data.frame(col = "posterior",
                             group = "Combined posterior",
                             x = x,
                             y = post_comb / sum(post_comb)))

# add default prior
y <- dbeta(x, shape1 = 1, shape2 = 9)
df1 <- rbind(df1, data.frame(col = "prior",
                             group = "Default prior",
                             x = x,
                             y = y / sum(y)))


# plot
group_names <- c("India: Odisha, 2013", "Peru: Loreto, 2011", "South Sudan: Yambio County, 2020",
                 "Ethiopia: Tigray, 2017", "Ethiopia: Benshangul-Gumaz, 2018", "Ethiopia: Amhara, 2017", 
                 "Burkina Faso: Centre-Ouest, 2020", "Eritrea: Gash Barka, 2019", "Ethiopia: Oromia, 2022",
                 "Ethiopia: Southern Nations, Nationalities, 2022")

plot1 <- df1 %>%
  filter(!group %in% c("Default prior", "Combined posterior")) |>
  mutate(group = factor(group, levels = group_names)) %>%
  ggplot() + theme_bw() +
  geom_ribbon(aes(x = x, ymin = 0, ymax = y, fill = group, col = group), alpha = 0.2) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.04), expand = c(0, 0)) +
  xlab("Intra-cluster correlation") + ylab("Posterior probability density") +
  theme(legend.title = element_blank(), legend.position = c(1, 1), legend.justification = c(1,1),
        legend.background = element_rect(linewidth = 0.3, colour = "black")) +
  theme(plot.margin = unit(c(.2,.5,.2,.2),"cm"))
plot1

plot2 <- df1 %>%
  filter(group %in% c("Default prior", "Combined posterior")) |>
  ggplot() + theme_bw() +
  geom_ribbon(aes(x = x, ymin = 0, ymax = y, fill = group, col = group), alpha = 0.5) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.04), expand = c(0, 0)) +
  xlab("Intra-cluster correlation") + ylab("Posterior probability density") +
  theme(legend.title = element_blank(), legend.position = c(1, 1), legend.justification = c(1,1),
        legend.background = element_rect(size = 0.3, colour = "black")) +
  theme(plot.margin = unit(c(.2,.5,.2,.2),"cm"))
plot2

# combine plots
plot3 <- cowplot::plot_grid(plot1 + labs(subtitle = "a)"),
                            plot2 + labs(subtitle = "b)"))
plot3

# get MAP estimate from combined posterior
x[which.max(post_comb)]

# save plot to file
ggsave(filename = "03.ICC_posteriors/outputs/historical_ICC.png", width = 8, height = 3.5,
       scale = 1.2)

