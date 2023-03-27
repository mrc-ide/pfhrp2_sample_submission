# 02.plot_ICC.R
#
# Author: Bob Verity
# Date: 2023-02-28
#
# Purpose:
# Read in ICC credible intervals and produce plot.
#
# ------------------------------------------------------------------

# read from file
df_CrI <- readRDS("historical_analysis/data/df_CrI.rds")

# plot
df_CrI %>%
  mutate(plot_name = factor(plot_name, levels = rev(plot_name))) %>%
  arrange(plot_name) %>%
  ggplot() + theme_bw() +
  geom_pointrange(aes(x = plot_name, y = median, ymin = lower, ymax = upper)) +
  ylim(c(0, 1)) +
  coord_flip() + xlab("") + ylab("Intra-cluster correlation") +
  theme(axis.text.y = element_text(face = c("bold", rep("plain", 5))))

# save plot to file
ggsave("historical_analysis/outputs/ICC_CrI.png", width = 7, height = 2.5)
