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
df_CrI_combined <- readRDS("historical_analysis/data/df_CrI_combined.rds")

# plot
df_CrI_combined %>%
  mutate(plot_name = factor(plot_name, levels = rev(plot_name))) %>%
  arrange(plot_name) %>%
  ggplot() + theme_bw() +
  geom_errorbar(aes(x = plot_name, ymin = CrI_lower, ymax = CrI_upper), width = 0.2) +
  geom_point(aes(x = plot_name, y = post_median)) +
  ylim(c(0, 1)) +
  coord_flip() + xlab("") + ylab("Intra-cluster correlation") +
  theme(axis.text.y = element_text(face = c("bold", rep("plain", 6))))

# save plot to file
ggsave("historical_analysis/outputs/ICC_CrI.png", width = 7, height = 2.5)
