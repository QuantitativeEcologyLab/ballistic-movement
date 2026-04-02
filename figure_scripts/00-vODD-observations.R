# script for the observation section of vODD

#..............................................................................
# load libraries ----

library(tidyverse)
library(patchwork)
library(scico)
library(viridis)

#..............................................................................

#load data

details_df <- bind_rows(readRDS("simulations/prey_results/patches/510p_prey_details.Rds"))

#identify evolutionary stable strategies
a <-
  ggplot() +
  geom_point(data = details_df, aes(x = generation, y = lv), col = "#CC9BA5", alpha = 0.3, size = 0.5, stroke = NA) +
  stat_summary(data = details_df, aes(x = generation, y = lv), fun = "mean", geom = "line", linewidth = 0.5) +
  labs(x = "Generation", y = expression(bold(l[v]))) +
  scale_x_continuous(expand = c(0.01,0.01)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=8, family = "sans", face = "bold"),
        axis.title.x = element_text(size=8, family = "sans", face = "bold"),
        axis.text.y = element_text(size=6, family = "sans"),
        axis.text.x  = element_text(size=6, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 10, family = "sans", face = "bold"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
        legend.position = "right",
        legend.text = element_text(size = 10, family = "sans"),
        legend.title = element_text(size = 12, family = "sans", face = "bold"),
        legend.key.size = unit(0.2, "cm"),
        legend.spacing.y = unit(0.1, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent"))

ggsave(a, file = "figures/maintext/vODD_figures/evo-stable.png", width = 5, height = 2.5, units = "in", bg = "white", dpi = 600)

details_df <- bind_rows(readRDS("simulations/prey_results/patches/250p_prey_details.Rds"))

# optimization of movement parameters
#observe calories gained by lv
b <-
  ggplot() +
  geom_point(data = details_df, aes(x = lv, y = cal_net, col = generation), size = 0.1) +
  labs(x = expression(bold(l[v])), y = "Net Calories Gained", col = "Generation") +
  scale_x_continuous(expand = c(0.01,0.01), limits = c(0, 750)) +
  scale_color_gradientn(colors = c("#401D1F", "#784748", "#95757cff", "#CC9BA5", "#E0C5CC")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=8, family = "sans", face = "bold"),
        axis.title.x = element_text(size=8, family = "sans", face = "bold"),
        axis.text.y = element_text(size=6, family = "sans"),
        axis.text.x  = element_text(size=6, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 10, family = "sans", face = "bold"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
        legend.position = "right",
        legend.text = element_text(size = 4, family = "sans"),
        legend.title = element_text(size = 6, family = "sans", face = "bold"),
        legend.key.size = unit(0.5, "cm"),
        legend.spacing.y = unit(0.2, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent"))

#observe calories gained by speed
c <-
  ggplot() +
  geom_point(data = details_df, aes(x = speed, y = cal_net, col = generation), size = 0.1) +
  labs(x = "Speed m/s", y = "Net Calories Gained", col = "Generation") +
  scale_x_continuous(expand = c(0.01,0.01), limits = c(0,4)) + 
  scale_color_gradientn(colors = c("#401D1F", "#784748", "#95757cff", "#CC9BA5", "#E0C5CC")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=8, family = "sans", face = "bold"),
        axis.title.x = element_text(size=8, family = "sans", face = "bold"),
        axis.text.y = element_text(size=6, family = "sans"),
        axis.text.x  = element_text(size=6, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 10, family = "sans", face = "bold"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
        legend.position = "right",
        legend.text = element_text(size = 4, family = "sans"),
        legend.title = element_text(size = 6, family = "sans", face = "bold"),
        legend.key.size = unit(0.5, "cm"),
        legend.spacing.y = unit(0.2, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent"))


FIG <- (b + c) + plot_layout(guides = "collect")

ggsave(FIG, file = "figures/maintext/vODD_figures/move-opt.png", width = 6.5, height = 2.5, units = "in", bg = "white", dpi = 600)
