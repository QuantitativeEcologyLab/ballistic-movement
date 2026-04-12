# script is designed to show how the ballistic motion model plays out
# for a single scenario

#..............................................................................
# load libraries ----

library(tidyverse)
library(patchwork)
library(scico)
library(viridis)

#..............................................................................

#load data

details_df <- bind_rows(readRDS("simulations/prey_results/patches/250p_prey_details.Rds"))

#observe lv over generation time 

#calculate stable lv
stablelv <- details_df %>% 
  filter(generation >= max(generation) - 100) %>% 
  summarise(mean_lv = mean(lv)) %>% 
  as.numeric()

#plot lv over course of simulation
a <-
  ggplot() +
  ggtitle("A") +
  geom_point(data = details_df, aes(x = generation, y = lv, col = generation), 
             # col = "#B5989E", 
             alpha = 0.1, size = 0.5, stroke = NA) +
  scale_color_viridis(option = "magma") +
  stat_summary(data = details_df, aes(x = generation, y = lv), fun = "mean", geom = "line", lwd = 0.5, col = "black") +
  stat_summary(data = details_df, aes(x = generation, y = lv), fun = "mean", geom = "line", lwd = 0.3, col = "white") +
  geom_hline(yintercept = stablelv, col = "black", lwd = 0.5) +
  geom_hline(yintercept = stablelv, lty = "dashed", col = "white", lwd = 0.3) +
  labs(x = "Generation", y = expression(bold(l[v] (m))), col = "Generation") +
  scale_x_continuous(expand = c(0.01,0.01)) +
  scale_y_continuous(limits = c(0, (max(details_df$lv)*0.8))) +
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
        legend.title = element_text(size = 8, family = "sans", face = "bold"),
        legend.text = element_text(size = 6, family = "sans"),
        legend.position = "right",
        legend.key.width = unit(0.3, 'cm'),
        legend.key.height = unit(1.3, 'cm'),
        panel.background = element_rect(fill = "transparent"))

#observe lv and speed relationship
b <-
  ggplot() +
  ggtitle("B") +
  geom_point(data = details_df, aes(x = lv, y = speed, color = generation), size = 0.5, alpha = 0.5, stroke = NA) +
  labs(y = "Speed (m/s)", x = expression(bold(l[v] (m))), col = "Generation") +
  scale_x_continuous(expand = c(0.01,0.01), limits = c(0, (max(details_df$lv)*0.8))) +
  scale_y_continuous(limits = c(0, (max(details_df$speed)*0.8))) +
  scale_color_viridis(option = "magma") +
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
        legend.title = element_text(size = 8, family = "sans", face = "bold"),
        legend.text = element_text(size = 6, family = "sans"),
        legend.position = "right",
        legend.key.width = unit(0.3, 'cm'),
        legend.key.height = unit(1.3, 'cm'),
        panel.background = element_rect(fill = "transparent"))

#observe calories gained by speed
c <-
  ggplot() +
  ggtitle("C") +
  geom_point(data = details_df, aes(x = speed, y = cal_net, color = generation), size = 0.5, alpha = 0.5, stroke = NA) +
  labs(x = "Speed (m/s)", y = "Net Calories", col = "Generation") +
  scale_x_continuous(expand = c(0.01,0.01), limits = c(0, (max(details_df$speed)*0.8))) + 
  scale_color_viridis(option = "magma") +
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
        legend.title = element_text(size = 8, family = "sans", face = "bold"),
        legend.text = element_text(size = 6, family = "sans"),
        legend.position = "right",
        legend.key.width = unit(0.3, 'cm'),
        legend.key.height = unit(1.3, 'cm'),
        panel.background = element_rect(fill = "transparent"))

#observe calories gained by lv
d <-
  ggplot() +
  ggtitle("D") +
  geom_point(data = details_df, aes(x = lv, y = cal_net, color = generation), size = 0.5, alpha = 0.5, stroke = NA) +
  geom_vline(xintercept = stablelv, lty = "dashed", col = "grey20", lwd = 0.3) +
  labs(x = expression(bold(l[v] (m))), y = "Net Calories", col = "Generation") +
  scale_x_continuous(expand = c(0.01,0.01), limits = c(0, (max(details_df$lv)*0.8))) +
  scale_color_viridis(option = "magma") +
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
        legend.title = element_text(size = 8, family = "sans", face = "bold"),
        legend.text = element_text(size = 6, family = "sans"),
        legend.position = "right",
        legend.key.width = unit(0.3, 'cm'),
        legend.key.height = unit(1.3, 'cm'),
        panel.background = element_rect(fill = "transparent"))


FIG <- (a + b + c + d) + plot_layout(guides = "collect") & 
  theme(legend.position = "right")

ggsave(FIG, file = "figures/maintext/simulation-overview.png", 
       width = 9, height = 4.5, units = "in", bg = "white", dpi = 600)
