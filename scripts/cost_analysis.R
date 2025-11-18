# Preamble ----

# Set the working directory
setwd("~/hdrive/GitHub/ballistic-movement")

library(ggplot2)
library(tidyverse)
library(dplyr)
library(mgcv)
library(scico)
library(gridExtra)
library(scales) # for scientific legend labels
library(ggspatial) # for scale bars 
library(grid)

source("scripts/functions.R")

#----------------------------------------------------------------------
# comparing with and without energetic costs----
#----------------------------------------------------------------------
## without costs associated with movement
load('simulations/supplementary/105500g_nocost_prey_details.Rda')
prey_details_df_nocosts <- do.call(rbind, prey_details)

## with costs associated with movement
load('simulations/prey_results/105500g_prey_details.Rda')
prey_details_df_costs <- do.call(rbind, prey_details)
### remove the costs variable from the data set so that the no costs/costs match
prey_details_df_costs <- subset(prey_details_df_costs, select = - costs)

# create figures ---------------------------------------------------------------

# NO ENERGETIC COSTS FIGURES (A-D)

## offspring versus speed
a <-
  ggplot() +
  ggtitle("A") +
  geom_point(data = prey_details_df_nocosts, aes(x = speed, y = offspring), size = 0.1, color = "#24262d") +
  labs(x = "Speed (m/s)", y = "Number of Offspring") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=9, family = "sans", face = "bold", margin = margin(r = 4)),
        axis.title.x = element_text(size=9, family = "sans", face = "bold"),
        axis.text.y = element_text(size=8, family = "sans"),
        axis.text.x  = element_text(size=8, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  theme(
    legend.position = "none",
    panel.grid = element_blank())

## offspring versus lv
b <-
  ggplot(data = prey_details_df_nocosts, aes(x = lv, y = offspring)) +
  ggtitle("B") +
  geom_point(size = 0.1, color = "#24262d") +
  xlab("Ballistic length scale (m)") +
  ylab("Number of Offspring") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=9, family = "sans", face = "bold", margin = margin(r = 4)),
        axis.title.x = element_text(size=9, family = "sans", face = "bold"),
        axis.text.y = element_text(size=8, family = "sans"),
        axis.text.x  = element_text(size=8, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  theme(
    legend.position = "none",
    panel.grid = element_blank())

## speed vs lv
c <-  
  ggplot() +
  ggtitle("C") +
  geom_point(data = prey_details_df_nocosts, aes(x = lv, y = speed), size = 0.1, color = "#24262d") +
  labs(x = "Ballistic length scale (m)", y = "Speed (m/s)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=9, family = "sans", face = "bold"),
        axis.title.x = element_text(size=9, family = "sans", face = "bold"),
        axis.text.y = element_text(size=8, family = "sans"),
        axis.text.x  = element_text(size=8, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  theme(
    legend.position = "none",
    panel.grid = element_blank())

# WITH ENERGETIC COSTS (E-H)

## offspring versus speed
d <-
ggplot() +
  ggtitle("D") +
  geom_point(data = prey_details_df_costs, aes(x = speed, y = offspring), size = 0.1, col = "#b02942e5") +
  labs(x = "Speed (m/s)", y = "Number of Offspring") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=9, family = "sans", face = "bold", margin = margin(r = 4)),
        axis.title.x = element_text(size=9, family = "sans", face = "bold"),
        axis.text.y = element_text(size=8, family = "sans"),
        axis.text.x  = element_text(size=8, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  theme(
    legend.position = "none",
    panel.grid = element_blank())

## offspring versus lv
e <-
ggplot() +
  ggtitle("E") +
  geom_point(data = prey_details_df_costs, aes(x = lv, y = offspring), size = 0.1, col = "#b02942e5") +
  xlab("Ballistic length scale (m)") +
  ylab("Number of Offspring") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=9, family = "sans", face = "bold", margin = margin(r = 4)),
        axis.title.x = element_text(size=9, family = "sans", face = "bold"),
        axis.text.y = element_text(size=8, family = "sans"),
        axis.text.x  = element_text(size=8, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  theme(
    legend.position = "none",
    panel.grid = element_blank())

## speed versus lv
f <-  
  ggplot() +
  ggtitle("F") +
  geom_point(data = prey_details_df_costs, aes(x = lv, y = speed), size = 0.1, col = "#b02942e5") +
  labs(x = "Ballistic length scale (m)", y = "Speed (m/s)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=9, family = "sans", face = "bold"),
        axis.title.x = element_text(size=9, family = "sans", face = "bold"),
        axis.text.y = element_text(size=8, family = "sans"),
        axis.text.x  = element_text(size=8, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  theme(
    legend.position = "none",
    panel.grid = element_blank())


# final plot -------------------------------------------------------------------

## combine figures
FIG <- grid.arrange(a,b,c,d,e,f, ncol = 3, nrow = 2)

## save 
ggsave(FIG,
       width = 9, height = 5, units = "in",
       dpi = 600,
       bg = "white",
       file="figures/supplementary/cost_analysis_offspring.png")

# compare evolution of lv with and without energetic costs----------------------

#remove clamped values
prey_details_df_nocosts <- prey_details_df_nocosts %>% 
  filter(abs(lv - 0.7432662) > 1e-6)

prey_details_df_costs <- prey_details_df_costs %>% 
  filter(abs(lv - 0.7475293) > 1e-6)

## calculate the mean for each generation (no energetic costs)
prey_summary_nocosts <- prey_details_df_nocosts %>%
  group_by(generation) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop")

## calculate the mean for each generation (energetic costs)
prey_summary_costs <- prey_details_df_costs %>%
  group_by(generation) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop")

#set seeds
seed <- 1:100

## generation 1 no costs movement track
prey_nocost1 <- prey_summary_nocosts %>% 
  filter(generation == 1)

nocost_tau_v1 <- prey_nocost1$tau_v
nocost_tau_p1 <- prey_nocost1$tau_p
nocost_sig1 <- prey_nocost1$sig

nocost_mod1 <- ctmm(tau = c(nocost_tau_p1, nocost_tau_v1), mu = c(0,0), sigma = nocost_sig1)

nocost_t1 <- sampling(105500)

nocost_t1 <- nocost_t1[nocost_t1 <= 3*nocost_tau_v1]

nocost_df1 <- list()
for(i in seed){
  nocost_tracks1 <- simulate(nocost_mod1, t = nocost_t1)
  
  nocost_df1[[i]] <- as.data.frame(nocost_tracks1)
}

plot_list_nocost <- lapply(nocost_df1, function(telem){

  df <- as.data.frame(telem)
  
  ggplot(df, aes(x = x, y = y)) +
    geom_path(linewidth = 0.5, col = "black") +
    coord_equal() +
    theme_void() +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x  = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_blank(),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent", color = NA),
          plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
    theme(legend.position = "none",
          panel.grid = element_blank())
})

output_dir <- "~/hdrive/GitHub/ballistic-movement/figures/supplementary/generation_moveplots/generation1_nocosts"
for(i in seq_along(plot_list_nocost)) {
  ggsave(
    filename = file.path(output_dir, paste0("nocost_plot_", sprintf("%03d", i), ".png")),
    plot = plot_list_nocost[[i]],
    width = 4 , height = 4, dpi = 300, bg = "white"
  )
}

## generation 1 costs movement track
prey_costs1 <- prey_summary_costs %>% 
  filter(generation == 1)

costs_tau_v1 <- prey_costs1$tau_v
costs_tau_p1 <- prey_costs1$tau_p
costs_sig1 <- prey_costs1$sig

costs_mod1 <- ctmm(tau = c(costs_tau_p1, costs_tau_v1), mu = c(0,0), sigma = costs_sig1)

costs_t1 <- sampling(105500)

cost_df1 <- list()
for(i in seed){
  cost_tracks1 <- simulate(costs_mod1, t = nocost_t1)
  
  cost_df1[[i]] <- as.data.frame(cost_tracks1)
}

plot_list_cost <- lapply(cost_df1, function(telem){
  
  df <- as.data.frame(telem)
  
  ggplot(df, aes(x = x, y = y)) +
    geom_path(linewidth = 0.5, col = "#b02942e5") +
    coord_equal() +
    theme_void() +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x  = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_blank(),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent", color = NA),
          plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
    theme(legend.position = "none",
          panel.grid = element_blank())
})

output_dir <- "~/hdrive/GitHub/ballistic-movement/figures/supplementary/generation_moveplots/generation1_costs"
for(i in seq_along(plot_list_cost)) {
  ggsave(
    filename = file.path(output_dir, paste0("cost_plot_", sprintf("%03d", i), ".png")),
    plot = plot_list_cost[[i]],
    width = 4 , height = 4, dpi = 300, bg = "white"
  )
}


nocost_gen1_tracks <-
  ggplot(nocost_df1[[23]], aes(x = x, y = y)) +
    geom_path(linewidth = 0.5, col = "black") +
    annotation_scale(plot_unit = "m", location = "tr", line_width = 0.5, style = "ticks", text_cex = 0.5) +
    ggplot2::annotate("text", x = min(nocost_df1[[23]]$x), y = min(nocost_df1[[23]]$y), label = "Gen. 1", hjust = 0, vjust = 0, size = 2, family = "sans") +
    coord_equal(xlim = range(nocost_df1[[23]]$x) + c(-0.05, 0.05) * diff(range(nocost_df1[[23]]$x)), ylim = range(nocost_df1[[23]]$y) + c(-0.05, 0.05) * diff(range(nocost_df1[[23]]$y))) +
    theme_void() +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x  = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_blank(),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent", color = NA),
          plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
    theme(legend.position = "none",
          panel.grid = element_blank())

final_nocost1 <- nocost_df1[[23]]
save(final_nocost1, file = "figures/supplementary/generation_moveplots/generation1_nocosts/final_tracks.Rda")

cost_gen1_tracks <-
  ggplot(cost_df1[[15]], aes(x = x, y = y)) +
  geom_path(linewidth = 0.5, col = "#b02942e5") +
  annotation_scale(plot_unit = "m", location = "tr", line_width = 0.5, style = "ticks", text_cex = 0.5) +
  ggplot2::annotate("text", x = min(cost_df1[[15]]$x), y = min(cost_df1[[15]]$y), label = "Gen. 1", hjust = 0, vjust = 0, size = 2, family = "sans") +
  coord_equal(xlim = range(cost_df1[[15]]$x) + c(-0.05, 0.05) * diff(range(cost_df1[[15]]$x)), ylim = range(cost_df1[[15]]$y) + c(-0.05, 0.05) * diff(range(cost_df1[[15]]$y))) +
  theme_void() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  theme(legend.position = "none",
        panel.grid = element_blank())

final_cost1 <- cost_df1[[15]]
save(final_nocost1, file = "figures/supplementary/generation_moveplots/generation1_costs/final_tracks.Rda")

nocost_gen1_tracks <- nocost_gen1_tracks + theme(aspect.ratio = 1)
cost_gen1_tracks <- cost_gen1_tracks + theme(aspect.ratio = 1)

## combine insets
inset_gen1 <- 
  wrap_plots(nocost_gen1_tracks, cost_gen1_tracks, ncol = 2, nrow = 1, heights = c(1, 1), widths = c(1,1))

## generation 1000 no costs movement track
prey_nocost1000 <- prey_summary_nocosts %>% 
  filter(generation == 1000)

nocost_tau_v1000 <- prey_nocost1000$tau_v
nocost_tau_p1000 <- prey_nocost1000$tau_p
nocost_sig1000 <- prey_nocost1000$sig

nocost_mod1000 <- ctmm(tau = c(nocost_tau_p1000, nocost_tau_v1000), mu = c(0,0), sigma = nocost_sig1000)

nocost_df1000 <- list()
for(i in seed){
  nocost_tracks1000 <- simulate(nocost_mod1000, t = nocost_t1)
  
  nocost_df1000[[i]] <- as.data.frame(nocost_tracks1000)
}

plot_list_nocost1000 <- lapply(nocost_df1000, function(telem){
  
  df <- as.data.frame(telem)
  
  ggplot(df, aes(x = x, y = y)) +
    geom_path(linewidth = 0.5, col = "black") +
    coord_equal() +
    theme_void() +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x  = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_blank(),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent", color = NA),
          plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
    theme(legend.position = "none",
          panel.grid = element_blank())
})

output_dir <- "~/hdrive/GitHub/ballistic-movement/figures/supplementary/generation_moveplots/generation1000_nocosts"
for(i in seq_along(plot_list_nocost1000)) {
  ggsave(
    filename = file.path(output_dir, paste0("nocost_plot_", sprintf("%03d", i), ".png")),
    plot = plot_list_nocost1000[[i]],
    width = 4 , height = 4, dpi = 300, bg = "white"
  )
}

## generation 1000 costs movement track
prey_costs1000 <- prey_summary_costs %>% 
  filter(generation == 1)

costs_tau_v1000 <- prey_costs1000$tau_v
costs_tau_p1000 <- prey_costs1000$tau_p
costs_sig1000 <- prey_costs1000$sig

costs_mod1000 <- ctmm(tau = c(costs_tau_p1000, costs_tau_v1000), mu = c(0,0), sigma = costs_sig1000)

cost_df1000 <- list()
for(i in seed){
  cost_tracks1000 <- simulate(costs_mod1000, t = nocost_t1)
  
  cost_df1000[[i]] <- as.data.frame(cost_tracks1000)
}

plot_list_cost1000 <- lapply(cost_df1000, function(telem){
  
  df <- as.data.frame(telem)
  
  ggplot(df, aes(x = x, y = y)) +
    geom_path(linewidth = 0.5, col = "#b02942e5") +
    coord_equal() +
    theme_void() +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x  = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_blank(),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent", color = NA),
          plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
    theme(legend.position = "none",
          panel.grid = element_blank())
})

output_dir <- "~/hdrive/GitHub/ballistic-movement/figures/supplementary/generation_moveplots/generation1000_costs"
for(i in seq_along(plot_list_cost1000)) {
  ggsave(
    filename = file.path(output_dir, paste0("cost_plot_", sprintf("%03d", i), ".png")),
    plot = plot_list_cost1000[[i]],
    width = 4 , height = 4, dpi = 300, bg = "white"
  )
}

nocost_gen1000_tracks <-
  ggplot(nocost_df1000[[3]], aes(x = x, y = y)) +
  geom_path(linewidth = 0.5, col = "black") +
  annotation_scale(plot_unit = "m", location = "tr", line_width = 0.5, style = "ticks", text_cex = 0.5) +
  ggplot2::annotate("text", x = min(nocost_df1000[[3]]$x), y = min(nocost_df1000[[3]]$y), label = "Gen. 1000", hjust = 0, vjust = 0, size = 2, family = "sans") +
  coord_equal(xlim = range(nocost_df1000[[3]]$x) + c(-0.05, 0.05) * diff(range(nocost_df1000[[3]]$x)), ylim = range(nocost_df1000[[3]]$y) + c(-0.05, 0.05) * diff(range(nocost_df1000[[3]]$y))) +
  theme_void() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  theme(legend.position = "none",
        panel.grid = element_blank())

final_nocost1000 <- nocost_df1000[[3]]
save(final_nocost1000, file = "figures/supplementary/generation_moveplots/generation1000_nocosts/final_tracks.Rda")

cost_gen1000_tracks <-
  ggplot(cost_df1000[[25]], aes(x = x, y = y)) +
  geom_path(linewidth = 0.5, col = "#b02942e5") +
  annotation_scale(plot_unit = "m", location = "tr", line_width = 0.5, style = "ticks", text_cex = 0.5) +
  ggplot2::annotate("text", x = min(cost_df1000[[25]]$x), y = min(cost_df1000[[25]]$y), label = "Gen. 1000", hjust = 0, vjust = 0, size = 2, family = "sans") +
  coord_equal(xlim = range(cost_df1000[[25]]$x) + c(-0.05, 0.05) * diff(range(cost_df1000[[25]]$x)), ylim = range(cost_df1000[[25]]$y) + c(-0.05, 0.05) * diff(range(cost_df1000[[25]]$y))) +
  theme_void() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  theme(legend.position = "none",
        panel.grid = element_blank())

final_cost1000 <- cost_df1000[[25]]
save(final_cost1000, file = "figures/supplementary/generation_moveplots/generation1000_costs/final_tracks.Rda")

nocost_gen1000_tracks <- nocost_gen1000_tracks + theme(aspect.ratio = 1)
cost_gen1000_tracks <- cost_gen1000_tracks + theme(aspect.ratio = 1)

## combine insets
inset_gen1000 <- 
  wrap_plots(nocost_gen1000_tracks, cost_gen1000_tracks, ncol = 2, nrow = 1, heights = c(1, 1), widths = c(1,1))

## create figure
p1 <-
  ggplot() +
  geom_point(dat = prey_details_df_costs, aes(x = generation, y = lv), col = "#efbdc7", alpha = 0.1, size = 0.1) +
  geom_point(dat = prey_details_df_nocosts, aes(x = generation, y = lv), col = "#c7c9d1", alpha = 0.1, size = 0.1) +
  geom_line(dat = prey_summary_costs, aes(x = generation, y = lv), col = "#741b2d", linewidth = 0.5) +
  geom_line(dat = prey_summary_nocosts, aes(x = generation, y = lv), col = "#24262d", linewidth = 0.5) +
  labs(y = "Ballistic length scale (m)", x = "Generation") +
  scale_y_continuous(expand = c(0.002,0.002), breaks = c(0, 10, 20, 30, 40, 50, 60)) +
  scale_x_continuous(expand = c(0.002,0.002), breaks = c(0, 250, 500, 750, 1000)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=9, family = "sans", face = "bold"),
        axis.title.x = element_text(size=9, family = "sans", face = "bold"),
        axis.text.y = element_text(size=8, family = "sans"),
        axis.text.x  = element_text(size=8, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  theme(
    legend.position = "none",
    panel.grid = element_blank())

p1 +
  inset_element(inset_gen1, left = 0.005, right = 0.305, bottom = 0.7, top = 1.005) +
  inset_element(inset_gen1000, left = 0.69, right = 0.99, bottom = 0.105, top = 0.6)

# save with movement tracks
ggsave(width = 12, height = 6, units = "in", dpi = 400, file = "figures/maintext/generation_withtracks.png")

## save figure
ggsave(p1, width = 10, height = 5, units = "in", dpi = 900, file = "figures/maintext/105kg_lv_gen.png")


#-------------------------
#No cost model
#-------------------------

#Plots of the trends
plot(offspring ~ sqrt(lv), data = combined_df[combined_df$condition == "no_cost",])
plot(offspring ~ sqrt(speed), data = combined_df[combined_df$condition == "no_cost",])

no_cost_model <- bam(offspring ~
                   s(sqrt(lv), k = 4) +
                   s(sqrt(speed), k = 3), #+
                   #ti(sqrt(lv), sqrt(speed), k = 5),
                 data = combined_df[combined_df$condition == "no_cost",],
                 family = poisson(link = "log"),
                 discrete = T,
                 method = "fREML")

#plot of the model
plot(no_cost_model, pages = 1, scheme = 2)


#-------------------------
#Cost model
#-------------------------

#Plots of the trends
plot(offspring ~ sqrt(lv), data = combined_df[combined_df$condition == "cost",])
plot(offspring ~ sqrt(speed), data = combined_df[combined_df$condition == "cost",])

cost_model <- bam(offspring ~
                   s(sqrt(lv), k = 4) +
                   s(sqrt(speed), k = 3),# +
                   #ti(sqrt(lv), sqrt(speed), k = 5),
                 data = combined_df[combined_df$condition == "cost",],
                 family = poisson(link = "log"),
                 discrete = T,
                 method = "fREML")

#plot of the model
plot(cost_model, pages = 1, scheme = 2)

