# Preamble ----

# Set the working directory
setwd("~/hdrive/GitHub/ballistic-movement")

library(ggplot2)
library(tidyverse)
library(dplyr)
library(mgcv)
library(scico)
library(gridExtra)

source("scripts/functions.R")

# without costs associated with movement ---------------------------------------
load('simulations/supplementary/105500g_nocost_prey_details.Rda')

prey_details_df_nocosts <- do.call(rbind, prey_details)

# with costs associated with movement ------------------------------------------

load('simulations/prey_results/105500g_prey_details.Rda')

prey_details_df_costs <- do.call(rbind, prey_details)

prey_details_df_costs <- subset(prey_details_df_costs, select = - costs)

# bam model --------------------------------------------------------------------

prey_details_df_nocosts$condition <- "no_cost"
prey_details_df_costs$condition <- "cost"
combined_df <- rbind(prey_details_df_nocosts, prey_details_df_costs)
combined_df$condition <- factor(combined_df$condition)

bam_model <- bam(cal_net ~ condition +
                   s(lv, by = condition, k = 10) +
                   s(speed, by = condition, k = 10) +
                   ti(lv, speed, by = condition, k = 5),
                 data = combined_df,
                 method = "fREML",
                 discrete = T,
                 nthreads = 5,
                 control = list(maxit = 800))

grid <- expand.grid(
  lv = seq(min(combined_df$lv), max(combined_df$lv), length.out = 100),
  speed= seq(min(combined_df$speed), max(combined_df$speed), length.out = 100),
  condition = factor(c("cost", "no_cost"), levels = c("no_cost", "cost"))
)

grid$cal_net_pred <- predict(bam_model, newdata = grid)

grid_cost <- subset(grid, condition == "cost")
grid_nocost <- subset(grid, condition == "no_cost")

ggplot(grid, aes(x = lv, y = speed, fill = cal_net_pred)) +
  geom_tile() +
  facet_wrap(~condition) + 
  scale_fill_scico(palette = "bamO", midpoint = 0) + 
  theme_bw()

# create figures ---------------------------------------------------------------

a <-
  ggplot() +
  ggtitle("A") +
  geom_point(data = prey_details_df_nocosts, aes(x = speed, y = cal_net), size = 0.1, color = "#24262d") +
  labs(x = "Speed (m/s)", y = "Net calories") +
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

b <-
  ggplot(data = prey_details_df_nocosts, aes(x = lv, y = cal_net)) +
  ggtitle("B") +
  geom_point(size = 0.1, color = "#24262d") +
  xlab(expression(bold(l[v])))+
  ylab("Net calories") +
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

c <-  
  ggplot() +
  ggtitle("C") +
  geom_point(data = prey_details_df_nocosts, aes(x = lv, y = speed), size = 0.1, color = "#24262d") +
  labs(x = expression(bold(l[v])), y = "Speed (m/s)") +
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

# observe the range and histogram of the prediction data
range(grid_nocost$cal_net_pred)
hist(grid_nocost$cal_net_pred) #heavily skewed, long tail

# due to the long tail, it may be helpful to cap the predictions with a quantile to limit the range and make the difference more noticeable
# set predictions to only include the 90 quantile
grid_nocost$cal_net_pred2 <-grid_nocost$cal_net_pred
grid_nocost[grid_nocost$cal_net_pred2 > quantile(grid_nocost$cal_net_pred2 , 0.95),"cal_net_pred2"] <- quantile(grid_nocost$cal_net_pred2 , 0.95)

# customizing the breaks
nocost_min_z <- min(grid_nocost$cal_net_pred2)
nocost_max_z <- max(grid_nocost$cal_net_pred2)

nocost_step <- (nocost_max_z - nocost_min_z) / 25 
nocost_breaks <- seq(nocost_min_z, nocost_max_z, nocost_step)

nocost_breaks <- nocost_breaks[nocost_breaks <= 6000000]

d <-  
  ggplot() +
  ggtitle("D") +
  geom_raster(data = grid_nocost, aes(x = lv, y = speed, fill = cal_net_pred2)) +
  #geom_point(data = prey_details_df_nocosts, aes(x = lv, y = speed, color = generation)) +
  geom_contour(data = grid_cost, aes(x = lv, y = speed, z = cal_net_pred2), color = 'black', linewidth = 0.1,
               na.rm = TRUE, breaks = nocost_breaks) +
  labs(x = expression(bold(l[v])), y = "Speed (m/s)", fill = "predicted\nnet calories") +
  scale_fill_scico(palette = 'vikO', midpoint = 0) +
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
    legend.text = element_text(size = 6, family = "sans"),
    legend.title = element_text(size = 6, family = "sans", face = "bold"),
    legend.key.size = unit(0.3, "cm"),
    legend.spacing.y = unit(0.01, "cm"),
    legend.margin = margin(0,0,0,0),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent"))

e <-
ggplot() +
  ggtitle("E") +
  geom_point(data = prey_details_df_costs, aes(x = speed, y = cal_net), size = 0.1, col = "#741b2d") +
  labs(x = "Speed (m/s)", y = "Net calories") +
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

f <-
ggplot() +
  ggtitle("F") +
  geom_point(data = prey_details_df_costs, aes(x = lv, y = cal_net), size = 0.1, col = "#741b2d") +
  xlab(expression(bold(l[v])))+
  ylab("Net calories") +
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

g <-  
  ggplot() +
  ggtitle("G") +
  geom_point(data = prey_details_df_costs, aes(x = lv, y = speed), size = 0.1, col = "#741b2d") +
  labs(x = expression(bold(l[v])), y = "Speed (m/s)") +
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

# observe the range and histogram of the prediction data
range(grid_cost$cal_net_pred)
hist(grid_cost$cal_net_pred) #heavily skewed, long tail

# due to the long tail, it may be helpful to cap the predictions with a quantile to limit the range and make the difference more noticeable
# set predictions to only include the 90 quantile
grid_cost$cal_net_pred2 <-grid_cost$cal_net_pred
grid_cost[grid_cost$cal_net_pred2 > quantile(grid_cost$cal_net_pred2 , 0.95),"cal_net_pred2"] <- quantile(grid_cost$cal_net_pred2 , 0.95)

# shows where the true negative patches are: 
# grid_cost$cal_net_pred3 <-grid_cost$cal_net_pred
# grid_cost[grid_cost$cal_net_pred3 > 0,"cal_net_pred3"] <- 1
# grid_cost[grid_cost$cal_net_pred3 < 0,"cal_net_pred3"] <- -1

# customizing the breaks
min_z <- min(grid_cost$cal_net_pred2)
max_z <- max(grid_cost$cal_net_pred2)

step <- (max_z - min_z) / 25 
breaks <- seq(min_z, max_z, step)

breaks <- breaks[breaks <= 6000000]

h <-  
ggplot() +
  ggtitle("H") +
  geom_raster(data = grid_cost, aes(x = lv, y = speed, fill = cal_net_pred2)) +
  #geom_point(data = prey_details_df_costs, aes(x = lv, y = speed, color = generation)) +
  geom_contour(data = grid_cost, aes(x = lv, y = speed, z = cal_net_pred2), color = 'black', linewidth = 0.1,
               na.rm = TRUE, breaks = breaks) +
  labs(x = expression(bold(l[v])), y = "Speed (m/s)", fill = "predicted\nnet calories") +
  scale_fill_scico(palette = 'vikO', midpoint = 0) +
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
    legend.text = element_text(size = 6, family = "sans"),
    legend.title = element_text(size = 6, family = "sans", face = "bold"),
    legend.key.size = unit(0.3, "cm"),
    legend.spacing.y = unit(0.01, "cm"),
    legend.margin = margin(0,0,0,0),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent"))

# final plot -------------------------------------------------------------------

FIG <- grid.arrange(a,b,c,d,e,f,g,h, ncol = 4, nrow = 2)

ggsave(FIG,
       width = 11, height = 5, units = "in",
       dpi = 600,
       bg = "white",
       file="figures/supplementary/cost_analysis.png")

# compare evolution of lv

prey_summary_nocosts <- prey_details_df_nocosts %>%
  group_by(generation) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop")

prey_summary_costs <- prey_details_df_costs %>%
  group_by(generation) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop")

p1 <-
  ggplot() +
  geom_point(dat = prey_details_df_costs, aes(x = generation, y = lv), col = "#efbdc7", alpha = 0.1, size = 0.1) +
  geom_point(dat = prey_details_df_nocosts, aes(x = generation, y = lv), col = "#c7c9d1", alpha = 0.1, size = 0.1) +
  geom_line(dat = prey_summary_costs, aes(x = generation, y = lv), col = "#741b2d", linewidth = 0.8) +
  geom_line(dat = prey_summary_nocosts, aes(x = generation, y = lv), col = "#24262d", linewidth = 0.8) +
  labs(y = expression(bold(l[v])), x = "Generation") +
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

ggsave(p1, width = 10, height = 5, units = "in", dpi = 900, file = "figures/maintext/105kg_lv_gen.png")



