library(ggplot2)
library(tidyverse)
library(dplyr)
library(mgcv)
library(scico)
library(gridExtra)

setwd("~/H/GitHub/ballistic-movement/sim_results/supporting_analysis")

# without costs associated with movement ---------------------------------------
load('June16_30000g_withoutcosts_prey_details.Rda')

prey_details_df_nocosts <- do.call(rbind, prey_details)

# ggplot(prey_details_df, aes(x = lv, y = speed, color = cal_net)) +
#   geom_point() +
#   labs(x = expression("ballistic length scale " * l[v]), y = "speed (m/s)") +
#   scale_color_scico(palette = 'romaO', midpoint = 0) +
#   theme_minimal()
# 
# ggplot(prey_details_df, aes(x = lv, y = cal_net)) +
#   geom_point() +
#   theme_minimal()
# 
# ggplot(prey_details_df, aes(x = generation, y = speed)) +
#   stat_summary(fun = mean, geom = "line", linewidth = 1) +
#   theme_minimal()
# 
# ggplot(prey_details_df, aes(x = speed, y = cal_net)) +
#   geom_point() +
#   theme_minimal()

# GAM model
gam_model_nocosts <- gam(cal_net ~ s(lv, k = 5) +
                         s(speed, k = 5) +
                         ti(lv, speed, k = 5),
                       family = gaussian(),
                       data = prey_details_df,
                       method = "REML")

summary(gam_model_nocosts)

# predict values with GAM model
lv_seq <- seq(min(prey_details_df$lv), max(prey_details_df$lv), length.out = 100)
speed_seq <- seq(min(prey_details_df$speed), max(prey_details_df$speed), length.out = 100)

pred_nocosts <- expand.grid(lv = lv_seq, speed = speed_seq)

pred_nocosts$cal_net_pred <- predict(gam_model_nocosts, newdata = pred_nocosts)

summary(pred_nocosts$cal_net_pred)

# plot the predictions
plot(predict(gam_model_nocosts), prey_details_df$cal_net)
abline(0, 1, col = "red")

# with costs associated with movement ------------------------------------------

load('June16_30000g_withcosts_prey_details.Rda')

# make data sets compatible
prey_details_df_costs <- do.call(rbind, prey_details)

# # cal_net ~ speed
# ggplot(prey_details_df_costs, aes(x = speed, y = cal_net)) +
#   geom_point() +
#   theme_minimal()
# 
# # cal_net ~ lv
# ggplot(prey_details_df_costs, aes(x = lv, y = cal_net)) +
#   geom_point() +
#   theme_minimal()
# 
# # speed ~ lv
# ggplot(prey_details_df_costs, aes(x = lv, y = speed)) +
#   geom_point() +
#   theme_minimal()

# GAM model
gam_model_costs <- gam(cal_net ~ s(lv, k = 5) +
                         s(speed, k = 5) +
                         ti(lv, speed, k = 5),
                       family = gaussian(),
                       data = prey_details_df_costs,
                       method = "REML")

summary(gam_model_final)

# predict values with GAM model
lv_seq <- seq(min(prey_details_df_costs$lv), max(prey_details_df_costs$lv), length.out = 100)
speed_seq <- seq(min(prey_details_df_costs$speed), max(prey_details_df_costs$speed), length.out = 100)

pred_costs <- expand.grid(lv = lv_seq, speed = speed_seq)

pred_costs$cal_net_pred <- predict(gam_model_costs, newdata = pred_costs)

summary(pred_costs$cal_net_pred)

plot(predict(gam_model_costs), prey_details_df_costs$cal_net)
abline(0, 1, col = "red")

cor(pred_costs$cal_net_pred, pred_nocosts$cal_net_pred)

# bam model --------------------------------------------------------------------

prey_details_df_nocosts$condition <- "no_cost"
prey_details_df_costs$condition <- "cost"
combined_df <- rbind(prey_details_df_nocosts, prey_details_df_costs)
combined_df$condition <- factor(combined_df$condition)

bam_model <- bam(cal_net ~ s(lv, by = condition, k = 5) +
                   s(speed, by = condition, k = 5) +
                   ti(lv, speed, by = condition, k = 5),
                 data = combined_df,
                 method = "fREML",
                 discrete = T,
                 nthreads = 5)

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
  scale_fill_scico(palette = "romaO", midpoint = 0) + 
  theme_bw()

# create figures ---------------------------------------------------------------

vmin <- min(c(grid_nocost$cal_net_pred, grid_cost$cal_net_pred))
vmax <- max(c(grid_nocost$cal_net_pred, grid_cost$cal_net_pred))

absmax <- max(abs(vmin), abs(vmax))

a <-
  ggplot() +
  ggtitle("A") +
  geom_point(data = prey_details_df, aes(x = speed, y = cal_net), size = 0.1) +
  labs(x = "speed (m/s)", y = "net calories") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=6, family = "sans", face = "bold"),
        axis.title.x = element_text(size=6, family = "sans", face = "bold"),
        axis.text.y = element_text(size=4, family = "sans"),
        axis.text.x  = element_text(size=4, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm"))
print(a)

b <-
  ggplot(data = prey_details_df, aes(x = lv, y = cal_net)) +
  ggtitle("B") +
  geom_point(size = 0.1) +
  xlab(expression(bold(l[v])))+
  ylab("net calories") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=6, family = "sans", face = "bold"),
        axis.title.x = element_text(size=6, family = "sans", face = "bold"),
        axis.text.y = element_text(size=4, family = "sans"),
        axis.text.x  = element_text(size=4, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) 
print(b)

c <-  
  ggplot() +
  ggtitle("C") +
  geom_tile(data = grid_nocost, aes(x = lv, y = speed, fill = cal_net_pred)) +
  labs(x = expression(bold(l[v])), y = "speed (m/s)", fill = "predicted\nnet calories") +
  scale_fill_scico(palette = 'romaO', midpoint = 0, limits = c(-absmax, absmax)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=6, family = "sans", face = "bold"),
        axis.title.x = element_text(size=6, family = "sans", face = "bold"),
        axis.text.y = element_text(size=4, family = "sans"),
        axis.text.x  = element_text(size=4, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        legend.position = "right",
        legend.text = element_text(size = 4, family = "sans"),
        legend.title = element_text(size = 5, family = "sans", face = "bold"),
        legend.key.size = unit(0.2, "cm"),
        legend.spacing.y = unit(0.1, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm"))
print(c)

e <-
ggplot() +
  ggtitle("E") +
  geom_point(data = prey_details_df_costs, aes(x = speed, y = cal_net), size = 0.1) +
  labs(x = "speed (m/s)", y = "net calories") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=6, family = "sans", face = "bold"),
        axis.title.x = element_text(size=6, family = "sans", face = "bold"),
        axis.text.y = element_text(size=4, family = "sans"),
        axis.text.x  = element_text(size=4, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm"))
print(e)

f <-
ggplot() +
  ggtitle("F") +
  geom_point(data = prey_details_df_costs, aes(x = lv, y = cal_net), size = 0.1) +
  xlab(expression(bold(l[v])))+
  ylab("net calories") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=6, family = "sans", face = "bold"),
        axis.title.x = element_text(size=6, family = "sans", face = "bold"),
        axis.text.y = element_text(size=4, family = "sans"),
        axis.text.x  = element_text(size=4, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) 
print(f)

g <-  
ggplot() +
  ggtitle("G") +
  geom_tile(data = grid_cost, aes(x = lv, y = speed, fill = cal_net_pred)) +
  labs(x = expression(bold(l[v])), y = "speed (m/s)", fill = "predicted\nnet calories") +
  scale_fill_scico(palette = 'romaO', midpoint = 0, limits= c(-absmax, absmax)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=6, family = "sans", face = "bold"),
        axis.title.x = element_text(size=6, family = "sans", face = "bold"),
        axis.text.y = element_text(size=4, family = "sans"),
        axis.text.x  = element_text(size=4, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        legend.position = "right",
        legend.text = element_text(size = 4, family = "sans"),
        legend.title = element_text(size = 5, family = "sans", face = "bold"),
        legend.key.size = unit(0.2, "cm"),
        legend.spacing.y = unit(0.1, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm"))
print(g)

# final plot -------------------------------------------------------------------

FIG <- grid.arrange(a,b,c,e,f,g, ncol = 3, nrow = 2)

ggsave(FIG,
       width = 6.86, height = 3.5, units = "in",
       dpi = 600,
       bg = "white",
       file="~/H/GitHub/ballistic-movement/figures/supporting_analysis/costsvsnocosts_modelresults.png")

# scrap model fitting ----------------------------------------------------------
gam <- gam(cal_net ~ ti(lv, speed), family = tw(), data = prey_details_df, method = "REML")
gam_model <- gam(cal_net ~ ti(lv, speed), family = gaussian(), data = prey_details_df, method = "REML")

par(mfrow = c(2,2))
qqnorm(residuals(gam_model))
qqline(residuals(gam_model))

hist(residuals(gam_model), breaks = 50)

plot(fitted(gam_model), residuals(gam_model))

plot(residuals(gam_model), type = "l")

plot(gam_model)

summary(gam_model)

plot(gam_model, pages = 1, residuals = TRUE)

par(mfrow = c(2,2))
gam.check(gam_model)

lv_seq <- seq(min(prey_details_df$lv), max(prey_details_df$lv), length.out = 200)
speed_seq <- seq(min(prey_details_df$speed), max(prey_details_df$speed), length.out = 200)

prediction_grid <- expand.grid(lv = lv_seq, speed = speed_seq)

prediction_grid$cal_net_pred <- predict(gam_model, newdata = prediction_grid)

summary(prediction_grid$cal_net_pred)

ggplot(prediction_grid, aes(x = lv, y = speed, fill = cal_net_pred)) +
  geom_tile() +
  labs(x = expression("ballistic length scale " * l[v]), y = "speed (m/s)") +
  scale_fill_viridis_c() +
  theme_minimal()

ggplot(prey_details_df, aes(x = lv, y = speed, z = cal_net)) +
  stat_summary_hex(
    fun = mean,
    bins = 50) +
  labs(x = expression("ballistic length scale " * l[v]), y = "speed (m/s)") +
  scale_fill_scico(palette = 'romaO', midpoint = 0) +
  theme_minimal()

gam_model2 <- gam(
  list(
    cal_net ~ s(lv, k = 5) + s(speed, k = 5) + ti(lv, speed, k = 5),
    ~ s(lv, k = 5) + s(speed, k = 5) + ti(lv, speed, k = 5)),
  family = gaulss(),
  data = prey_details_df,
  method = 'REML')

summary(gam_model2)

combined_df <- pred_costs
combined_df$cal_net_pred_nocosts <- pred_nocosts$cal_net_pred
combined_df$cal_net_diff <- 
  combined_df$cal_net_pred - combined_df$cal_net_pred_nocosts

ggplot(combined_df, aes(x = lv, y = speed, fill = cal_net_diff)) +
  geom_tile() +
  scale_fill_scico(palette = "vik", midpoint = 0)

