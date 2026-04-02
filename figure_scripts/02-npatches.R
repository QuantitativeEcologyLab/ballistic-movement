# script for making plots for simulations testing the number of patches in landscape

#..............................................................................
# load libraries ----

library(tidyverse)
library(patchwork)
library(scico)
library(viridis)

#..............................................................................


cmbnd <- list.files(path = "simulations/prey_results/patches", 
                    pattern = "prey_details\\.Rds$", 
                    full.names = TRUE)  %>% 
  map(~{
    data_list <- readRDS(.x) %>% 
      bind_rows() %>% 
      na.omit() %>% 
      group_by(tot_patch) %>% 
      filter(generation >= max(generation) - 10) %>% 
      summarise(tot_patch = mean(tot_patch),
                tot_cal = mean(tot_patch * cal_per_patch),
                mean_lv = mean(lv),
                mean_speed = mean(speed))
    return(data_list)
  }) %>% 
  list_rbind()

#model lv
mod_lv <- glm(mean_lv ~ tot_cal,
              data = cmbnd, 
              family = Gamma(link = "log"))

New_Data_lv <- data.frame(tot_cal = seq(min(cmbnd$tot_cal) - 90000, max(cmbnd$tot_cal) + 90000, length.out = 100))

#generate predictions from GLM
preds_lv <- predict(mod_lv, newdata = New_Data_lv, type = "link", se = TRUE)

New_Data_lv$fit <- exp(preds_lv$fit)
New_Data_lv$lowerci <- exp(preds_lv$fit - preds_lv$se.fit * 1.96)
New_Data_lv$upperci <- exp(preds_lv$fit + preds_lv$se.fit * 1.96)

p1 <-
  ggplot() +
  ggtitle("A") +
  geom_point(data = cmbnd, aes(x = tot_cal, y = mean_lv)) +
  geom_ribbon(data = New_Data_lv, aes(ymin = lowerci, ymax = upperci, x = tot_cal), alpha = 0.3, fill = "#CC9BA5") +
  geom_line(data = New_Data_lv, aes(x = tot_cal, y = fit), col = "#401D1F") +
  labs(x = "Number of Calories in Landscape", y = expression(bold(l[v]))) +
  scale_x_continuous(expand = c(0,0), limits = c(min(New_Data_lv$tot_cal), max(New_Data_lv$tot_cal))) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
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
        legend.text = element_text(size = 6, family = "sans"),
        legend.title = element_text(size = 8, family = "sans", face = "bold"),
        # legend.key.size = unit(0.2, "cm"),
        # legend.spacing.y = unit(0.1, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent"))

ggsave(p1, file = "figures/maintext/patches-vs-lv.png", width = 6, height = 3, units = "in", bg = "white", dpi = 600)

#model speed
mod_speed <- glm(mean_speed ~ tot_cal,
                 data = cmbnd, 
                 family = Gamma(link = "log"))

New_Data_speed <- data.frame(tot_cal = seq(min(cmbnd$tot_cal) - 90000, max(cmbnd$tot_cal) + 90000, length.out = 100))

#generate predictions from GLM
preds_speed <- predict(mod_speed, newdata = New_Data_speed, type = "link", se = TRUE)

New_Data_speed$fit <- exp(preds_speed$fit)
New_Data_speed$lowerci <- exp(preds_speed$fit - preds_speed$se.fit * 1.96)
New_Data_speed$upperci <- exp(preds_speed$fit + preds_speed$se.fit * 1.96)

p2 <-
  ggplot() +
  ggtitle("B") +
  geom_point(data = cmbnd, aes(x = tot_cal, y = mean_speed)) +
  geom_ribbon(data = New_Data_speed, aes(ymin = lowerci, ymax = upperci, x = tot_cal), alpha = 0.3, fill = "#CC9BA5") +
  geom_line(data = New_Data_speed, aes(x = tot_cal, y = fit), col = "#401D1F") +
  labs(x = "Number of Calories in Landscape", y = "Speed (m/s)") +
  scale_x_continuous(expand = c(0, 0), limits = c(min(New_Data_speed$tot_cal), max(New_Data_speed$tot_cal))) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
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
        legend.text = element_text(size = 6, family = "sans"),
        legend.title = element_text(size = 8, family = "sans", face = "bold"),
        # legend.key.size = unit(0.2, "cm"),
        # legend.spacing.y = unit(0.1, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent"))

ggsave(p2, file = "figures/maintext/patches-vs-speed.png", width = 6, height = 3, units = "in", bg = "white", dpi = 600)

#combine plots
final <- grid.arrange(p1, p2, ncol = 2, top = textGrob("Landscape Density", gp = gpar(fontsize = 16, family = "sans", face = "bold")))

ggsave(final, file = "presentations/poster-components/figures/density-analysis.png", width = 9, height = 4, units = "in", bg = "white", dpi = 600)
