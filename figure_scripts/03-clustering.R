# script for making plot of ballistic length scale vs degree of patch clustering in the habitat

#..............................................................................
# load libraries ----

library(tidyverse)
library(patchwork)
library(scico)
library(viridis)

#..............................................................................


clust <- list.files(path = "simulations/prey_results/clustering", 
                    pattern = "prey_details\\.Rds$", 
                    full.names = TRUE) %>% 
  map(~ {
    id_val <- as.numeric(str_extract(basename(.x), "\\d+"))
    
    data_list <- readRDS(.x) %>% map_dfr(., ~ mutate(.x, mu = id_val))
    
    data_list <- data_list %>%
      group_by(mu) %>%
      filter(generation >= max(generation) - 10) %>% 
      summarise(tot_patch = mean(tot_patch),
                tot_cal = mean(tot_patch * cal_per_patch),
                mu = mean(as.numeric(mu)),
                mean_lv = mean(lv),
                mean_speed = mean(speed))
    
    return(data_list)
  }) %>%
  list_rbind()

mod_lv <- glm(mean_lv ~ mu,
              data = clust, 
              family = Gamma(link = "log"))

New_Data_lv <- data.frame(mu = seq(min(clust$mu) - 2, max(clust$mu) + 2, length.out = 100))

#generate predictions from GLM
preds_lv <- predict(mod_lv, newdata = New_Data_lv, type = "link", se = TRUE)

New_Data_lv$fit <- exp(preds_lv$fit)
New_Data_lv$lowerci <- exp(preds_lv$fit - preds_lv$se.fit * 1.96)
New_Data_lv$upperci <- exp(preds_lv$fit + preds_lv$se.fit * 1.96)

p1 <-
  ggplot() +
  ggtitle("A") +
  geom_point(data = clust, aes(x = mu, y = mean_lv)) +
  geom_ribbon(data = New_Data_lv, aes(ymin = lowerci, ymax = upperci, x = mu), alpha = 0.3, fill = "#CC9BA5") +
  geom_line(data = New_Data_lv, aes(x = mu, y = fit), col = "#401D1F") +
  labs(x = "Degree of Clustering", y = expression(bold(l[v]))) +
  scale_x_continuous(expand = c(0,0), limits = c(min(New_Data_lv$mu), max(New_Data_lv$mu))) +
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
        legend.margin = margin(0,0,0,0),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent"))

ggsave(p1, file = "figures/maintext/patches-vs-lv.png", width = 6, height = 3, units = "in", bg = "white", dpi = 600)

#model speed
mod_speed <- glm(mean_speed ~ mu,
                 data = clust, 
                 family = Gamma(link = "log"))

New_Data_speed <- data.frame(mu = seq(min(clust$mu) - 2, max(clust$mu) + 2, length.out = 100))

#generate predictions from GLM
preds_speed <- predict(mod_speed, newdata = New_Data_speed, type = "link", se = TRUE)

New_Data_speed$fit <- exp(preds_speed$fit)
New_Data_speed$lowerci <- exp(preds_speed$fit - preds_speed$se.fit * 1.96)
New_Data_speed$upperci <- exp(preds_speed$fit + preds_speed$se.fit * 1.96)

p2 <-
  ggplot() +
  ggtitle("B") +
  geom_point(data = clust, aes(x = mu, y = mean_speed)) +
  geom_ribbon(data = New_Data_speed, aes(ymin = lowerci, ymax = upperci, x = mu), alpha = 0.3, fill = "#CC9BA5") +
  geom_line(data = New_Data_speed, aes(x = mu, y = fit), col = "#401D1F") +
  labs(x = "Degree of Clustering", y = "Speed (m/s)") +
  scale_x_continuous(expand = c(0, 0), limits = c(min(New_Data_speed$mu), max(New_Data_speed$mu))) +
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
        legend.margin = margin(0,0,0,0),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent"))

ggsave(p2, file = "figures/maintext/patches-vs-speed.png", width = 6, height = 3, units = "in", bg = "white", dpi = 600)

#create combined plot
final <- grid.arrange(p1, p2, ncol = 2, top = textGrob("Landscape Clustering", gp = gpar(fontsize = 16, family = "sans", face = "bold")))

ggsave(final, file = "presentations/poster-components/figures/clust-analysis.png", width = 9, height = 4, units = "in", bg = "white", dpi = 600)
