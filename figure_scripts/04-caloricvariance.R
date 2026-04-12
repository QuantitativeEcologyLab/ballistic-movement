# script for making plot of ballistic length scale vs caloric distribution in the habitat

#..............................................................................
# load libraries ----

library(tidyverse)
library(patchwork)
library(scico)
library(viridis)
library(mgcv)
library(gridExtra)

source("simulation_scripts/01-prey-functions.R")
#..............................................................................

#load data from folder
calories <- list.files(path = "simulations/prey_results/calorie-variance", 
                    pattern = "prey_details\\.Rds$", 
                    full.names = TRUE) %>% 
  map(~ {
    data_list <- readRDS(.x) %>% 
      bind_rows() %>% 
      # na.omit() %>% 
      group_by(cal_var) %>% 
      filter(generation >= max(generation) - 10) %>% 
      summarise(cal_var = mean(cal_var),
                mean_lv = mean(lv),
                mean_speed = mean(speed))
    return(data_list)
  }) %>%
  list_rbind()

# add column for labeling
calories <- calories %>% 
  filter(cal_var == 0 | cal_var >= 1000) %>%
  mutate(label = case_when(
    cal_var == 0 ~ "label",
    TRUE ~ "other"
  ))

#fit model to data
calories_lv <- glm(mean_lv ~ cal_var,
                data = calories, 
                family = Gamma(link = "log"))

#predict data from model
calories_lv_data <- 
  data.frame(cal_var = seq(min(calories$cal_var)-1000, max(calories$cal_var)+1000, length.out = 100)) %>% 
  mutate(pred = as.data.frame(predict(calories_lv, newdata = ., type = "link", se = TRUE)),
         fit = exp(pred$fit),
         lowerci = exp(pred$fit - pred$se.fit * 1.96),
         upperci = exp(pred$fit + pred$se.fit * 1.96)) %>% 
  select(!pred)

p1 <-
  ggplot() +
  ggtitle("A") +
  geom_ribbon(data = calories_lv_data, aes(ymin = lowerci, ymax = upperci, x = cal_var), alpha = 0.3, fill = "#CC9BA5") +
  geom_line(data = calories_lv_data, aes(x = cal_var, y = fit), col = "#401D1F") +
  geom_point(data = calories, aes(x = cal_var, y = mean_lv, col = label)) +
  # geom_segment(aes(x = 200, y = subset(calories, cal_var == 0)$mean_lv +75, 
  #                  xend = 0.5, yend = subset(calories, cal_var == 0)$mean_lv+10), 
  #              arrow = arrow(length = unit(0.2, "cm")), col = "#0062b8") +
  labs(x = "Caloric Variance", y = expression(bold(l[v] (m)))) +
  scale_color_manual(values = c("label" = "#0062b8", "other" = "grey20")) +
  scale_x_continuous(expand = c(0,0), limits = c(min(calories_lv_data$cal_var), max(calories_lv_data$cal_var))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = NA, linewidth = 1.2),
        axis.title.y = element_text(size=9, family = "sans", face = "bold"),
        axis.title.x = element_text(size=9, family = "sans", face = "bold"),
        axis.text.y = element_text(size=7, family = "sans"),
        axis.text.x  = element_text(size=7, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 10, family = "sans", face = "bold"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
        legend.position = "none")

# ggsave(p1, file = "figures/maintext/patches-vs-lv.png", width = 6, height = 3, units = "in", bg = "white", dpi = 600)

#model speed
calories_speed <- glm(mean_speed ~ cal_var,
                   data = calories, 
                   family = Gamma(link = "log"))

calories_speed_data <- 
  data.frame(cal_var = seq(min(calories$cal_var)-1000, max(calories$cal_var)+1000, length.out = 100)) %>% 
  mutate(pred = as.data.frame(predict(calories_speed, newdata = ., type = "link", se = TRUE)),
         fit = exp(pred$fit),
         lowerci = exp(pred$fit - pred$se.fit * 1.96),
         upperci = exp(pred$fit + pred$se.fit * 1.96)) %>% 
  select(!pred)

p2 <-
  ggplot() +
  ggtitle("B") +
  geom_ribbon(data = calories_speed_data, aes(ymin = lowerci, ymax = upperci, x = cal_var), alpha = 0.3, fill = "#CC9BA5") +
  geom_line(data = calories_speed_data, aes(x = cal_var, y = fit), col = "#401D1F") +
  geom_point(data = calories, aes(x = cal_var, y = mean_speed, col = label)) +
  scale_color_manual(values = c("label" = "#0062b8", "other" = "grey20")) +
  # geom_segment(aes(x = 200, y = subset(calories, cal_var == 0)$mean_speed*1.03, 
  #                  xend = 0.5, yend = subset(calories, cal_var == 0)$mean_speed*1.005), 
  #              arrow = arrow(length = unit(0.2, "cm")), col = "#0062b8") +
  labs(x = "Caloric Variance", y = "Speed (m/s)") +
  scale_x_continuous(expand = c(0, 0), limits = c(min(calories_speed_data$cal_var), max(calories_speed_data$cal_var))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = NA, linewidth = 1.2),
        axis.title.y = element_text(size=9, family = "sans", face = "bold"),
        axis.title.x = element_text(size=9, family = "sans", face = "bold"),
        axis.text.y = element_text(size=7, family = "sans"),
        axis.text.x  = element_text(size=7, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 10, family = "sans", face = "bold"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
        legend.position = "none")

# ggsave(p2, file = "figures/maintext/patches-vs-speed.png", width = 6, height = 3, units = "in", bg = "white", dpi = 600)

#create combined plot
final <- grid.arrange(p1, p2, ncol = 2)

# ggsave(final, file = "presentations/poster-components/figures/calories-analysis.png", width = 9, height = 4, units = "in", bg = "white", dpi = 600)

#habitat caloriesering
base_food <- readRDS("simulations/prey_results/patches/habitats/500_patches.Rds")

mass_prey <- 105500

caloriesvals <- c(0, 2500, 5000, 7500, 10000)

calories_res <- vector("list", length(caloriesvals))
for(idx in seq_along(calories_res)) {
  i <- caloriesvals[idx]
  
  if(i == 1) {
    food <- base_food   # <-- force identical landscape
  } else {
    success <- FALSE
    
    set.seed(123 + idx)
    
    while(!success){
      FOOD <- try({makeHabitat(mass_prey,
                               r = 1,
                               var = i, 
                               mu = 1,
                               target_n = 500,
                               cal = 4000)},
                  silent = TRUE)
      
      success <- !inherits(FOOD, "try-error")
    }
    
    food <- as.data.frame(FOOD)
  }
  
  calories_res[[idx]] <- data.frame(x = food$x,
                                    y = food$y,
                                    cals = food$marks,
                                    cal_var = i)
}


calories_res <- bind_rows(calories_res)

calories_habitats <-
  ggplot(calories_res, aes(x = x, y = y, col = cals)) +
  geom_rect(data = filter(calories_res, cal_var == 0),
            xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf,
            fill = "#ebf6ff",
            inherit.aes = FALSE) +
  geom_point(size = 0.5) +
  facet_wrap(~cal_var, ncol = 5, labeller = as_labeller(function(x) paste("Calorie Variance", x))) +
  scale_color_scico(palette = "berlin",
                    midpoint = 4000
                    )+
  labs(color = "Calories") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = -0.05, size = 10, family = "sans", face = "bold"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 8, family = "sans", face = "bold"),
        legend.position = "left",
        legend.text = element_text(size = 4, family = "sans"),
        legend.title = element_text(size = 6, family = "sans", face = "bold"),
        legend.key.size = unit(0.3, "cm"),
        legend.spacing.y = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent"))
  

# ggsave(habitats, file = "presentations/poster-components/figures/habitat-caloriesering.png", width = 10, height = 1.9, units = "in", bg = "white", dpi = 600)

FIG <- grid.arrange(calories_habitats, final, nrow = 2, heights = c(1,2))

ggsave(FIG, file = "presentations/poster-components/figures/landscape-calories.png", 
       width = 9, height = 5, units = "in", dpi = 600, bg = "white")
