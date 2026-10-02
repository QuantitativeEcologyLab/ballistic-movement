# script for making plot of ballistic length scale vs caloric distribution in the habitat

#..............................................................................
# load libraries ----

library(tidyverse)
library(patchwork)
library(scico)
library(viridis)
library(gridExtra)

source("simulation_scripts/01-prey-functions.R")

#..............................................................................

#load data from folder ----

## original simulation files
calories <- list.files(path = "simulations/prey_results/calorie-variance/originals/prey_details", 
                    pattern = "prey_details\\.Rds$", 
                    full.names = TRUE) %>% 
  map(~ {
    data_list <- readRDS(.x) %>% 
      bind_rows() %>% 
      # na.omit() %>% 
      group_by(cal_var) %>% 
      filter(generation >= max(generation) - 10) %>% 
      summarise(cal = mean(cal_per_patch),
                cv = mean(cal_var),
                cov = mean(cov),
                mean_lv = mean(lv),
                var_lv = var(lv),
                mean_speed = mean(speed),
                var_speed = mean(speed),
                patches = mean(patches),
                mean_offspring = mean(offspring),
                offspring_var = var(offspring))
    return(data_list)
  }) %>%
  list_rbind()

## first round of duplicates
dups <- list.files(path = "simulations/prey_results/calorie-variance/duplicates", 
                       pattern = "prey_details\\.Rds$", 
                       full.names = TRUE) %>% 
  map(~ {
    data_list <- readRDS(.x) %>% 
      bind_rows() %>% 
      # na.omit() %>% 
      group_by(cal_var) %>% 
      filter(generation >= max(generation) - 10) %>% 
      summarise(cal = mean(cal_per_patch),
                cv = mean(cal_var),
                mean_lv = mean(lv),
                var_lv = var(lv),
                mean_speed = mean(speed),
                var_speed = mean(speed),
                patches = mean(patches),
                mean_offspring = mean(offspring),
                offspring_var = var(offspring))
    return(data_list)
  }) %>%
  list_rbind()
  
## combine datasets to plot
## add column for labelling repeated simulation (500 patches, no clustering, no variance)
calories <- calories %>% 
  mutate(
    label = case_when(
      cv == 0 ~ "label",
      TRUE ~ "other"
    ))

# modelling ----

## fit lv ~ variance model
calories_lv <- glm(mean_lv ~ cov,
                   data = calories, 
                   family = Gamma(link = "log"))

## predict data for lv
calories_lv_data <- 
  data.frame(cov = seq(min(calories$cov)-0.1, 
                      max(calories$cov)+0.1, 
                      length.out = 100)) %>% 
  mutate(pred = as.data.frame(predict(calories_lv, newdata = ., type = "link", se = TRUE)),
         fit = exp(pred$fit),
         lowerci = exp(pred$fit - pred$se.fit * 1.96),
         upperci = exp(pred$fit + pred$se.fit * 1.96)) %>% 
  select(!pred)

## fit lv ~ variance model for duplicated simulations
dups_lv <- glm(mean_lv ~ cv,
                   data = dups, 
                   family = Gamma(link = "log"))

## predict data for lv
dups_lv_data <- 
  data.frame(cv = seq(min(calories$cv)-0.1, 
                      max(calories$cv)+0.1, 
                      length.out = 100)) %>% 
  mutate(pred = as.data.frame(predict(dups_lv, newdata = ., type = "link", se = TRUE)),
         fit = exp(pred$fit),
         lowerci = exp(pred$fit - pred$se.fit * 1.96),
         upperci = exp(pred$fit + pred$se.fit * 1.96)) %>% 
  select(!pred)

## fit lv ~ variance model for all simulations
comb_lv <- glm(mean_lv ~ cv,
                   data = comb, 
                   family = Gamma(link = "log"))

## predict data for lv
comb_lv_data <- 
  data.frame(cv = seq(min(calories$cv)-0.1, 
                      max(calories$cv)+0.1, 
                      length.out = 100)) %>% 
  mutate(pred = as.data.frame(predict(comb_lv, newdata = ., type = "link", se = TRUE)),
         fit = exp(pred$fit),
         lowerci = exp(pred$fit - pred$se.fit * 1.96),
         upperci = exp(pred$fit + pred$se.fit * 1.96)) %>% 
  select(!pred)

## fit speed ~ variance model
dups_speed <- glm(mean_speed ~ cv,
                      data = calories, 
                      family = Gamma(link = "log"))

## predict data for speed
dups_speed_data <- 
  data.frame(cv = seq(min(calories$cv)-0.1, max(calories$cv)+0.1, length.out = 100)) %>% 
  mutate(pred = as.data.frame(predict(dups_speed, newdata = ., type = "link", se = TRUE)),
         fit = exp(pred$fit),
         lowerci = exp(pred$fit - pred$se.fit * 1.96),
         upperci = exp(pred$fit + pred$se.fit * 1.96)) %>% 
  select(!pred)

## fit speed ~ variance model
comb_speed <- glm(mean_speed ~ cv,
                      data = comb, 
                      family = Gamma(link = "log"))

## predict data for speed
comb_speed_data <- 
  data.frame(cv = seq(min(calories$cv)-0.1, max(calories$cv)+0.1, length.out = 100)) %>% 
  mutate(pred = as.data.frame(predict(comb_speed, newdata = ., type = "link", se = TRUE)),
         fit = exp(pred$fit),
         lowerci = exp(pred$fit - pred$se.fit * 1.96),
         upperci = exp(pred$fit + pred$se.fit * 1.96)) %>% 
  select(!pred)



## fit speed ~ variance model
calories_speed <- glm(mean_speed ~ cv,
                      data = calories, 
                      family = Gamma(link = "log"))

## predict data for speed
calories_speed_data <- 
  data.frame(cv = seq(min(calories$cv)-0.1, max(calories$cv)+0.1, length.out = 100)) %>% 
  mutate(pred = as.data.frame(predict(calories_speed, newdata = ., type = "link", se = TRUE)),
         fit = exp(pred$fit),
         lowerci = exp(pred$fit - pred$se.fit * 1.96),
         upperci = exp(pred$fit + pred$se.fit * 1.96)) %>% 
  select(!pred)

## fit offspring ~ variance model
calories_offspring <- glm(mean_offspring ~ cv,
                          data = calories, 
                          family = Gamma(link = "log"))

## predict data for offspring
calories_offspring_data <- 
  data.frame(cv = seq(min(calories$cv)-0.1, max(calories$cv)+0.1, length.out = 100)) %>% 
  mutate(pred = as.data.frame(predict(calories_offspring, newdata = ., type = "link", se = TRUE)),
         fit = exp(pred$fit),
         lowerci = exp(pred$fit - pred$se.fit * 1.96),
         upperci = exp(pred$fit + pred$se.fit * 1.96)) %>% 
  select(!pred)

## fit offspring variance ~ caloric variance model
calories_offspring_var <- glm(offspring_var ~ cv,
                          data = calories, 
                          family = Gamma(link = "log"))

## predict data for offspring
calories_offspring_var_data <- 
  data.frame(cv = seq(min(calories$cv)-0.1, max(calories$cv)+0.1, length.out = 100)) %>% 
  mutate(pred = as.data.frame(predict(calories_offspring_var, newdata = ., type = "link", se = TRUE)),
         fit = exp(pred$fit),
         lowerci = exp(pred$fit - pred$se.fit * 1.96),
         upperci = exp(pred$fit + pred$se.fit * 1.96)) %>% 
  select(!pred)

#plotting ---- 

## plot lv ~ variance
p1 <-
  ggplot() +
  ggtitle("A") +
  geom_ribbon(data = calories_lv_data, 
              aes(ymin = lowerci, ymax = upperci, x = cov), 
              alpha = 0.3, fill = "#CC9BA5") +
  geom_line(data = calories_lv_data, aes(x = cov, y = fit), col = "#401D1F") +
  geom_point(data = calories, aes(x = cov, y = mean_lv, col = label)) +
  labs(x = "Caloric CoV", y = expression(bold(l[v] (m)))) +
  scale_color_manual(values = c("label" = "#0062b8", "other" = "grey20")) +
  scale_x_continuous(expand = c(0,0), limits = c(min(calories_lv_data$cov), 
                                                 max(calories_lv_data$cov))) +
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


p1 <-
  ggplot() +
  ggtitle("A") +
  #original data
  geom_point(data = calories, aes(x = cv, y = mean_lv, col = "Original Sims"), alpha = 0.5) +
  #repeated simulations
  geom_point(data = dups, aes(x = cv, y = mean_lv, col = "Repeated Sims"), alpha = 0.5) +
  scale_color_manual(name = "Data",
                     values = c("Original Sims" = "black",
                                "Repeated Sims" = "forestgreen")) +
  new_scale_color() + #new color scale for fit lines
  #original data
  geom_ribbon(data = calories_lv_data, 
              aes(ymin = lowerci, ymax = upperci, x = cv, fill = "Originals"), 
              alpha = 0.3) +
  geom_line(data = calories_lv_data, aes(x = cv, y = fit, col = "Originals")) +
  #duplicate runs
  geom_ribbon(data = dups_lv_data, 
              aes(ymin = lowerci, ymax = upperci, x = cv, fill = "Duplicates"), 
              alpha = 0.3) +
  geom_line(data = dups_lv_data, aes(x = cv, y = fit, col = "Duplicates")) +
  #combined data
  geom_ribbon(data = comb_lv_data, 
              aes(ymin = lowerci, ymax = upperci, x = cv, fill = "Combined"), 
              alpha = 0.3) +
  geom_line(data = comb_lv_data, aes(x = cv, y = fit, col = "Combined")) +
  scale_color_manual(name = "Fit",
                     values = c("Originals" = "#401D1F",
                                "Duplicates" = "darkolivegreen",
                                "Combined" = "steelblue4")) +
  scale_fill_manual(name = "Fit",
                    values = c("Originals" = "grey50",
                               "Duplicates" = "darkgreen",
                               "Combined" = "steelblue")) +
  labs(x = "cv", y = expression(bold(l[v] (m)))) +
  scale_x_continuous(expand = c(0,0), limits = c(min(calories_lv_data$cv), 
                                                 max(calories_lv_data$cv))) +
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
        legend.position = "right")

## save independently
# ggsave(p1, file = "figures/maintext/patches-vs-lv.png", 
#        width = 6, height = 3, units = "in", bg = "white", dpi = 600)

## plot speed ~ variance
p2 <-
  ggplot() +
  ggtitle("B") +
  geom_ribbon(data = calories_speed_data, 
              aes(ymin = lowerci, ymax = upperci, x = cv), 
              alpha = 0.3, fill = "#CC9BA5") +
  geom_line(data = calories_speed_data, aes(x = cv, y = fit), col = "#401D1F") +
  geom_point(data = calories, aes(x = cv, y = mean_speed, col = label)) +
  scale_color_manual(values = c("label" = "#0062b8", "other" = "grey20")) +
  labs(x = "Caloric CoV", y = "Speed (m/s)") +
  scale_x_continuous(expand = c(0, 0), limits = c(min(calories_speed_data$cv), 
                                                  max(calories_speed_data$cv))) +
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

p2 <-
  ggplot() +
  ggtitle("B") +
  #original simulations
  geom_point(data = calories, aes(x = cv, y = mean_speed, col = "Original Sims"), alpha = 0.5) +
  #repeated simulations
  geom_point(data = dups, aes(x = cv, y = mean_speed, col = "Repeated Sims"), alpha = 0.5) +
  scale_color_manual(name = "Data",
                     values = c("Original Sims" = "black",
                                "Repeated Sims" = "forestgreen")) +
  new_scale_color() + #new color scale for fit lines
  #original data
  geom_ribbon(data = calories_speed_data, 
              aes(ymin = lowerci, ymax = upperci, x = cv, fill = "Originals"), 
              alpha = 0.3) +
  geom_line(data = calories_speed_data, aes(x = cv, y = fit, col = "Originals")) +
  #duplicates
  geom_ribbon(data = dups_speed_data, 
              aes(ymin = lowerci, ymax = upperci, x = cv, fill = "Duplicates"), 
              alpha = 0.3) +
  geom_line(data = dups_speed_data, aes(x = cv, y = fit, col = "Duplicates")) +
  #combined
  geom_ribbon(data = comb_speed_data, 
              aes(ymin = lowerci, ymax = upperci, x = cv, fill = "Combined"), 
              alpha = 0.3) +
  geom_line(data = comb_speed_data, aes(x = cv, y = fit, col = "Combined")) +
  scale_color_manual(name = "Fit",
                     values = c("Originals" = "#401D1F",
                                "Duplicates" = "darkolivegreen",
                                "Combined" = "steelblue4")) +
  scale_fill_manual(name = "Fit",
                    values = c("Originals" = "grey50",
                               "Duplicates" = "darkgreen",
                               "Combined" = "steelblue")) +
  labs(x = "Number of Patches in Landscape", y = "Speed (m/s)") +
  scale_x_continuous(expand = c(0,0), limits = c(min(calories$cv), max(calories$cv))) +
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
        legend.position = "right")

## save independently
# ggsave(p2, file = "figures/maintext/patches-vs-speed.png", 
#        width = 6, height = 3, units = "in", bg = "white", dpi = 600)

## combine p1 + p2 
final <- grid.arrange(p1, p2, ncol = 2)

## save the combined figures
ggsave(final, file = "project_updates/26-09sept/caloricvariance.png", 
       width = 9, height = 4, units = "in", bg = "white", dpi = 600)

## plot offspring ~ variance
p3 <-
  ggplot() +
  ggtitle("C") +
  geom_ribbon(data = calories_offspring_data, 
              aes(ymin = lowerci, ymax = upperci, x = cv), 
              alpha = 0.3, fill = "#CC9BA5") +
  geom_line(data = calories_offspring_data, aes(x = cv, y = fit), col = "#401D1F") +
  geom_point(data = calories, aes(x = cv, y = mean_offspring, col = label)) +
  scale_color_manual(values = c("label" = "#0062b8", "other" = "grey20")) +
  labs(x = "Caloric CoV", y = "Mean Number of Offspring") +
  scale_x_continuous(expand = c(0, 0), limits = c(min(calories_offspring_data$cv), 
                                                  max(calories_offspring_data$cv))) +
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

## plot offspring variance ~ caloric variance
p4 <-
  ggplot() +
  ggtitle("D") +
  geom_ribbon(data = calories_offspring_var_data, 
              aes(ymin = lowerci, ymax = upperci, x = cv), 
              alpha = 0.3, fill = "#CC9BA5") +
  geom_line(data = calories_offspring_var_data, aes(x = cv, y = fit), col = "#401D1F") +
  geom_point(data = calories, aes(x = cv, y = offspring_var, col = label)) +
  scale_color_manual(values = c("label" = "#0062b8", "other" = "grey20")) +
  labs(x = "Caloric CoV", y = "Variance in Mean Number of Offspring") +
  scale_x_continuous(expand = c(0, 0), limits = c(min(calories_offspring_var_data$cv), 
                                                  max(calories_offspring_var_data$cv))) +
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

# habitat figures ----

## habitat calories for poster
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
                               n_points = 500,
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
  facet_wrap(~cal_var, ncol = 5, labeller = as_labeller(function(x) paste("Calorie CoV", x))) +
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


ggplot() +
  geom_point(data = calories, aes(x = cv, y = patches)) +
  labs(x = "Caloric CoV", y = "Number of Patches Consumed") +
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

caloric_var <- list.files(path = "simulations/sensitivity/variance", 
                          pattern = "prey_details\\.Rds$", 
                          full.names = TRUE) %>% 
  map(~ {
    data_list <- readRDS(.x) %>% 
      bind_rows() %>% 
      # na.omit() %>% 
      group_by(cal_var) %>% 
      filter(generation >= max(generation) - 10) %>% 
      select(cal_var,
             patches,
             offspring,
             lv,
             speed)
    return(data_list)
  }) %>%
  list_rbind()

cal_VAR <- caloric_var %>% 
  group_by(cal_var) %>% 
  filter(!cal_var %in% c(2.55, 2.65, 2.75, 2.85, 2.95, 3.05, 3.15, 3.25,  3.35, 3.45)) %>% 
  summarise(patches_cv = var(patches),
            offspring_cv = var(offspring),
            lv_cv = var(lv),
            speed_cv = var(speed))

ggplot() +
  geom_point(data = cal_VAR, aes(x = cal_var, y = offspring_cv)) + 
  theme_bw()


df_long <- cal_VAR %>%
  pivot_longer(
    cols = c(patches_cv, offspring_cv, lv_cv, speed_cv),
    names_to = "variable",
    values_to = "value"
  )

ggplot(df_long, aes(x = cal_var, y = value)) +
  geom_point() +
  labs(y = "Variance", y = "Caloric Variance") +
  facet_wrap(~ variable, scales = "free_y", labeller = as_labeller(c(
    patches_cv = "Number of Consumed Patches",
    offspring_cv = "Number of Offspring",
    lv_cv = "Ballistic Length Scale",
    speed_cv = "Speed"))) +
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