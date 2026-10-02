# script for making plots for simulations testing the number of patches in landscape

#..............................................................................
# load libraries ----

library(tidyverse)
library(patchwork)
library(viridis)
library(gridExtra)
library(propagate)
library(investr)
library(ggnewscale)

source("simulation_scripts/01-prey-functions.R")

#..............................................................................

patches <- list.files(path = "simulations/prey_results/patches", 
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
                var_lv = var(lv),
                mean_speed = mean(speed),
                var_speed = var(speed),
                babies = mean(offspring),
                offspring_var = var(offspring))
    return(data_list)
  }) %>% 
  list_rbind()

# repeated simulations
dups <- list.files(path = "simulations/prey_results/patches/duplicates", 
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
                var_lv = var(lv),
                mean_speed = mean(speed),
                var_speed = var(speed),
                babies = mean(offspring),
                offspring_var = var(offspring))
    return(data_list)
  }) %>% 
  list_rbind()

# combine originals and repeats
comb <- bind_rows(patches, dups)

# add column for labeling
# patches <- patches %>% 
#   mutate(label = case_when(
#     tot_patch == 500 ~ "label",
#     TRUE ~ "other"
#   ))

#model lv from original simulations
patches_lv <- glm(mean_lv ~ tot_patch,
              data = patches, 
              family = Gamma(link = "log"))

patches_lv_data <- 
  data.frame(tot_patch = seq(min(patches$tot_patch)-100, 
                             max(patches$tot_patch)+100, 
                             length.out = 100)) %>% 
  mutate(pred = as.data.frame(predict(patches_lv, newdata = ., type = "link", se = TRUE)),
         fit = exp(pred$fit),
         lowerci = exp(pred$fit - pred$se.fit * 1.96),
         upperci = exp(pred$fit + pred$se.fit * 1.96)) %>% 
  select(!pred)

#model lv from duplicated simulations
dups_lv <- glm(mean_lv ~ tot_patch,
                  data = dups, 
                  family = Gamma(link = "log"))

dups_lv_data <- 
  data.frame(tot_patch = seq(min(patches$tot_patch)-100, 
                             max(patches$tot_patch)+100, 
                             length.out = 100)) %>% 
  mutate(pred = as.data.frame(predict(dups_lv, newdata = ., type = "link", se = TRUE)),
         fit = exp(pred$fit),
         lowerci = exp(pred$fit - pred$se.fit * 1.96),
         upperci = exp(pred$fit + pred$se.fit * 1.96)) %>% 
  select(!pred)

#model lv for all simulations
comb_lv <- glm(mean_lv ~ tot_patch,
                  data = comb, 
                  family = Gamma(link = "log"))

comb_lv_data <- 
  data.frame(tot_patch = seq(min(comb$tot_patch)-100, 
                             max(comb$tot_patch)+100, 
                             length.out = 100)) %>% 
  mutate(pred = as.data.frame(predict(comb_lv, newdata = ., type = "link", se = TRUE)),
         fit = exp(pred$fit),
         lowerci = exp(pred$fit - pred$se.fit * 1.96),
         upperci = exp(pred$fit + pred$se.fit * 1.96)) %>% 
  select(!pred)

# patches versus lv
p1 <-
  ggplot() +
    ggtitle("A") +
    #original data
    geom_point(data = patches, aes(x = tot_patch, y = mean_lv, col = "Original Sims"), alpha = 0.5) +
    #repeated simulations
    geom_point(data = dups, aes(x = tot_patch, y = mean_lv, col = "Repeated Sims"), alpha = 0.5) +
    scale_color_manual(name = "Data",
                       values = c("Original Sims" = "black",
                                  "Repeated Sims" = "forestgreen")) +
    new_scale_color() + #new color scale for fit lines
    #original data
    geom_ribbon(data = patches_lv_data, 
                aes(ymin = lowerci, ymax = upperci, x = tot_patch, fill = "Originals"), 
                alpha = 0.3) +
    geom_line(data = patches_lv_data, aes(x = tot_patch, y = fit, col = "Originals")) +
    #duplicate runs
    geom_ribbon(data = dups_lv_data, 
                  aes(ymin = lowerci, ymax = upperci, x = tot_patch, fill = "Duplicates"), 
                  alpha = 0.3) +
    geom_line(data = dups_lv_data, aes(x = tot_patch, y = fit, col = "Duplicates")) +
    #combined data
    geom_ribbon(data = comb_lv_data, 
                  aes(ymin = lowerci, ymax = upperci, x = tot_patch, fill = "Combined"), 
                  alpha = 0.3) +
    geom_line(data = comb_lv_data, aes(x = tot_patch, y = fit, col = "Combined")) +
    scale_color_manual(name = "Fit",
                       values = c("Originals" = "#401D1F",
                                  "Duplicates" = "darkolivegreen",
                                  "Combined" = "steelblue4")) +
    scale_fill_manual(name = "Fit",
                       values = c("Originals" = "grey50",
                                  "Duplicates" = "darkgreen",
                                  "Combined" = "steelblue")) +
    labs(x = "Number of Patches in Landscape", y = expression(bold(l[v] (m)))) +
    scale_x_continuous(expand = c(0,0), limits = c(min(patches_lv_data$tot_patch), 
                                                   max(patches_lv_data$tot_patch))) +
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

# ggsave(p, file = "project_updates/26-04apr/cal-lvspeed.png", 
#        width = 6, height = 3, units = "in", bg = "white", dpi = 600)

# model speed

#original data
patches_nls <- nls(mean_speed ~ (a * tot_patch) / (b + tot_patch),
               start = list(a = 1, b = 1),
               data = patches)

#investr method (less computation power)
preds <- predFit(patches_nls, 
                 newdata = data.frame(
                            tot_patch = 
                              seq(min(patches$tot_patch)-100, 
                                  max(patches$tot_patch)+100, 
                                  length.out = 100)),
                 interval = "confidence", level = 0.95)

cis <- data.frame(tot_patch = 
                  seq(min(patches$tot_patch)-100, 
                      max(patches$tot_patch)+100, 
                      length.out = 100),
                  fit = preds[, "fit"],
                  lowerci = preds[, "lwr"],
                  upperci = preds[, "upr"])

#duplicates
dups_nls <- nls(mean_speed ~ (a * tot_patch) / (b + tot_patch),
                   start = list(a = 1, b = 1),
                   data = dups)

dups_preds <- predFit(dups_nls, 
                 newdata = data.frame(
                   tot_patch = 
                     seq(min(patches$tot_patch)-100, 
                         max(patches$tot_patch)+100, 
                         length.out = 100)),
                 interval = "confidence", level = 0.95)

dups_cis <- data.frame(tot_patch = 
                    seq(min(patches$tot_patch)-100, 
                        max(patches$tot_patch)+100, 
                        length.out = 100),
                  fit = dups_preds[, "fit"],
                  lowerci = dups_preds[, "lwr"],
                  upperci = dups_preds[, "upr"])

#combined
comb_nls <- nls(mean_speed ~ (a * tot_patch) / (b + tot_patch),
                start = list(a = 1, b = 1),
                data = comb)

comb_preds <- predFit(comb_nls, 
                      newdata = data.frame(
                        tot_patch = 
                          seq(min(patches$tot_patch)-100, 
                              max(patches$tot_patch)+100, 
                              length.out = 100)),
                      interval = "confidence", level = 0.95)

comb_cis <- data.frame(tot_patch = 
                         seq(min(patches$tot_patch)-100, 
                             max(patches$tot_patch)+100, 
                             length.out = 100),
                       fit = comb_preds[, "fit"],
                       lowerci = comb_preds[, "lwr"],
                       upperci = comb_preds[, "upr"])

# patches vs speed
p2 <-
  ggplot() +
    ggtitle("B") +
    #original simulations
    geom_point(data = patches, aes(x = tot_patch, y = mean_speed, col = "Original Sims"), alpha = 0.5) +
    #repeated simulations
    geom_point(data = dups, aes(x = tot_patch, y = mean_speed, col = "Repeated Sims"), alpha = 0.5) +
    scale_color_manual(name = "Data",
                       values = c("Original Sims" = "black",
                                  "Repeated Sims" = "forestgreen")) +
    new_scale_color() + #new color scale for fit lines
    #original data
    geom_ribbon(data = cis, 
               aes(ymin = lowerci, ymax = upperci, x = tot_patch, fill = "Originals"), 
               alpha = 0.3) +
    geom_line(data = cis, aes(x = tot_patch, y = fit, col = "Originals")) +
    #duplicates
    geom_ribbon(data = dups_cis, 
                aes(ymin = lowerci, ymax = upperci, x = tot_patch, fill = "Duplicates"), 
                alpha = 0.3) +
    geom_line(data = dups_cis, aes(x = tot_patch, y = fit, col = "Duplicates")) +
    #combined
    geom_ribbon(data = comb_cis, 
                aes(ymin = lowerci, ymax = upperci, x = tot_patch, fill = "Combined"), 
                alpha = 0.3) +
    geom_line(data = comb_cis, aes(x = tot_patch, y = fit, col = "Combined")) +
    scale_color_manual(name = "Fit",
                       values = c("Originals" = "#401D1F",
                                  "Duplicates" = "darkolivegreen",
                                  "Combined" = "steelblue4")) +
    scale_fill_manual(name = "Fit",
                      values = c("Originals" = "grey50",
                                 "Duplicates" = "darkgreen",
                                 "Combined" = "steelblue")) +
    labs(x = "Number of Patches in Landscape", y = "Speed (m/s)") +
    scale_x_continuous(expand = c(0,0), limits = c(min(cis$tot_patch), max(cis$tot_patch))) +
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

# ggsave(p2, file = "figures/maintext/patches-vs-speed.png", 
#        width = 6, height = 3, units = "in", bg = "white", dpi = 600)

#combine plots
#final <- grid.arrange(p1, p2, ncol = 2)
#final <- (p1 + p2) + plot_layout(guides = "collect")

# ggsave(final, file = "figures/02-npatches.png", 
#        width = 9, height = 4, units = "in", bg = "white", dpi = 600)

# patches vs lv variance
p3 <- 
  ggplot() +
  geom_point(data = comb, aes(x = tot_patch, y = var_lv)) + 
  labs(x = "Number of Patches in Landscape", y = expression(bold("Variance in" ~ l[v] ~ "(m)"))) +
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

# patches vs offspring variance
p4 <- 
  ggplot() +
  geom_point(data = comb, aes(x = tot_patch, y = var_speed)) + 
  labs(x = "Number of Patches in Landscape", y = "Variance in Speed") +
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

# patches vs mean offspring
p5 <- 
  ggplot() +
  geom_point(data = comb, aes(x = tot_patch, y = babies)) + 
  labs(x = "Number of Patches in Landscape", y = "Mean Offspring") +
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

# patches vs offspring variance
p6 <- 
  ggplot() +
  geom_point(data = comb, aes(x = tot_patch, y = offspring_var)) + 
  labs(x = "Number of Patches in Landscape", y = "Mean Offspring") +
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

#representative habitats for poster
base_food <- readRDS("simulations/prey_results/patches/habitats/500_patches.Rds")

mass_prey <- 105500

den <- c(250, 500, 750, 900, 1050)

res <- vector("list", length(den))
for(idx in seq_along(den)) {
  i <- den[idx]
  
  if(i == 500) {
    food <- base_food   
  } else {
    success <- FALSE
    
    set.seed(123 + idx)
    
    while(!success){
      FOOD <- try({makeHabitat(mass_prey,
                               r = 1,
                               mu = 1,
                               n_points = i,
                               cal = 4000)},
                  silent = TRUE)
      
      success <- !inherits(FOOD, "try-error")
    }
    
    food <- as.data.frame(FOOD)
  }
  
  res[[idx]] <- data.frame(x = food$x,
                           y = food$y,
                           target_n = i)
}
res <- bind_rows(res)

habitats <- res %>% 
  ggplot(aes(x = x, y = y)) +
  geom_rect(data = filter(res, target_n == 500),
            xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf,
            fill = "#ebf6ff",
            inherit.aes = FALSE) +
  geom_point(col = "#461300", size = 0.5) +
  facet_wrap(~target_n, ncol = 8, labeller = as_labeller(function(x) paste(x, "patches"))) +
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
        panel.background = element_rect(fill = "transparent"))

# ggsave(habitats, file = "presentations/poster-components/figures/habitat-density.png", 
#        width = 10, height = 1.9, units = "in", bg = "white", dpi = 600)

#combine all
# FIG <- grid.arrange(habitats, final, nrow = 2, heights = c(1,2))
# 
# ggsave(FIG, file = "presentations/poster-components/figures/landscape-density.png", 
#        width = 9, height = 5, units = "in", dpi = 600, bg = "white")
