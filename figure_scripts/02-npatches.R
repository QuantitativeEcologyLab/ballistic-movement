# script for making plots for simulations testing the number of patches in landscape

#..............................................................................
# load libraries ----

library(tidyverse)
library(patchwork)
library(viridis)
library(gridExtra)
library(propagate)

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
                mean_speed = mean(speed))
    return(data_list)
  }) %>% 
  list_rbind()

# add column for labeling
patches <- patches %>% 
  mutate(label = case_when(
    tot_patch == 500 ~ "label",
    TRUE ~ "other"
  ))

#model lv
patches_lv <- glm(mean_lv ~ tot_patch,
              data = patches, 
              family = Gamma(link = "log"))

patches_lv_data <- 
  data.frame(tot_patch = seq(min(patches$tot_patch)+100, 
                             max(patches$tot_patch)-100, 
                             length.out = 100)) %>% 
  mutate(pred = as.data.frame(predict(patches_lv, newdata = ., type = "link", se = TRUE)),
         fit = exp(pred$fit),
         lowerci = exp(pred$fit - pred$se.fit * 1.96),
         upperci = exp(pred$fit + pred$se.fit * 1.96)) %>% 
  select(!pred)

p1 <-
  ggplot() +
  ggtitle("A") +
  geom_ribbon(data = patches_lv_data, 
              aes(ymin = lowerci, ymax = upperci, x = tot_patch), 
              alpha = 0.3, fill = "#CC9BA5") +
  geom_line(data = patches_lv_data, aes(x = tot_patch, y = fit), col = "#401D1F") +
  geom_point(data = patches, aes(x = tot_patch, y = mean_lv, col = label)) +
  scale_color_manual(values = c("label" = "#0062b8", "other" = "grey20")) +
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
        legend.position = "none")

# ggsave(p, file = "project_updates/26-04apr/cal-lvspeed.png", 
#        width = 6, height = 3, units = "in", bg = "white", dpi = 600)

# model speed

patches_nls <- nls(mean_speed ~ (a * tot_patch) / (b + tot_patch),
               start = list(a = 1, b = 1),
               data = patches)

#makes output silent
suppressMessages(
  invisible(capture.output(
  #uses simulation to get confidence intervals
    preds <- predictNLS(patches_nls, 
                        newdata = data.frame(
                                  tot_patch = 
                                   seq(min(patches$tot_patch)-100, 
                                       max(patches$tot_patch)+100, 
                                       length.out = 100)))
  )))

cis <- data.frame(tot_patch = 
                    seq(min(patches$tot_patch)-100, 
                        max(patches$tot_patch)+100, 
                        length.out = 100),
                  fit = preds$summary[, "mean.1"],
                  lowerci = preds$summary[, "2.5%"],
                  upperci = preds$summary[, "97.5%"])

p2 <-
  ggplot() +
  ggtitle("B") +
  geom_ribbon(data = cis, 
              aes(ymin = lowerci, ymax = upperci, x = tot_patch), 
              alpha = 0.3, fill = "#CC9BA5") +
  geom_line(data = cis, aes(x = tot_patch, y = fit), col = "#401D1F") +
  geom_point(data = patches, aes(x = tot_patch, y = mean_speed, col = label)) +
  scale_color_manual(values = c("label" = "#0062b8", "other" = "grey20")) +
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
        legend.position = "none")

# ggsave(p2, file = "figures/maintext/patches-vs-speed.png", 
#        width = 6, height = 3, units = "in", bg = "white", dpi = 600)

#combine plots
final <- grid.arrange(p1, p2, ncol = 2)

ggsave(final, file = "figures/02-npatches.png", 
       width = 9, height = 4, units = "in", bg = "white", dpi = 600)

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
