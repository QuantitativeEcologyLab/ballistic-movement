# script for making plot of ballistic length scale vs degree of patch clustering in the habitat

#..............................................................................
# load libraries ----

library(tidyverse)
library(patchwork)
library(viridis)
library(gridExtra)

source("simulation_scripts/01-prey-functions.R")

#..............................................................................

#load data from folder
clust <- list.files(path = "simulations/prey_results/clustering", 
                    pattern = "prey_details\\.Rds$", 
                    full.names = TRUE) %>% 
  map(~ {
    #add id for clustering level from file name
    id_val <- as.numeric(str_extract(basename(.x), "\\d+"))
    
    data_list <- readRDS(.x) %>% map_dfr(., ~ mutate(.x, mu = id_val))
    
    #summarise all the data 
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

# add column for labeling
clust <- clust %>% 
  mutate(label = case_when(
    mu == 1 ~ "label",
    TRUE ~ "other"
  ))

#fit model to data
clust_lv <- glm(mean_lv ~ mu,
              data = clust, 
              family = Gamma(link = "log"))

#predict data from model
clust_lv_data <- 
  data.frame(mu = seq(min(clust$mu) - 2, max(clust$mu) + 2, length.out = 100)) %>% 
  mutate(pred = as.data.frame(predict(clust_lv, newdata = ., type = "link", se = TRUE)),
         fit = exp(pred$fit),
         lowerci = exp(pred$fit - pred$se.fit * 1.96),
         upperci = exp(pred$fit + pred$se.fit * 1.96)) %>% 
  select(!pred)

p1 <-
  ggplot() +
  ggtitle("A") +
  geom_ribbon(data = clust_lv_data, 
              aes(ymin = lowerci, ymax = upperci, x = mu), 
              alpha = 0.3, fill = "#CC9BA5") +
  geom_line(data = clust_lv_data, aes(x = mu, y = fit), col = "#401D1F") +
  geom_point(data = clust, aes(x = mu, y = mean_lv, col = label)) +
  labs(x = "Degree of Clustering", y = expression(bold(l[v] (m)))) +
  scale_color_manual(values = c("label" = "#0062b8", "other" = "grey20")) +
  scale_x_continuous(expand = c(0,0), limits = c(min(clust_lv_data$mu), 
                                                 max(clust_lv_data$mu))) +
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

# ggsave(p1, file = "figures/maintext/patches-vs-lv.png", 
#        width = 6, height = 3, units = "in", bg = "white", dpi = 600)

#model speed
clust_speed <- glm(mean_speed ~ mu,
                 data = clust, 
                 family = Gamma(link = "log"))

clust_speed_data <- 
  data.frame(mu = seq(min(clust$mu) - 2, max(clust$mu) + 2, length.out = 100)) %>% 
  mutate(pred = as.data.frame(predict(clust_speed, newdata = ., type = "link", se = TRUE)),
         fit = exp(pred$fit),
         lowerci = exp(pred$fit - pred$se.fit * 1.96),
         upperci = exp(pred$fit + pred$se.fit * 1.96)) %>% 
  select(!pred)

p2 <-
  ggplot() +
  ggtitle("B") +
  geom_ribbon(data = clust_speed_data, 
              aes(ymin = lowerci, ymax = upperci, x = mu), 
              alpha = 0.3, fill = "#CC9BA5") +
  geom_line(data = clust_speed_data, aes(x = mu, y = fit), col = "#401D1F") +
  geom_point(data = clust, aes(x = mu, y = mean_speed, col = label)) +
  scale_color_manual(values = c("label" = "#0062b8", "other" = "grey20")) +
  labs(x = "Degree of Clustering", y = "Speed (m/s)") +
  scale_x_continuous(expand = c(0,0), limits = c(min(clust_speed_data$mu), max(clust_speed_data$mu))) +
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

#create combined plot
final <- grid.arrange(p1, p2, ncol = 2)

ggsave(final, file = "figures/03-clustering.png", 
       width = 9, height = 4, units = "in", bg = "white", dpi = 600)

#habitat clustering for poster
base_food <- readRDS("simulations/prey_results/patches/habitats/500_patches.Rds")

mass_prey <- 105500

clustvals <- c(1, 5, 13, 21, 30)

clust_res <- vector("list", length(clustvals))
for(idx in seq_along(clustvals)) {
  i <- clustvals[idx]
  
  if(i == 1) {
    food <- base_food   # <-- force identical landscape
  } else {
    success <- FALSE
    
    set.seed(123 + idx)
    
    while(!success){
      FOOD <- try({makeHabitat(mass_prey,
                               r = 1,
                               mu = i,
                               target_n = 500,
                               cal = 4000)},
                  silent = TRUE)
      
      success <- !inherits(FOOD, "try-error")
    }
    
    food <- as.data.frame(FOOD)
  }
  
  clust_res[[idx]] <- data.frame(x = food$x,
                                    y = food$y,
                                    mu = i)
}

clust_res <- bind_rows(clust_res)

clust_habitats <-
  clust_res %>% 
  ggplot(aes(x = x, y = y)) +
  geom_rect(data = filter(clust_res, mu == 1),
            xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf,
            fill = "#ebf6ff",
            inherit.aes = FALSE) +
  geom_point(col = "#461300",size = 0.5) +
  facet_wrap(~mu, ncol = 8, labeller = as_labeller(function(x) paste("Clustering Coefficient", x))) +
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

# ggsave(habitats, file = "presentations/poster-components/figures/habitat-clustering.png", 
#        width = 10, height = 1.9, units = "in", bg = "white", dpi = 600)

# FIG <- grid.arrange(clust_habitats, final, nrow = 2, heights = c(1,2))
# 
# ggsave(FIG, file = "presentations/poster-components/figures/landscape-clustering.png", 
#        width = 9, height = 5, units = "in", dpi = 600, bg = "white")
