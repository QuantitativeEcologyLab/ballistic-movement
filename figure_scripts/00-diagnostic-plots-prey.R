# preamble

# Set the working directory
setwd("~/hdrive/GitHub/ballistic-movement")

# load packages
library(tidyverse)
library(data.table)
library(mgcv)
#library(khroma)

# Source the functions 
source("figure_scripts/00-diagnostic-functions.R")

#..............................................................................
# quick checks ----
#..............................................................................

#load in your data
details_df <- bind_rows(readRDS("simulations/sensitivity/variance/0.5cv_prey_details.Rds"))

p1 <-
  explore.gen(details_df)

# save figures 
ggsave(p1, file = "project_updates/26-08aug/0.5cv_gen2.png", width = 14.8, height = 6.74, units = "in", dpi = 600)

p2 <-
  explore.var(details_df) 

ggsave(p2, file = "project_updates/26-08aug/0.5cv_var2.png", width = 20, height = 5, units = "in", dpi = 600)


# compare different simulations
df1 <- bind_rows(readRDS("simulations/sensitivity/variance/0.1cv_prey_details.Rds")) 

df2 <- bind_rows(readRDS("simulations/sensitivity/variance/0.5cv_prey_details.Rds")) 

ggplot() +
  stat_summary(data = df1, aes(x = generation, y = lv), fun = "mean", geom = "point", col = "blue", alpha = 0.3, size = 1) +
  stat_summary(data = df2, aes(x = generation, y = lv), fun = "mean", geom = "point", col = "red", alpha = 0.3, size = 1)


ggplot() + 
    geom_point(data = df1, aes(x = tau_v, y = tau_p), col = "blue") + 
    geom_point(data = df2, aes(x = tau_v, y = tau_p), col = "red")


df1 <- bind_rows(readRDS("simulations/prey_results/patches/250p_prey_details.Rds")) 

df2 <- bind_rows(readRDS("simulations/prey_results/patches/800p_prey_details.Rds"))

ggplot() +
  geom_point(data = df1, aes(x = lv, y = speed, col = generation), size = 0.5, col = "steelblue") +
  geom_point(data = df2, aes(x = lv, y = speed, col = generation), size = 0.5, col = "palevioletred")


ggplot(details_df) + 
  geom_point(aes(x = speed, y = scale(patches)), col = "palevioletred") + 
  geom_point(aes(x = speed, y = scale(cal_net)), col = "steelblue")

ggplot(details_df) + 
  geom_point(aes(x = lv, y = scale(patches)), col = "palevioletred") + 
  geom_point(aes(x = lv, y = scale(cal_net)), col = "steelblue")

#save calorie variance plots
calories <- list.files(path = "simulations/prey_results/calorie-variance", 
                       pattern = "prey_details\\.Rds$", 
                       full.names = TRUE) %>% 
  map(~ {
    id_val <- as.numeric(str_extract(basename(.x), "\\d+"))
    
    fig <- readRDS(.x) %>% 
      bind_rows() %>% 
      explore.gen()
    
    ggsave(fig, file = paste0("project_updates/26-04apr/calorie-variance/", paste(id_val, sep = "_"), "var_genplot.png"),
           width = 14.8, height = 6.74, units = "in")
    
    fig2 <- readRDS(.x) %>% 
      bind_rows() %>% 
      explore.var()
    
    ggsave(fig2, file = paste0("project_updates/26-04apr/calorie-variance/", paste(id_val, sep = "_"), "var_varplot.png"),
           width = 14.8, height = 3.4, units = "in")
  })

food1 <- as.data.frame(readRDS("simulations/sensitivity/variance/habitats/0.1cv.Rds"))

range(food1$marks)

ggplot(food1, aes(x = x, y = y, color = marks)) +
  geom_point(size = 3) +
  scale_color_viridis_c(option = "magma") +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  labs(color = "Calories") +
  theme_void() +
  theme() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        plot.margin = margin(5,5,5,5))

food2 <- as.data.frame(readRDS("simulations/sensitivity/variance/habitats/0.5cv.Rds"))

range(food2$marks)

ggplot(food1, aes(x = x, y = y, color = marks)) +
  geom_point(size = 3) +
  scale_color_viridis_c(option = "magma") +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  labs(color = "Calories") +
  theme_void() +
  theme() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        plot.margin = margin(5,5,5,5))

food3 <- as.data.frame(makeHabitat(mass = 105500,
                                   mu = 1,
                                   cv = 0.15,
                                   n_points = 500,
                                   cal = 4000,
                                   seed = 123))

range(food3$marks)
