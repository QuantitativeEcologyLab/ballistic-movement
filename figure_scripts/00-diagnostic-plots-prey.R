# preamble

# Set the working directory
setwd("~/hdrive/GitHub/ballistic-movement")

# load packages
library(tidyverse)
library(data.table)
library(mgcv)

# Source the functions 
source("figure_scripts/01-custom-ggplot-theme.R")
source("figure_scripts/00-diagnostic-functions.R")

#..............................................................................
# quick checks ----
#..............................................................................

#load in your data
prey_details <- readRDS('simulations/prey_results/num_patches/30800p_prey_details.Rds')

# make data sets compatible with ggplot2
prey_details_df <- bind_rows(prey_details)

#FIG <- 
  simple.diag(prey_details_df)

# save figures 
ggsave(FIG, file = "project_updates/26-03mar/40400_patches.png", width = 20, height = 10, units = "in", dpi = 600)


#plot food landscape
FOOD <- readRDS('simulations/prey_results/num_patches/habitats/50000_patches.Rds')

food <- data.frame(x = FOOD$x,
                   y = FOOD$y,
                   cal = FOOD$marks)

p1 <- 
  ggplot() +
  geom_point(data = food, aes(x = x, y = y), size = 0.5) +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  theme.qel()

ggsave(p1, file = "project_updates/26-02feb/50000patches_habitat.png", width = 5, height = 5, units = "in", dpi = 600, bg = "white")

#..............................................................................
# patches vs bls ----
#..............................................................................

combined_df <- list.files(path = "simulations/prey_results/num_patches", 
                          pattern = "prey_details\\.Rds$", 
                          full.names = TRUE) %>% 
  map(~ {
    id_val <- as.numeric(str_extract(basename(.x), "\\d+"))
    
    data_list <- readRDS(.x)
    
    map_dfr(data_list, ~ mutate(.x, file_id = id_val))
  }) %>%
  list_rbind()

#p1 <- 
  combined_df %>% 
  group_by(file_id) %>% 
  filter(file_id <= 50000) %>%
  filter(!file_id %in% c(40400,1000)) %>%
  filter(generation >= max(generation - 500)) %>% 
  summarise(mean_lv = mean(lv),
            npatches = mean(file_id)) %>% 
  ggplot(aes(x = file_id, y = mean_lv)) +
  geom_point() + 
  geom_smooth() +
  labs(x = "Number of Patches in Landscape", y = "Mean ballistic length scale") +
  theme_bw()

  df1 <- bind_rows(readRDS("simulations/prey_results/num_patches/40400p_prey_details.Rds")) %>% 
    filter(generation <= 4000)
  
  df2 <- bind_rows(readRDS("simulations/prey_results/num_patches/40400p_prey_details2.Rds")) 
  
  ggplot() +
    stat_summary(data = df1, aes(x = generation, y = lv), fun = "mean", geom = "point", col = "blue", alpha = 0.3, size = 1) +
    stat_summary(data = df2, aes(x = generation, y = lv), fun = "mean", geom = "point", col = "red", alpha = 0.3, size = 1)

#ggsave(p1, file = "project_updates/26-02feb/calsvslv.png", width = 6, height = 4, units = "in", bg = "white")

total <- combined_df %>% 
  filter(file_id <= 50000) %>% 
  mutate(
    patchcal = case_when(
      file_id == "50000" ~ 45, 
      file_id == "45200" ~ 47, 
      file_id == "35600" ~ 61,
      file_id == "30800" ~ 67,
      file_id == "21200" ~ 98,
      file_id == "16400" ~ 122,
      file_id == "11600" ~ 196,
      file_id == "6800" ~ 309,
      file_id == "4400" ~ 476,
      file_id == "2000" ~ 1013),
    totcal = patchcal*file_id) 

#p1 <- 
  total %>% 
    group_by(totcal) %>% 
    filter(generation >= max(generation - 500)) %>% 
    select(lv, totcal) %>% 
    mutate(mean_lv = mean(lv)) %>% 
    ggplot(aes(x = totcal, y = mean_lv)) +
    geom_point() + 
    labs(x = "Total Calories in Landscape", y = "Mean ballistic length scale") +
    theme_bw() 

#..............................................................................
# modelling cal/patch versus number of patches ----
#..............................................................................

#creating this model will help predict the number of calories needed for simulations 
#with differing number of patches 
# ONLY FOR MASS = 105500

cal_list <- c(45,47,61,67,78,98,121.5,196,309,476,1013)

patches <- c(50000,45200,35600,30800,26000,21200,16400,11600,6800,4400,2000)

cal_df <- bind_cols(cal = cal_list, npatches = patches)

FIT <- gam(log10(cal) ~ log10(npatches), family = Gamma(link = "log"), data = cal_df)

gam.check(FIT)
summary(FIT)
plot(fitted.values(FIT), residuals(FIT))

ggplot(cal_df, aes(x = npatches, y = cal)) +
  geom_point() +
  geom_smooth(method = "lm") + 
  scale_y_log10() +
  scale_x_log10() +
  theme_classic()

newdata <- data.frame(npatches = seq(1000, 50000, length.out = 12))

newdata <- newdata %>% 
  filter(!(npatches %in% patches))

preds <- predict(FIT, newdata = newdata, type = "response") 

plot(newdata$npatches, 10^preds)

newdata %>% 
  mutate(pred_cal = 10^preds)

