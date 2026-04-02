# preamble

# Set the working directory
setwd("~/hdrive/GitHub/ballistic-movement")

# load packages
library(tidyverse)
library(data.table)
library(mgcv)

# Source the functions 
source("figure_scripts/00-diagnostic-functions.R")

#..............................................................................
# quick checks ----
#..............................................................................

#load in your data
details_df <- bind_rows(readRDS("simulations/prey_results/clustering/3mu_prey_details.Rds"))

#p1 <-
  explore.gen(details_df)

s# save figures 
ggsave(p1, file = "project_updates/26-03mar/550p.png", width = 20, height = 10, units = "in", dpi = 600)

p2 <-
  explore.var(details_df) 

ggsave(p2, file = "project_updates/26-03mar/250p.png", width = 20, height = 5, units = "in", dpi = 600)


# compare different simulations
df1 <- bind_rows(readRDS("simulations/prey_results/num_patches/500p_prey_details.Rds")) 

df2 <- bind_rows(readRDS("simulations/prey_results/num_patches/475p_prey_details.Rds")) 

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
