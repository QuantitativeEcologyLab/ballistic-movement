#preamble

# Set the working directory
setwd("~/hdrive/GitHub/ballistic-movement")

# Import necessary packages
library(ggplot2)
library(dplyr)
library(gridExtra)
library(patchwork)
library(purrr)

# Source the functions (ensure 'functions.R' is available in the working directory)
source("scripts/functions.R")

#----------------------------------------------------------------------
# make figures for results section----
#----------------------------------------------------------------------

#load in your data
load('prey_results/189500g_prey_res.Rda')
load('prey_results/189500g_prey_details.Rda')

prey_res_df <- do.call(rbind, prey_res)
prey_details_df <- do.call(rbind, prey_details)

# relative change in lv ~ gen
PREY_LV <- prey_res_df$lv[1]
prey_res_df$rel.lv <- prey_res_df$lv/PREY_LV
prey_res_df$rel_var <- prey_res_df$var / (PREY_LV^2)
prey_res_df$rel_sd <- sqrt(prey_res_df$rel_var)

prey_summary <- prey_details_df %>%
  group_by(generation) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop")

prey_details_df <- prey_details_df %>% 
  group_by(generation) %>% 
  mutate(lv_mean = mean(lv))

lastgens <- prey_details_df %>% 
  filter(generation >= max(generation) - 199)

p1 <-
  ggplot(prey_res_df, aes(x = generation, y = rel.lv)) +
  geom_line(color = "#4278aaff", linewidth = 0.5) +
  geom_hline(yintercept = 1, color = 'grey30', linetype = "dashed") +
  geom_ribbon(aes(ymin = rel.lv - rel_sd,
                  ymax = rel.lv + rel_sd),
              fill = "#4278aaff", alpha = 0.3) +
  labs(y = "Relative change in ballistic lengthscale (m)",x = "Generation") +
  theme.qel()

p2 <-
  ggplot() +
  geom_point(dat = prey_details_df, aes(x = generation, y = lv), col = "grey70", alpha = 0.8, size = 0.1) +
  geom_line(dat = prey_details_df, aes(x = generation, y = lv_mean), col = "#4278aaff", linewidth = 0.5) +
  labs(y = "Ballistic lengthscale (m)", x = "Generation") +
  theme.qel()

p3 <-
  ggplot(lastgens, aes(x = lv, y = patches)) +
  geom_point(color = "#4278aaff") +
  labs(y = "Patches visited", x = "Ballistic lengthscale (m)") +
  theme.qel()
print(p3)

p4 <- prey_details_df %>%
  filter(generation >= max(generation) - 199) %>% 
  ggplot(aes(x = lv, y = speed)) +
  geom_point(color = "#4278aaff") +
  labs(y = "Speed (m/s)", x = "Ballistic lengthscale (m)") +
  theme.qel()
print(p4)

p5 <- prey_details_df %>% 
  filter(generation >= max(generation) - 20) %>% 
  ggplot(aes(x = lv, y = offspring)) +
  geom_point(color = "#4278aaff") +
  theme.qel()
print(p5)

p6 <-
  ggplot(prey_details_df, aes(x = speed, y = offspring)) +
  geom_point(color = "#4278aaff") +
  theme.qel()

FIG <- grid.arrange(p1, p2, p3, p4, p5, p6)

#ggsave(FIG, width = 6.86, height = 7.36, units = "in", dpi = 600, bg = "white", file="prey_results/figures/189500g_results_test.png")

binned_lv <- prey_details_df %>% 
  mutate(lv_bin = floor(lv/0.5)*0.5) %>% 
  group_by(lv_bin) %>% 
  mutate(mean_patch = mean(patches), .groups = "drop") %>% 
  ggplot(aes(x = lv, y = mean_patch)) +
  geom_point(color = "#4278aaff") +
  labs(y = "Offspring", x = "Ballistic lengthscale (m)") +
  theme.qel()
print(binned_lv)


## analysis across mass spectrum
load.prey.details <- function(file_path) {
  #load the .Rda file into a new environment to avoid cluttering global environment
  env <- new.env() 
  
  load(file_path, envir = env)
  
  #assuming the loaded object is a list with generations as elements
  #find the first object in env (usually only one)
  obj_name <- ls(env)[1]
  gen_list <- env[[obj_name]]
  
  #extract prey_details data frames from each generation and row bind
  df <- purrr::map_dfr(gen_list, ~ .x)
  return(df)
}

files <- list.files(path = "~/hdrive/GitHub/ballistic-movement/prey_results", pattern = "^[0-9]+g_prey_details\\.Rda$", full.names = TRUE)

all_prey_details <- map_dfr(files, load.prey.details)

summary <- all_prey_details %>% 
  filter(generation > max(generation - 100)) %>% 
  group_by(mass) %>% 
  summarise(    mean_lv = mean(lv, na.rm = TRUE),
                sd_lv = sd(lv, na.rm = TRUE),
                n = n(),
                se_lv = sd_lv / sqrt(n))

ggplot() +
  geom_point(data = summary, aes(x = mass, y = mean_lv, size = se_lv)) +
  scale_size_continuous(range = c(2, 8)) +
  labs(x = "Prey body mass (g)", y = "Ballistic lengthscale (m)", title = "Mean ballistic lengthscale (m) in final 100 generations by body mass") +
  theme.qel()

ggsave(file = "prey_results/figures/mean_lv_by_mass.png", width = 6.86, height = 3.5, units = "in",
       dpi = 600)
