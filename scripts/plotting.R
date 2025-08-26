#preamble

# Set the working directory
setwd("~/hdrive/GitHub/ballistic-movement")

# Import necessary packages
library(ggplot2)
library(dplyr)
library(gridExtra)
library(patchwork)
library(purrr)
library(minpack.lm)
library(mgcv) # for gam model
library(scico)

# Source the functions (ensure 'functions.R' is available in the working directory)
source("scripts/functions.R")

#----------------------------------------------------------------------
# make figures for results section----
#----------------------------------------------------------------------

#load in your data
load('simulations/prey_results/21500g_prey_res.Rda')
load('simulations/prey_results/21500g_prey_details.Rda')

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

#use palette "berlin" for small body sizes, "lipari" for large body sizes

#105500 grams
#p1 <-
  ggplot() +
  ggtitle("A") +
  geom_point(dat = prey_details_df, aes(x = generation, y = lv), col = "grey70", alpha = 0.8, size = 0.1) +
  geom_line(dat = prey_summary, aes(x = generation, y = lv), col = "#000", linewidth = 0.5) +
  labs(y = expression(bold(l[v])), x = "Generation") +
  theme.qel()

ggsave(p1, width = 10, height = 5, units = "in", dpi = 900, file = "figures/maintext/105kg_lv_gen.png")

#p2 <-
  ggplot(prey_details_df, aes(x = lv, y = cal_net, color = generation)) +
  ggtitle("B") +
  geom_point(size = 0.5) +
  scale_colour_scico(palette = 'lipari') +
  labs(y = "Net calories", x = "Ballistic lengthscale (m)") +
  theme.qel()

#p3 <- 
  ggplot(prey_details_df, aes(x = lv, y = speed, color = generation)) +
  ggtitle("C") +
  geom_point(size = 0.5) +
  scale_colour_scico(palette = "lipari") +
  labs(x = "Ballistic lengthscale (m)", y = "Speed (m/s)") +
  theme.qel()

#200000 grams
p4 <-
  ggplot() +
  ggtitle("A") +
  geom_point(dat = prey_details_df, aes(x = generation, y = lv), col = "#c7c9d1", alpha = 0.1, size = 0.1) +
  geom_line(dat = prey_summary, aes(x = generation, y = lv), col = "#24262d", linewidth = 0.5) +
  labs(y = expression(bold(l[v])), x = "Generation") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=9, family = "sans", face = "bold"),
        axis.title.x = element_text(size=9, family = "sans", face = "bold"),
        axis.text.y = element_text(size=8, family = "sans"),
        axis.text.x  = element_text(size=8, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  theme(
    legend.position = "none",
    panel.grid = element_blank())

p5 <-
  ggplot(prey_details_df, aes(x = lv, y = cal_net, color = generation)) +
  ggtitle("B") +
  geom_point(size = 0.5) +
  scale_colour_scico(palette = 'lipari') +
  labs(x = expression(bold(l[v])), y = "Net calories") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=9, family = "sans", face = "bold"),
        axis.title.x = element_text(size=9, family = "sans", face = "bold"),
        axis.text.y = element_text(size=8, family = "sans"),
        axis.text.x  = element_text(size=8, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  theme(
        legend.text = element_text(size = 6, family = "sans"),
        legend.title = element_text(size = 6, family = "sans", face = "bold"),
        legend.key.size = unit(0.3, "cm"),
        legend.spacing.y = unit(0.01, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent")) +
  theme(legend.position = c(0.85, 0.15)) 

p6 <- 
  ggplot(prey_details_df, aes(x = lv, y = speed, color = generation)) +
  ggtitle("C") +
  geom_point(size = 0.5) +
  scale_colour_scico(palette = "lipari") +
  labs(x = expression(bold(l[v])), y = "Speed (m/s)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=9, family = "sans", face = "bold"),
        axis.title.x = element_text(size=9, family = "sans", face = "bold"),
        axis.text.y = element_text(size=8, family = "sans"),
        axis.text.x  = element_text(size=8, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  theme(
    legend.text = element_text(size = 6, family = "sans"),
    legend.title = element_text(size = 6, family = "sans", face = "bold"),
    legend.key.size = unit(0.3, "cm"),
    legend.spacing.y = unit(0.01, "cm"),
    legend.margin = margin(0,0,0,0),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent")) +
  theme(legend.position = c(0.85, 0.85))

FIG2 <- grid.arrange(p4, p5, p6, ncol = 3)

ggsave(FIG2, width = 12, height = 4, units = "in", dpi = 900, bg = "white", file="figures/maintext/200kg_poster.png")

# restrict to final generations

lastgens <- prey_details_df %>%
  filter(generation >= max(generation) - 199) %>% 
  group_by(generation) %>% 
  mutate(mean_lv = mean(lv))

p4 <-
  ggplot() +
  ggtitle("D") +
  geom_point(dat = lastgens, aes(x = generation, y = lv), col = "grey70", alpha = 0.8, size = 0.1) +
  geom_line(dat = lastgens, aes(x = generation, y = mean_lv), col = "#000", linewidth = 0.5) +
  labs(y = "Ballistic lengthscale (m)", x = "Generation") +
  theme.qel()

p5 <-
  ggplot(lastgens, aes(x = lv, y = patches)) +
  ggtitle("E") +
  geom_point(color = "#000") +
  labs(y = "Patches visited", x = "Ballistic lengthscale (m)") +
  theme.qel()

p6 <- 
  ggplot(lastgens, aes(x = speed, y = lv)) +
  ggtitle("F") +
  geom_point(color = "#000") +
  labs(y = "Ballistic lengthscale (m)", x = "Speed (m/s)") +
  theme.qel()

FIG2 <- grid.arrange(p4, p5, p6, nrow = 1, ncol = 3)

FIG.final <- grid.arrange(FIG, FIG2, ncol = 1, nrow = 2)

# ------------------------------------------------------------------------------
## analysis across mass spectrum ----
# ------------------------------------------------------------------------------
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

files <- list.files(path = "~/hdrive/GitHub/ballistic-movement/simulations/prey_results", pattern = "^[0-9]+g_prey_details\\.Rda$", full.names = TRUE)

all_prey_details <- map_dfr(files, load.prey.details)

summary <- all_prey_details %>% 
  filter(generation > max(generation - 100)) %>% 
  group_by(mass) %>% 
  summarise(    mean_lv = mean(lv, na.rm = TRUE),
                sd_lv = sd(lv, na.rm = TRUE),
                n = n(),
                se_lv = sd_lv / sqrt(n))

# summary500g <- prey_details_df %>% 
#   filter(generation > max(generation - 100)) %>% 
#   group_by(mass) %>% 
#   summarise(mean_lv = mean(lv, na.rm = TRUE),
#             sd_lv = sd(lv, na.rm = TRUE),
#             n = n(),
#             se_lv = sd_lv / sqrt(n))
# 
# summary_comb <- bind_rows(summary, summary500g)


FIT_linear_log <- gam(log10(mean_lv) ~ log10(mass), family = tw(), data = summary) 

FIT_log <- gam(log10(mean_lv) ~ s(log10(mass), k = 3), family = tw(), data = summary) #lowest AIC

FIT_linear <- gam((mean_lv) ~ (mass), family = tw(), data = summary) 

FIT <- gam((mean_lv) ~ s((mass), k = 3), family = tw(), data = summary)

AIC(FIT, FIT_linear, FIT_log, FIT_linear_log)
summary(FIT_log)

ggplot(data = summary, aes(x = mass/1000, y = mean_lv)) +
  geom_smooth(method = gam, formula = y ~ s(log10(x), k = 4), method.args = list(family = "tw"), color = "#5b557bff") +
  geom_point(aes(size = se_lv), alpha = 1, color = "grey20") +
  scale_size_continuous(range = c(2, 8)) +
  labs(x = "Prey body mass (kg)", y = expression(bold(l[v]))) +
  coord_cartesian(ylim = c(-5,250)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=9, family = "sans", face = "bold"),
        axis.title.x = element_text(size=9, family = "sans", face = "bold"),
        axis.text.y = element_text(size=8, family = "sans"),
        axis.text.x  = element_text(size=8, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  theme(legend.position = "none",
        panel.grid = element_blank())


ggsave(file = "figures/maintext/massspectrum_model.png", width = 6.86, height = 3, units = "in",
       dpi = 600)
