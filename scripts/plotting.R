# preamble
# this code was used to generate the main text figures for the project
# "A computational approach to identifying evolutionarily stable strategies in mammalian search behaviour"

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
library(cowplot) # for combining figures

# Source the functions (ensure 'functions.R' is available in the working directory)
source("scripts/functions.R")

#----------------------------------------------------------------------
# main figures----
#----------------------------------------------------------------------

## load in  data and convert to data frame
load('simulations/prey_results/200000g_prey_details.Rda')
prey_details_df <- do.call(rbind, prey_details)

## calculate mean of each parameters for every generation
prey_summary <- prey_details_df %>%
  group_by(generation) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop")

## lv versus generation
p1 <-
  ggplot() +
  ggtitle("A") +
  geom_point(dat = prey_details_df, aes(x = generation, y = lv), col = "#c7c9d1", alpha = 0.1, size = 0.1) +
  geom_line(dat = prey_summary, aes(x = generation, y = lv), col = "#24262d", linewidth = 0.5) +
  labs(y = "Ballistic length scale (m)", x = "Generation") +
  ylim(0, 400) +
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

## lv versus net calories
p2 <-
  ggplot(prey_details_df, aes(x = lv, y = cal_net, color = generation)) +
  ggtitle("B") +
  geom_point(size = 0.5) +
  scale_colour_scico(palette = 'lipari') +
  labs(x = "Ballistic length scale (m)", y = "Net calories") +
  ylim(-250000,1416405) +
  xlim(0, 400) +
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

## lv versus speed
p3 <- 
  ggplot(prey_details_df, aes(x = lv, y = speed, color = generation)) +
  ggtitle("C") +
  geom_point(size = 0.5) +
  scale_colour_scico(palette = "lipari") +
  labs(x = "Ballistic length scale (m)", y = "Speed (m/s)") +
  ylim(0, 20) +
  xlim(0, 400) +
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

## combine figure
FIG2 <- grid.arrange(p1, p2, p3, ncol = 3)

## save figure
ggsave(FIG2, width = 12, height = 4, units = "in", dpi = 900, bg = "white", file="figures/maintext/200kg_poster.png")

# ------------------------------------------------------------------------------
## analysis across mass spectrum ----
# ------------------------------------------------------------------------------

# create function which loads all files from directory

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

# load files 

## assign what files to load
files <- list.files(path = "~/hdrive/GitHub/ballistic-movement/simulations/prey_results", pattern = "^[0-9]+g_prey_details\\.Rda$", full.names = TRUE)
## load all files to single data frame
all_prey_details <- map_dfr(files, load.prey.details)

## data wrangling
summary <- all_prey_details %>% 
  filter(generation > max(generation - 100)) %>% 
  group_by(mass) %>% 
  summarise(    mean_lv = mean(lv, na.rm = TRUE),
                sd_lv = sd(lv, na.rm = TRUE),
                n = n(),
                se_lv = sd_lv / sqrt(n))

# adding simulations that are partially completed

# summary500g <- prey_details_df %>% 
#   filter(generation > max(generation - 100)) %>% 
#   group_by(mass) %>% 
#   summarise(mean_lv = mean(lv, na.rm = TRUE),
#             sd_lv = sd(lv, na.rm = TRUE),
#             n = n(),
#             se_lv = sd_lv / sqrt(n))
# 
# summary_comb <- bind_rows(summary, summary500g)

# TESTING DIFFERENT GAM MODELS

FIT_linear_log <- gam(log10(mean_lv) ~ log10(mass), family = tw(), data = summary) 

FIT_log <- gam(log10(mean_lv) ~ s(log10(mass), k = 3), family = tw(), data = summary) #lowest AIC

FIT_linear <- gam((mean_lv) ~ (mass), family = tw(), data = summary) 

FIT <- gam((mean_lv) ~ s((mass), k = 3), family = tw(), data = summary)

AIC(FIT, FIT_linear, FIT_log, FIT_linear_log)
summary(FIT_log)

# create figure 

## remove 21.5kg because it didn't stabilize
summary <- summary %>% 
  filter(mass > 21500)

## add colour column for highlighting mass sizes
summary <- summary %>% 
  mutate(label = case_when(
    mass == 42500 ~ "small",
    mass == 200000 ~ "big",
    TRUE ~ "other"
  ))

## generate plot
mass_spec <- 
  ggplot(data = summary, aes(x = mass/1000, y = mean_lv)) +
  geom_smooth(method = gam, formula = y ~ s(log10(x), k = 4), method.args = list(family = "tw"), color = "grey20") +
  geom_point(aes(size = se_lv, color = label), alpha = 1) +
  scale_color_manual(values = c("small" = "#5d6da7",
                                "big" = "#ff6d6b",
                                "other" = "grey20")) + 
  scale_size_continuous(range = c(2, 8)) +
  labs(x = "Prey body mass (kg)", y = "Ballistic length scale (m)") +
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
        plot.margin = unit(c(1,1,1,1), "cm")) +
  theme(legend.position = "none",
        panel.grid = element_blank())

## add label A
mass_spec <- ggdraw(mass_spec) +
  draw_plot_label("A", x = 0.02, y = 0.99, size = 14)

## save figure without insets
ggsave(mass_spec, 
       file = "figures/maintext/massspectrum_model.png", 
       width = 6.86, height = 3, units = "in",
       dpi = 600)


# create insets

## smallest mass
mass1 <- 42500
t <- sampling(mass1)

## load in data
load('simulations/prey_results/42500g_prey_details.Rda')
## convert to data frame
prey_details_df <- do.call(rbind, prey_details)
## data wrangling
prey_summary <- prey_details_df %>%
  group_by(generation) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop") %>% 
  filter(generation == 1000)

## extract movement parameters
tau_v <- prey_summary$tau_v
tau_p <- prey_summary$tau_p
sig <- prey_summary$sig

## create ctmm model
mod <- ctmm(tau = c(tau_p, tau_v), mu = c(0,0), sigma = sig)

## simulate movement tracks
tracks <- simulate(mod, t = t)

## convert to data frame
df <- as.data.frame(tracks)

## take subset of data to generate plot (increase visibility)
split_point1 <- ceiling(nrow(df)/25)
df <- df[1:split_point1, ]

## save movement tracks 
save(df, file = "figures/maintext/42kg_movetracks.Rda")

## create first inset
a <- 
  ggplot(df, aes(x = x, y = y)) +
  geom_path(linewidth = 0.3, col = "#5d6da7") +
  ggplot2::annotate("text", x = max(df$x), y = min(df$y), label = "42.5 kg", hjust = 1.1, vjust = 0.3, size = 4, family = "sans") +
  ggplot2::annotate("text", x = min(df$x), y = max(df$y), label = "(i)", hjust = 0, vjust = 1, size = 4, fontface = "bold", family = "sans") +
  coord_equal() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  theme(legend.position = "none",
        panel.grid = element_blank())

## largest mass

mass2 <- 200000
t2 <- sampling(mass2)

## load in data
load('simulations/prey_results/200000g_prey_details.Rda')
## convert to data frame
prey_details_df2 <- do.call(rbind, prey_details)
## data wrangling
prey_summary2 <- prey_details_df2 %>%
  group_by(generation) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop") %>% 
  filter(generation == 1000)

## extract movement parameters
tau_v2 <- prey_summary2$tau_v
tau_p2 <- prey_summary2$tau_p
sig2 <- prey_summary2$sig

## create ctmm model
mod2 <- ctmm(tau = c(tau_p2, tau_v2), mu = c(0,0), sigma = sig2)

## simulate movement
tracks2 <- simulate(mod2, t = t2)

## convert to data frame
df2 <- as.data.frame(tracks2)

## take subset of data to generate plot (increase visibility)
split_point2 <- ceiling(nrow(df2)/25)
df2 <- df2[1:split_point2, ]

## save movement tracks
save(df2, file = "figures/maintext/200kg_movetracks.Rda")

b <- 
  ggplot(df2, aes(x = x, y = y)) +
  geom_path(linewidth = 0.5, col = "#ff6d6b") +
  ggplot2::annotate("text", x = max(df2$x), y = min(df2$y), label = "200 kg", hjust = 1.1, vjust = 0.3, size = 4, family = "sans") +
  ggplot2::annotate("text", x = min(df2$x), y = max(df2$y), label = "(ii)", hjust = 0, vjust = 1, size = 4, fontface = "bold", family = "sans") +
  theme_void() +
  coord_equal() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  theme(legend.position = "none",
        panel.grid = element_blank())

## make both inserts have the same aspect ratio
a <- a + theme(aspect.ratio = 1)
b <- b + theme(aspect.ratio = 1)

## combine insets
inset <- 
  wrap_plots(a, b, ncol = 2, nrow = 1, heights = c(1, 1), widths = c(1,1))

## add insets to model figure
mass_spec + 
  inset_element(inset, left = 0.09, right = 0.55, bottom = 0.3, top = 1.07, align_to = "panel") 

## save figure
ggsave(file = "figures/maintext/mass_spec_insets.png", width = 10, height = 5.3, dpi = 900)

