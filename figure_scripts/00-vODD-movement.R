setwd("H:/GitHub/ballistic-movement/")

library(tidyverse)
library(ctmm)
library(scico)
library(paletteer)
library(ggforce)

source("simulation_scripts/01-prey-functions.R")

track <- readRDS("simulations/figures/track.Rds")
food <- readRDS("simulations/figures/food.Rds")
#...............................................................................
## with consumed patches as a different colour
#...............................................................................

mass <- 105500

tau_p <- prey.tau_p(mass, variance = TRUE)
tau_v <- prey.tau_v(mass, variance = TRUE)
sig <- prey.SIG(mass)

mod <- ctmm(tau = c(tau_p, tau_v), mu = c(0,0), sigma = sig)

lifespan <- (4.88*mass^0.153) * 31536000 # years to seconds
time_total <- lifespan * 0.002 # 1/500 of a lifespan

#sampling interval (tau_v) in seconds, max prevents tau_v < 1
#increasing x decreases interval, making sampling more frequent
interval <- max(1, round(prey.tau_v(mass))) / 10

#lifespan and sampling interval for simulations
t <- seq(0,
         time_total,
         interval)

food <- makeHabitat(mass, 
                    r = 1, 
                    mu = 1,
                    target_n = 500, 
                    cal = 4000,
                    var = 0)

track <- simulate(mod, t = t)

## start here if using pre-saved landscape and track

feed <- grazing(mass, track, food)

consumed <- attr(feed, "consumed")

consumed <- data.frame(id = rownames(consumed), consumed = consumed)

consumed.true <- consumed[consumed$consumed == TRUE, ]

food_df <- data.frame(x = food$x, y = food$y, consumed = consumed$consumed)

track_df <- data.frame(track)

p1 <-
  ggplot() +
    geom_path(data = track_df, aes(x = x, y = y), color = "black", linewidth = 0.4, alpha = 0.8) +
    geom_point(data = food_df, aes(x = x, y = y, colour = consumed), size = 1.3, alpha = 1, stroke = NA) +
    scale_color_manual(values = c("TRUE" = "#e18297", "FALSE" = "#2a3b2b")) +
    labs(color = "Consumed") +
    xlim(-1500,4000) +
    ylim(-4500,1000) +
    coord_equal() +
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


ggsave(p1, file = "figures/ODD/foraging.png", width = 3, height = 3, units = "in", bg = "transparent", dpi = 600)

saveRDS(track, file = "figures/maintext/method_fig_files/track.Rds")
saveRDS(food, file = "figures/maintext/method_fig_files/food.Rds")

