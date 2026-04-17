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
time_total <- lifespan # 1/500 of a lifespan

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

feed <- grazing(mass, track, food)

consumed <- attr(feed, "consumed")

consumed <- data.frame(id = rownames(consumed), consumed = consumed)

consumed.true <- consumed[consumed$consumed == TRUE, ]

food_df <- data.frame(x = food$x, y = food$y, consumed = consumed$consumed)

track_df <- data.frame(track)

#p1 <-
  ggplot() +
  geom_point(data = food_df, aes(x = x, y = y, colour = consumed), size = 1.3, alpha = 0.6, stroke = NA) +
  scale_color_manual(values = c("TRUE" = "#e18297", "FALSE" = "#2a3b2b")) +
  geom_path(data = track_df, aes(x=x, y=y), color = "black", linewidth = 0.4, alpha = 0.9) +
  labs(color = "Consumed") +
  # xlim(-1500,4000) +
  # ylim(-5000,1000) +
  coord_equal() +
  theme_void() +
  theme(legend.position = "none")


ggsave(p1, file = "figures/maintext/vODD/components/movetrack.png", width = 3, height = 3, units = "in", bg = "transparent", dpi = 600)

saveRDS(track, file = "figures/maintext/method_fig_files/track.Rds")
saveRDS(food, file = "figures/maintext/method_fig_files/food.Rds")

