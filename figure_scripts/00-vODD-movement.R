setwd("H:/GitHub/ballistic-movement/")

library(tidyverse)
library(ctmm)
library(scico)
library(paletteer)
library(ggforce)

source("simulation_scripts/01-prey-functions.R")
source("figure_scripts/01-custom-ggplot-theme.R")

#...............................................................................
## with consumed patches as a different colour

mass <- 105500

tau_p <- prey.tau_p(mass)
tau_v <- prey.tau_v(mass)
sig <- prey.SIG(mass)

mod <- ctmm(tau = c(tau_p, tau_v), mu = c(0,0), sigma = sig)

t <- sampling(mass, x = 10)

food <- makeHabitat(mass, 
                    r = 1, 
                    mu = 1,
                    target_n = 2000, 
                    cal = 1013,
                    var = 10)

track <- simulate(mod, t = t)

feed <- grazing(mass, track, food)

consumed <- attr(feed, "consumed")

consumed <- data.frame(id = rownames(consumed), consumed = consumed)

consumed.true <- consumed[consumed$consumed == TRUE, ]

food_df <- data.frame(x = food$x, y = food$y, consumed = consumed$consumed)

track_df <- data.frame(track)

pr_df <- data.frame(
  x0 = track_df$x[14550],
  y0 = track_df$y[14550],
  pr = sqrt(prey.SIG(mass))*0.05
)

circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

p1 <- 
  ggplot() +
  geom_point(data = food_df, aes(x = x, y = y, colour = consumed), size = 0.8) +
  scale_color_manual(values = c("TRUE" = "#e18297", "FALSE" = "#2a3b2b")) +
  geom_path(data = track_df, aes(x=x, y=y), color = "grey25", linewidth = 0.3) +
 # geom_path(data = circleFun(center = c(pr_df$x0, pr_df$y0), diameter = pr_df$pr), aes(x=x,y=y),color = "green", linewidth = 1) +
  labs(color = "Consumed") +
  xlim(-1500,4000) +
  ylim(-5000,1000) +
  coord_equal() +
  theme_void() +
  theme(legend.position = "none")


ggsave(p1, file = "project_updates/26-02feb/methods.png", width = 3, height = 3, units = "in", bg = "transparent", dpi = 600)

saveRDS(track, file = "figures/maintext/method_fig_files/track.Rds")
saveRDS(food, file = "figures/maintext/method_fig_files/food.Rds")

track <- readRDS("figures/maintext/method_fig_files/track.Rds")
food <- readRDS("figures/maintext/method_fig_files/food.Rds")
