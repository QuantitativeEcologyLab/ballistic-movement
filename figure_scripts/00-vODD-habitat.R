# habitat figures for vODD

#..............................................................................
# load libraries ----
#..............................................................................
library(spatstat)
library(tidyverse)
library(gridExtra)
library(patchwork)
library(viridis)

# load custom functions
source("simulation_scripts/01-prey-functions.R")
#..............................................................................

# set animal body mass
mass_prey <- 105500

#...............................................................................
# comparing density----
#...............................................................................

## sparse habitat
p1 <- as.data.frame(
  makeHabitat(mass_prey,
                   r = 1,
                   mu = 1,
                   target_n = 450,
                   cal = 1)
  )

## dense habitat
p2 <- as.data.frame(
  makeHabitat(mass_prey,
                   r = 1,
                   mu = 1,
                   target_n = 1000,
                   cal = 1)
  )

plot1 <- 
  ggplot(p1, aes(x = x, y = y)) +
  geom_point(size = 2, color = "#B63679FF") +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  theme_void() +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        plot.margin = margin(5,5,5,5))

plot2 <- 
  ggplot(p2, aes(x = x, y = y)) +
  geom_point(size = 2, color = "#B63679FF") +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  theme_void() +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        plot.margin = margin(5,5,5,5))

FIG <- grid.arrange(plot1, 
                    plot2, nrow = 1)

ggsave(FIG, file = "figures/maintext/vODD_figures/density.png", width = 10, height = 5, units = "in", bg = "white")

#...............................................................................
# comparing clustering----
#...............................................................................

# no clustering
pp1 <- makeHabitat(mass_prey,
                   r = 1,
                   mu = 1,
                   target_n = 500,
                   cal = 1)

#max clustering
pp2 <- makeHabitat(mass_prey,
                   r = 1,
                   mu = 200,
                   target_n = 500,
                   cal = 1)

p1 <- as.data.frame(pp1)
p2 <- as.data.frame(pp2)

plot1 <- 
  ggplot(p1, aes(x = x, y = y)) +
  geom_raster(data = as.data.frame(density(pp1)), aes(x = x, y = y, fill = value)) +
  geom_point(size = 2, color = "black") +
  geom_point(size = 1.6, color = "white") +
  scale_fill_viridis_c(option = "magma") +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  theme_void() +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        plot.margin = margin(5,5,5,5))

plot2 <-
  ggplot(p2, aes(x=x,y=y))+
  geom_raster(data = as.data.frame(density(pp2)), aes(x = x, y = y, fill = value)) +
  geom_point(size = 2, color = "black") +
  geom_point(size = 1.6, color = "white") +
  scale_fill_viridis_c(option = "magma") +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  theme_void() +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        plot.margin = margin(5,5,5,5))

FIG <- grid.arrange(plot1, 
                    plot2, nrow = 1)

ggsave(FIG, file = "figures/maintext/vODD_figures/clustering.png", width = 10, height = 5, units = "in", bg = "white")

#...............................................................................
#comparing caloric variation----
#...............................................................................

p1 <- as.data.frame(
  makeHabitat(mass_prey,
               r = 1,
               mu = 1,
               target_n = 500,
               cal = 4000)
  )

p2 <- as.data.frame(
  makeHabitat(mass_prey,
               r = 1,
               mu = 1,
               target_n = 500,
               cal = 4000,
               var = 50)
  )

plot1 <- 
  ggplot(p1, aes(x = x, y = y, color = marks)) +
  geom_point(size = 3) +
  scale_color_viridis_c(option = "magma") +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  theme_void() +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        plot.margin = margin(5,5,5,5))

plot2 <- 
  ggplot(p2, aes(x = x, y = y, color = marks)) +
  geom_point(size = 3) +
  scale_color_viridis_c(option = "magma") +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  labs(color = "Calories") +
  theme_void() +
  theme() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        plot.margin = margin(5,5,5,5))

FIG <- plot1 + plot2

ggsave(FIG, file = "figures/maintext/vODD_figures/caloricvariance.png", width = 10, height = 5, units = "in", bg = "white")
