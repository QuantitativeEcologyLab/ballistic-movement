
# habitat figures for vODD

library(spatstat)

source("simulation_scripts/01-prey-functions.R")

mass_prey <- 105500

#comparing density
pp1 <- makeHabitat(mass_prey,
                   r = 1,
                   mu = 1,
                   target_n = 2000,
                   cal = 1)

pp2 <- makeHabitat(mass_prey,
                   r = 1,
                   mu = 1,
                   target_n = 50000,
                   cal = 1)


p1 <- as.data.frame(pp1)
p2 <- as.data.frame(pp2)

plot1 <- 
  ggplot(p1, aes(x = x, y = y)) +
  geom_point(size = 0.5, color = "#95757cff") +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  theme_void() +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        plot.margin = margin(5,5,5,5))

plot2 <- 
  ggplot(p2, aes(x = x, y = y)) +
  geom_point(size = 0.5, color = "#95757cff") +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  theme_void() +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        plot.margin = margin(5,5,5,5))

FIG <- grid.arrange(plot1, 
                    plot2, nrow = 1)

ggsave(FIG, file = "figures/maintext/method_fid_files/density.png", width = 25, height = 6, units = "in", bg = "white")

# comparing clustering
pp1 <- makeHabitat(mass_prey,
                   r = 1,
                   mu = 1,
                   target_n = 2000,
                   cal = 1)

pp2 <- makeHabitat(mass_prey,
                   r = 1,
                   mu = 200,
                   target_n = 2000,
                   cal = 1)

p1 <- as.data.frame(pp1)
p2 <- as.data.frame(pp2)

plot1 <- 
  ggplot(p1, aes(x = x, y = y)) +
  geom_raster(data = as.data.frame(density(pp1)), aes(x = x, y = y, fill = value)) +
  geom_point(size = 1, color = "black") +
  geom_point(size = 0.6, color = "white") +
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
  geom_point(size = 1, color = "black") +
  geom_point(size = 0.6, color = "white") +
  scale_fill_viridis_c(option = "magma") +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  theme_void() +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        plot.margin = margin(5,5,5,5))

FIG <- grid.arrange(plot1, 
                    plot2, nrow = 1)

ggsave(FIG, file = "figures/maintext/method_fid_files/clustering.png", width = 25, height = 6, units = "in", bg = "white")

#comparing caloric variation
pp1 <- makeHabitat(mass_prey,
                   r = 1,
                   mu = 1,
                   target_n = 1000,
                   cal = 309)

pp2 <- makeHabitat(mass_prey,
                   r = 1,
                   mu = 1,
                   target_n = 1000,
                   cal = 309,
                   var = 50)

p1 <- as.data.frame(pp1)
p2 <- as.data.frame(pp2)


plot1 <- 
  ggplot(p1, aes(x = x, y = y, color = marks)) +
  geom_point(size = 1) +
  scale_colour_scico(palette = "lipari") +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  theme_void() +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        plot.margin = margin(5,5,5,5))

plot2 <- 
  ggplot(p2, aes(x = x, y = y, color = marks)) +
  geom_point(size = 1) +
  scale_colour_scico(palette = "lipari") +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  theme_void() +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        plot.margin = margin(5,5,5,5))

FIG <- grid.arrange(plot1, 
                    plot2, nrow = 1)

ggsave(FIG, file = "figures/maintext/method_fid_files/caloricvariance.png", width = 25, height = 6, units = "in", bg = "white")