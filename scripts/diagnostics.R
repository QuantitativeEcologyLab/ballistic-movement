# preamble

# Set the working directory
setwd("~/hdrive/GitHub/ballistic-movement")

# Import necessary packages
library(ggplot2)
library(dplyr)
library(gridExtra)
library(patchwork)

# Source the functions (ensure 'functions.R' is available in the working directory)
source("scripts/functions.R")
source("scripts/figure_functions.R")

###-----------------------------------------------------------------------------
#### prey diagnostics ----
###-----------------------------------------------------------------------------

# load data
load("simulations/prey_results/105500g_prey_res.Rda")
load("simulations/prey_results/105500g_prey_details.Rda")

# make data sets compatiable with ggplot2
prey_res_df <- do.call(rbind, prey_res)
prey_details_df <- do.call(rbind, prey_details)

# summary of traits over generation
prey.gen.summary(prey_res_df, prey_details_df)

# comparing variables to each other
prey.var.summary(prey_details_df)

#plot raster with track
#create data frame from food raster
df_raster <- as.data.frame(FOOD, xy = TRUE)
colnames(df_raster) <- c("x", "y", "calories")

#collect movement tracks from all individuals
prey_list <- lapply(seq_along(PREY_tracks), function(i){
  df <- as.data.frame(PREY_tracks[[i]])
  colnames(df)[2:3] <- c("x", "y")
  return(df)
})

#get x and y lines to show the patches in the raster
xlines <- unique(df_raster$x)
ylines <- unique(df_raster$y)

#custome function used to draw circles with ggplot
circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

#calculate the 95% (HR) and 99.9% (EXT) home range area
HR <- round(sqrt((-2*log(0.05)*pi)*prey.SIG(mass)))
EXT <- round(sqrt((-2*log(0.0001)*pi)* prey.SIG(mass)))

#create circles to add to ggplot
HR_area <- circleFun(diameter = 2*HR)
EXT_area <- circleFun(diameter = 2*EXT)

#plot it all together
ggplot() +
  geom_raster(data = df_raster, aes(x = x, y = y, fill = calories)) +
  geom_vline(xintercept = xlines, color = "white", alpha = 0.5) +
  geom_hline(yintercept = ylines, color = "white", alpha = 0.5) +
  scale_fill_viridis_c() +
  geom_path(data = PREY_tracks[[1]], aes(x = x, y = y), color = "black", linewidth = 0.7, alpha = 0.8) +
  geom_path(dat = HR_area, aes(x,y), color = "#467378") +
  geom_path(dat = EXT_area, aes(x,y), color = "#68855C") + 
  coord_equal() +
  theme_minimal() 

#save the figure
#ggsave(file = 'sim_results/july16/figures/30000g_1000thlifespan_movepath.png', width = 8, height = 8, dpi = 900, bg = "white")

###-----------------------------------------------------------------------------
#### predator diagnostics ----
###-----------------------------------------------------------------------------

# load data
load("simulations/pred_results/100000g_pred_res_stationary_prey.Rda")
load("simulations/pred_results/100000g_pred_details_stationary_prey.Rda")

# make data sets compatible with ggplot2
pred_res_df <- do.call(rbind, pred_res)
pred_details_df <- do.call(rbind, pred_details)

# summary of traits over generations
pred.gen.summary(pred_res_df, pred_details_df)

# comparing variables to each other
pred.var.summary(pred_details_df)

# plotting movement of predators and stationary prey

pred_tracks <- as.data.frame(PRED_tracks)
centres <- do.call(rbind, CENTRES)

ggplot() +
  geom_path(data = pred_tracks, aes(x = x, y = y), color = "steelblue", linewidth = 0.7, alpha = 0.8) +
  geom_point(data = centres, aes(x = x, y = y)) +
  coord_equal() +
  # xlim(-600,200) +
  # ylim(-600,200) +
  theme_minimal() 

pred_details_df <- do.call(rbind, pred_details)

prey_details_df <- do.call(rbind, prey_details)

#-------------------------------------------------------------------------------
# project updates ----
#-------------------------------------------------------------------------------

#load in your data
load('simulations/prey_results/200000g_prey_res_lowerlv_mutation.Rda')
load('simulations/prey_results/200000g_prey_details_lowerlv_mutation.Rda')

# make data sets compatible with ggplot2
prey_res_df <- do.call(rbind, prey_res)
prey_details_df <- do.call(rbind, prey_details)

prey.proj.update(prey_res_df, prey_details_df)

# save figures 
ggsave(FIG1, file = "project_updates/july/63500g_genplots.png", width = 6, height = 5, units = "in", dpi = 600)
ggsave(FIG2, file = "project_updates/july/63500g_variableplots.png", width = 6, height = 5, units = "in", dpi = 600)
