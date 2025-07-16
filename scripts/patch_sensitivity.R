
#preamble ----

# this script contains the code needed to investigate the sensitivity of patch width 
# and sampling interval on the number of patches that an individual encounters.

#set the working directory
setwd("~/hdrive/GitHub/ballistic-movement")

#import necessary packages
library(ctmm)
library(terra)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(patchwork)
library(data.table)
library(parallel)

#source the functions (ensure 'functions.R' is available in the working directory)
source("scripts/functions.R")

#--------------------------------------------------------------------------
# create function for sensitivity analysis  -------------------------------
#--------------------------------------------------------------------------

#for adding circles for 95% and 99% HR area
circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

#for generating movement tracks for patch per 95% HR area analysis
get.tracks <- function(mass_prey, k) {
  
  t <- sampling(mass_prey, x = 1)
  
  interval <- attr(t, "interval")
  
  # Create or load FOOD raster based on parameters
  FOOD <- createFoodRaster(mass_prey, k,
                           calories = 1)
  
  prey_tau_p <- prey.tau_p(mass_prey, variance = TRUE)
  prey_tau_v <- prey.tau_v(mass_prey, variance = TRUE)
  prey_sig <- prey.SIG(mass_prey)
  prey_lv <- sqrt((prey_tau_v/prey_tau_p) * prey_sig)
        
  PREY_mods <- ctmm(tau = c(prey_tau_p, prey_tau_v),
                            mu = c(0,0),
                            sigma = prey_sig)
    

  PREY_tracks <- simulate(PREY_mods, t = t)
    
  benefits_prey <- grazing(PREY_tracks, FOOD)
  
  patches <- attr(benefits_prey, "patches")

  prey <- data.frame(patches = unlist(patches),
                       k = k,
                       interval = interval)
  
  return(list(prey_details = prey,
              tracks = PREY_tracks,
              food = FOOD))
    
  } # closes over function

#for generating a single, high resolution movement track for
#sampling interval sensitivity analysis
get.interval.tracks <- function(mass_prey, x, seed, mod) {
  
  t <- sampling(mass_prey, x = x)
  
  set.seed(seed)
  track <- simulate(mod, t = t)
  
  attr(track, "seed") <- seed
  
  return(track)
}

#create function to thin the track for sampling interval sensitivity analysis
thin <- function(track, food, t_thin){
  df <- as.data.table(track)
  
  setkey(df, t)
  
  df_thin <- data.table(t = t_thin)
  
  df_thin <- df[df_thin, roll = "nearest"]
  
  pts <- data.frame(x = df_thin$x, y = df_thin$y)
  cell_ids <- suppressWarnings(cellFromXY(food, pts))
  
  patches <- sum(c(TRUE, diff(cell_ids) != 0), na.rm = TRUE)
  
  return(list(
    summary = data.frame(
      interval = attr(t_thin, "interval"),
      patches = patches
    ),
    data = df_thin
  ))
}

#create function to generate tracks and calculate offspring for kcal/patch analysis
get.kcal.tracks <- function(mass_prey, k, cal) {
  
  t <- sampling(mass_prey, x = 35)
  
  interval <- attr(t, "interval")
  
  # Create or load FOOD raster based on parameters
  FOOD <- createFoodRaster(mass_prey, k,
                           calories = cal)
  
  prey_tau_p <- prey.tau_p(mass_prey, variance = TRUE)
  prey_tau_v <- prey.tau_v(mass_prey, variance = TRUE)
  prey_sig <- prey.SIG(mass_prey)
  prey_lv <- sqrt((prey_tau_v/prey_tau_p) * prey_sig)
  
  PREY_mods <- ctmm(tau = c(prey_tau_p, prey_tau_v),
                    mu = c(0,0),
                    sigma = prey_sig)
  
  
  PREY_tracks <- simulate(PREY_mods, t = t)
  
  benefits_prey <- grazing(PREY_tracks, FOOD)
  
  patches <- attr(benefits_prey, "patches")
  
  prey_speed <- get.speed(models = PREY_mods)
  
  prey_gains <- prey.cals.net(IDs = benefits_prey,
                              mass = mass_prey,
                              speed = prey_speed,
                              t = t)
  
  offspring <- prey.fitness(mass = mass_prey,
                            cal_net = prey_gains)
  
  prey <- data.frame(patches = unlist(patches),
                     k = k,
                     interval = interval,
                     calories = cal,
                     offspring = offspring)
  
  return(list(prey_details = prey,
              tracks = PREY_tracks,
              food = FOOD)) 
}

run.track <- function(seed) {
  get.interval.tracks(mass_prey, x = 100, seed = seed, mod = mod)
}

#--------------------------------------------------------------------------
# section 1: patch resolution sensitivity ----
#--------------------------------------------------------------------------

#prey mass
mass_prey <- 100

#create grid of values to test
param_grid <- expand.grid(
  k = seq(1, 1000001, 500))

#create lists for storing results
patch_results <- list()
# tracks <- list()
# food <- list()

#begin sensitivity analysis
for(i in 1:nrow(param_grid)){
  set.seed(5) # set seed at beginning to make movement paths identical
  
  k <- param_grid$k[i] # assign the k value for each test

  #progress reports
  cat("\n--- Running Scenario", i, "of", nrow(param_grid), "---\n")
  cat(sprintf("k: %f\n", k))
  
  #run the simulation 
  res <- get.tracks(mass_prey = mass_prey,
                    k = k)
  
  #store results
  patch_results[[i]] <- res$prey_details
  # tracks[[i]] <- res$tracks
  # food[[i]] <- res$food
  
  #save results
  save(patch_results, file = "sim_results/sensitivity/supplementary/data/100g_vary_k.Rda")
}

#convert results to a data frame (analysis 1)
patch_results_df <- bind_rows(patch_results)
# save(patch_results_df, file = 'sim_results/sensitivity/supplementary/data/30000g_vary_patches_details.Rda')
# save(tracks, file = 'sim_results/sensitivity/supplementary/data/30000g_vary_patches_tracks.Rda')
# save(food, file = "sim_results/sensitivity/supplementary/data/30000g_vary_tracks_food.Rda")

#generate plot for number of patches visited versus number of patches in 95% HR area
patch_fig <- ggplot(patch_results_df, aes(x = k, y = patches)) +
  geom_point(size = 0.5, alpha = 0.7) +
  labs(y = "patches visited", x = 'patches in 95% HR') +
  theme.qel(legend = FALSE)
print(patch_fig)

# ggsave(patch_fig,
#        width = 6.86, height = 3.5, units = "in",
#        dpi = 600,
#        bg = "white",
#        file="~/hdrive/GitHub/ballistic-movement/figures/supporting_analysis/patchvisits_patchesperHR.png")

#plot raster with track (for all analyses)
df_raster <- as.data.frame(food[[2]], xy = TRUE)
colnames(df_raster) <- c("x", "y", "calories")

df <- as.data.frame(tracks[[2]])
colnames(df)[2:3] <- c("x", "y")

xlines <- unique(df_raster$x)
ylines <- unique(df_raster$y)

HR <- round(sqrt((-2*log(0.05)*pi)*prey.SIG(mass_prey)))
EXT <- round(sqrt((-2*log(0.0001)*pi)* prey.SIG(mass_prey)))

HR_area <- circleFun(diameter = 2*HR)
EXT_area <- circleFun(diameter = 2*EXT)

ggplot() +
  geom_raster(data = df_raster, aes(x = x, y = y, fill = calories)) +
  geom_vline(xintercept = xlines, color = "white", alpha = 0.5) +
  geom_hline(yintercept = ylines, color = "white", alpha = 0.5) +
  geom_path(data = df, aes(x = x, y = y), color = "black", linewidth = 0.7, alpha = 0.8) +
  geom_path(dat = HR_area, aes(x,y), color = "grey60") +
  geom_path(dat = EXT_area, aes(x,y), color = "grey90") +
  coord_equal() +
  theme_minimal() +
  theme(legend.position = "none") 

#--------------------------------------------------------------------------
# section 2: sampling interval sensitivity ----
#--------------------------------------------------------------------------

#one core ---- 
mass_prey <- 100

food <- createFoodRaster(mass_prey, k = 240000)

tau_p <- prey.tau_p(mass_prey)
tau_v <- prey.tau_v(mass_prey)
sig <- prey.SIG(mass_prey)

mod <- ctmm(tau = c(tau_p, tau_v),
            mu = c(0,0),
            sigma = sig)

seeds <- c(1, 2, 3)

tracks <- list()

for(i in seq_along(seeds)) {
  seed <- seeds[i]
  
  results <- get.interval.tracks(mass_prey, x = 500, seed, mod)
  
  tracks[[i]] <- results
}
#add seed column and convert to data frame before saving
tracks_seeds <- Map(function(df, seed) {
  df$seed <- seed
  df
}, tracks, seeds)

tracks_comb <- do.call(rbind, tracks_seeds)

save(tracks_comb, file = "sim_results/sensitivity/supplementary/data/5000000g_vary_interval_tracks.Rda")

load('sim_results/sensitivity/supplementary/data/100g_vary_interval_tracks.Rda')

tracks <- split(tracks_comb, tracks_comb$seed)

#thin the tracks
x_vals <- seq(100, 0.5, -0.5)

thinned <- list()
thin_tracks <- list()

for(j in seq_along(tracks)){
  track_j <- tracks[[j]]
  seed <- unique(track_j$seed)
  
  for(i in seq_along(x_vals)){
  x <- x_vals[i]
  
  t_thin <- sampling(mass_prey, x)
  
  result <- thin(track_j, food, t_thin)
  
  result$summary$x <- x
  result$summary$seed <- seed
  
  thinned[[length(thinned) + 1]] <- result$summary
  # thin_tracks[[length(thin_tracks) + 1]] <- result$data
  }
}

interval_res <- do.call(rbind, thinned)

#convert results to a data frame (analysis 2)
save(interval_res, file = 'sim_results/sensitivity/supplementary/data/100g_vary_interval_details.Rda')

#paralleled for 100 grams to increase speed ----
#x values from 0 to 300 instead of 500 due to run speeds/crashing 

Ncores <- 3

mass_prey <- 100

food <- createFoodRaster(mass_prey, k = 240000)

tau_p <- prey.tau_p(mass_prey)
tau_v <- prey.tau_v(mass_prey)
sig <- prey.SIG(mass_prey)

mod <- ctmm(tau = c(tau_p, tau_v),
            mu = c(0,0),
            sigma = sig)

seeds <- c(1, 2, 3)

tracks <- list()

tracks <- mclapply(seeds, run.track, mc.cores = Ncores)

tracks_seeds <- Map(function(df, seed) {
  df$seed <- seed
  df
}, tracks, seeds)

tracks_comb <- do.call(rbind, tracks_seeds)

save(tracks_comb, file = "sim_results/sensitivity/supplementary/data/100g_vary_interval_tracks.Rda")

# thin the track with the code for one core (above)

#--------------------------------------------------------------------------
# section 3: kcal per patch sensitivity ----
#--------------------------------------------------------------------------

#prey mass
mass_prey <- 30000

#create grid of values to test
param_grid <- expand.grid(
  cal = seq(0, 10, 0.5))

#set number of patches in 95% HR area
k <- 100000

#create lists for storing results
patch_results <- list()
tracks <- list()
food <- list()

#begin sensitivity analysis
for(i in 1:nrow(param_grid)){
  set.seed(5) # set seed at beginning to make movement paths identical
  
  cal <- param_grid$cal[i] # assign the k value for each test
  
  #progress reports
  cat("\n--- Running Scenario", i, "of", nrow(param_grid), "---\n")
  cat(sprintf("kcal/patch: %f\n", cal))
  
  #run the simulation 
  res <- get.kcal.tracks(mass_prey = mass_prey, k = k, cal = cal)
  
  #store results
  patch_results[[i]] <- res$prey_details
  tracks[[i]] <- res$tracks
  food[[i]] <- res$food
  
  #save results
  save(patch_results, file = "sim_results/sensitivity/supplementary/data/30000g_vary_kcal_per_patch.Rda")
}

#convert to data frame to use in ggplot
patch_results_df <- bind_rows(patch_results)

ggplot(patch_results_df, aes(x = calories, y = offspring)) +
  geom_point(size = 1.5, alpha = 0.7) +
  labs(y = "offspring", x = 'kcal per patch') +
  theme_minimal()

#confirm that movement tracks are identical
#plot raster with track (for all analyses)
df_raster <- as.data.frame(food[[5]], xy = TRUE)
colnames(df_raster) <- c("x", "y", "calories")

save(df_raster, file = "~/hdrive/GitHub/ballistic-movement/sim_results/sensitivity/supplementary/data/30000g_vary_kcal_per_patch_foodraster.Rda")

df <- as.data.frame(tracks[[5]])
colnames(df)[2:3] <- c("x", "y")
save(df, file = "~/hdrive/GitHub/ballistic-movement/sim_results/sensitivity/supplementary/data/30000g_vary_kcal_per_patch_tracks.Rda")

xlines <- unique(df_raster$x)
ylines <- unique(df_raster$y)

HR <- round(sqrt((-2*log(0.05)*pi)*prey.SIG(mass_prey)))
EXT <- round(sqrt((-2*log(0.0001)*pi)* prey.SIG(mass_prey)))

HR_area <- circleFun(diameter = 2*HR)
EXT_area <- circleFun(diameter = 2*EXT)

ggplot() +
  geom_raster(data = df_raster, aes(x = x, y = y, fill = calories)) +
  geom_vline(xintercept = xlines, color = "white", linewidth = 0.01, alpha = 0.5) +
  geom_hline(yintercept = ylines, color = "white", linewidth = 0.01, alpha = 0.5) +
  geom_path(data = df, aes(x = x, y = y), color = "black", linewidth = 0.7, alpha = 0.8) +
  geom_path(dat = HR_area, aes(x,y), color = "grey60") +
  geom_path(dat = EXT_area, aes(x,y), color = "grey90") +
  coord_equal() +
  theme_minimal() +
  theme(legend.position = "none") 
