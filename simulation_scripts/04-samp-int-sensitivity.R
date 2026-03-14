#preamble ----

# this script contains the code needed to investigate the sensitivity of patch width 
# and sampling interval on the number of patches that an individual encounters.

#set the working directory
setwd("~/hdrive/GitHub/ballistic-movement")

#import necessary packages
library(extraDistr) # for bivariate poisson distribution of centres
library(parallel) # for parallel computing
library(tictoc) # for timing simulations
library(tidyverse)
library(data.table)
library(paletteer)

#source the functions (ensure 'functions.R' is available in the working directory)
source("simulation_scripts/01-prey-functions.R")
source("figure_scripts/01-custom-ggplot-theme.R")

#..........................................................................
# create function for sensitivity analysis  -------------------------------
#..........................................................................

get.interval.tracks <- function(mass_prey, x, seed, mod) {
  
  t <- sampling(mass_prey, x = x)
  
  set.seed(seed)
  track <- simulate(mod, t = t)
  
  attr(track, "seed") <- seed
  
  return(track)
}

#create function to thin the track for sampling interval sensitivity analysis
thin <- function(mass, track, food, t_thin){
  df <- as.data.table(track)
  
  setkey(df, t)
  
  df_thin <- data.table(t = t_thin)
  
  df_thin <- df[df_thin, roll = "nearest"]
  
  pts <- data.frame(x = df_thin$x, y = df_thin$y)
  
  calories <- grazing(mass, pts, food)
  
  consumed <- attr(calories, "patches")
  
  return(list(
    summary = data.frame(
      interval = attr(t_thin, "interval"),
      patches = consumed
    ),
    data = df_thin
  ))
}

#set x to the highest resolution you want to test (big value means smaller interval)
#interval = tau_v / x
run.track <- function(seed) {
  get.interval.tracks(mass_prey, x = 100, seed = seed, mod = mod)
}

#..........................................................................
# sampling interval sensitivity  -------------------------------
#..........................................................................

#set number of cores
Ncores <- 1

#set mass of prey to test
mass_prey <- 105500

#create food landscape point process
food <- makeHabitat(mass_prey,
                    r = 1,
                    mu = 1, 
                    target_n = 50000,
                    cal = 1)

#set model for simulation
tau_p <- prey.tau_p(mass_prey)
tau_v <- prey.tau_v(mass_prey)
sig <- prey.SIG(mass_prey)

mod <- ctmm(tau = c(tau_p, tau_v),
            mu = c(0,0),
            sigma = sig)

#set seeds
seeds <- c(1, 2, 3)

#create empty list to store results
tracks <- list()

#create movement tracks (highest resolution)
tracks <- mclapply(seeds, run.track, mc.cores = Ncores)

#assign seed value to each track
tracks_seeds <- Map(function(df, seed) {
  df$seed <- seed
  df
}, tracks, seeds)

tracks_comb <- do.call(rbind, tracks_seeds)

#save the tracks 
save(tracks_comb, file = "simulations/sensitivity/sampling-int-tracks.Rda")

tracks <- split(tracks_comb, tracks_comb$seed)

#thin the tracks 
#set intervals to test (evenly spaced on log2 scale)
#interval is calculated as tau_v / x
x_vals <- 2^seq(log2(100), log2(0.5), length.out = 200)

#create empty lists to store results
thinned <- list()
thin_tracks <- list()

for(j in seq_along(tracks)){
  track_j <- tracks[[j]]
  seed <- unique(track_j$seed)
  
  for(i in seq_along(x_vals)){
    x <- x_vals[i]
    
    t_thin <- sampling(mass_prey, x)
    
    result <- thin(mass_prey, track_j, food, t_thin)
    
    result$summary$x <- x
    result$summary$seed <- seed
    
    thinned[[length(thinned) + 1]] <- result$summary
  }
}

interval_res <- do.call(rbind, thinned)

#save results
save(interval_res, file = 'simulations/sensitivity/sampling-int-res.Rda')

load("simulations/sensitivity/sampling-int-res.Rda")

# data cleaning and plotting
interval_res <- interval_res %>% 
  group_by(seed) %>% 
  mutate(frac = patches / max(patches), 
         interval_95 = max(interval[frac >= 1]), 
         x_95 = x[interval == interval_95],#calculate the proportion of encounters with patches
         seed = as.factor(seed))
  
# plot results
#p <- 
  ggplot(interval_res, aes(x = log2(interval), y = frac)) + #log2 to show the effect of doubling the sampling interval
  geom_point(size = 0.5) +
  scale_x_continuous(labels = function(x)2^x) +
  facet_grid(seed ~ .) +
  geom_smooth(aes(color = seed), method = 'gam', formula = y ~s(x)) +
  labs(y = "Patches visited", x = expression(bold(log2("Sampling Interval")))) +
  geom_vline(xintercept = log2(tau_v/min(interval_res$x_95))) +
  scale_color_paletteer_d("ggsci::nrc_npg") +
  theme.qel() +
  theme(strip.text.y = element_blank()) #removes facet_grid labels

ggsave(p, file = "figures/sensitivity/samplingint_105500g.png",  width = 6, height = 4.5, units = "in", dpi = 600, bg = "white")



food %>% 
  as.data.frame() %>% 
  ggplot(aes(x = x, y = y)) +
  geom_point() + 
  theme_bw()
