
#preamble ----

#set the working directory
setwd("~/ballisticmovement/ballistic-movement")
#set the random seed
set.seed(1)

#import necessary packages
library(extraDistr)
library(ctmm)
library(terra)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(patchwork)
library(tictoc)

#source the functions (ensure 'functions.R' is available in the working directory)
source("scripts/functions.R")

#----------------------------------#
# create function for simulation   #
#----------------------------------#

run_simulation <- function(mass_prey,
                                  GENS,
                                  REPS,
                                  n_prey,
                                  patch_width,
                                  biomass,
                                  calories,
                                  save_prefix) {
  
  # Create or load FOOD raster based on parameters
  FOOD <- createFoodRaster(mass_prey, 
                           patch_width = patch_width, 
                           calories = calories, 
                           dry_biomass = biomass)
  
  # Other initializations here
  prey_res <- list()
  prey_details <- list()
  
  t <- sampling(mass_prey)
  interval <- sampling(mass_prey, metric = "interval")
  
  PREY_tau_p <- numeric(0)
  PREY_tau_v <- numeric(0)
  PREY_sig <- numeric(0)
  
  for(G in 1:GENS) {
    
    prey <- list()
    
    for(R in 1:REPS) {
      
      if (G == 1) {
        CENTRES <- rbvpois(n = n_prey,
                           a = prey.SIG(mass_prey)*0.75,
                           b = prey.SIG(mass_prey)*0.75,
                           c = 0)
        CENTRES <- scale(CENTRES, scale = FALSE)
        
        PREY_mods <- list()
        for(i in 1:n_prey){
          prey_tau_p <- prey.tau_p(mass_prey, variance = FALSE)
          prey_tau_v <- prey.tau_v(mass_prey, variance = FALSE)
          prey_sig <- prey.SIG(mass_prey)
          
          PREY_mods[[i]] <- ctmm(tau = c(prey_tau_p, prey_tau_v),
                                 mu = c(CENTRES[i,1], CENTRES[i,2]),
                                 sigma = prey_sig)
        }
      }
      
      if (G != 1) {
        CENTRES <- rbvpois(n = n_prey,
                           a = prey.SIG(mass_prey)*0.75,
                           b = prey.SIG(mass_prey)*0.75,
                           c = 0)
        CENTRES <- scale(CENTRES, scale = FALSE)
        
        PREY_mods <- list()
        for(i in 1:n_prey) {
          prey_tau_p <- sample(PREY_tau_p, 1) + rnorm(1, 0, 10)
          prey_tau_p <- ctmm:::clamp(prey_tau_p, min = 0.1, max = Inf)
          prey_tau_v <- sample(PREY_tau_v, 1) + rnorm(1, 0, 2)
          prey_tau_v <- ctmm:::clamp(prey_tau_v, min = 0.1, max = Inf)
          prey_sig <- sample(PREY_sig, 1)
          
          PREY_mods[[i]] <- ctmm(tau = c(prey_tau_p, prey_tau_v),
                                 mu = c(CENTRES[i,1], CENTRES[i,2]),
                                 sigma = prey_sig)
        }
      }
      
      PREY_tracks <- list()
      for(i in 1:n_prey){
        PREY_tracks[[i]] <- simulate(PREY_mods[[i]], t = t)
      }
      
      benefits_prey <- vector("list", n_prey)
      for(i in 1:n_prey){
        benefits_prey[[i]] <- grazing(PREY_tracks[[i]], FOOD, metric = "ids")
      }
      
      patches <- vector("list", n_prey)
      for(i in 1:n_prey){
        patches[[i]] <- grazing(PREY_tracks[[i]], FOOD, metric = "patches")
      }
      
      speed <- numeric(n_prey)
      for(i in 1:n_prey){
        speed[[i]] <- speed_val(models = PREY_mods[[i]])
      }
      
      cal_list <- vector("list", n_prey)
      cal_net <- numeric(n_prey)
      cal_max <- numeric(n_prey)
      for(i in 1:n_prey){
        mass <- if(length(mass_prey) == 1) mass_prey else mass_prey[i]
        
        cal_list[[i]] <- cals_net(IDs = benefits_prey[[i]], 
                                  habitat = FOOD, 
                                  mass = mass, 
                                  models = PREY_mods[[i]],
                                  speed = speed[[i]],
                                  interval = interval)
        if (is.list(cal_list[[i]]) &&
            all(c("cal_net", "cal_max") %in% names(cal_list[[i]]))) {
          cal_net[i] <- cal_list[[i]]$cal_net
          cal_max[i] <- cal_list[[i]]$cal_max
        } else {
          cal_net[i] <- NA
          cal_max[i] <- NA
          warning(sprintf("Invalid result from cals_net for individual %d", i))
        }
      }
      
      results <- prey.fitness(mass = mass_prey,
                              cal_net = cal_net)
      
      offspring_prey <- results$offspring
      mass_update_prey <- results$mass_update
      
      prey_lvs <- numeric(n_prey)
      prey_TAU_V <- numeric(n_prey)
      prey_TAU_P <- numeric(n_prey)
      prey_SIGMA <- numeric(n_prey)
      for(i in 1:n_prey){
        prey_TAU_V[i] <- PREY_mods[[i]]$tau["velocity"]
        prey_TAU_P[i] <- PREY_mods[[i]]$tau["position"]
        prey_SIGMA[i] <- ctmm:::area.covm(PREY_mods[[i]]$sigma)
        prey_lvs[i] <- sqrt((prey_TAU_V[i]/prey_TAU_P[i]) * prey_SIGMA[i])
      }
      
      prey[[R]] <- data.frame(generation = G,
                              tau_p = prey_TAU_P,
                              tau_v = prey_TAU_V,
                              sig = prey_SIGMA,
                              lv = prey_lvs,
                              patches = unlist(patches),
                              cal = cal_net,
                              cal_max = cal_max,
                              speed = unlist(speed),
                              offspring = unlist(offspring_prey),
                              mass = mass_prey,
                              mass_update = unlist(mass_update_prey),
                              patch_width = unlist(patch_width),
                              biomass = unlist(biomass),
                              calories = unlist(calories))
    }
    
    prey <- do.call(rbind, prey)
    
    prey_details[[G]] <- prey
    
    PREY_tau_p <- numeric(0)
    PREY_tau_v <- numeric(0)
    PREY_sig <- numeric(0)
    
    for(i in 1:nrow(prey)){
      if(prey$offspring[i] > 0){
        PREY_tau_p <- c(PREY_tau_p, rep(prey$tau_p[i], prey$offspring[i]))
        PREY_tau_v <- c(PREY_tau_v, rep(prey$tau_v[i], prey$offspring[i]))
        PREY_sig <- c(PREY_sig, rep(prey$sig[i], prey$offspring[i]))
      }
    }
    
    if(length(PREY_tau_p) == 0 || length(PREY_tau_v) == 0 || length(PREY_sig) == 0){
      warning(sprintf("Simulation stopped early at generation %d due to extinction (no offspring)", G))
      break
    }
    
    # Save results each generation to disk
    save(prey_res, file = paste0(save_prefix, "_res.Rda"))
    save(prey_details, file = paste0(save_prefix, "_details.Rda"))
    
    cat("Completed generation", G, "\n")
  }
  return(prey_details)
}


#-------------------------------------------------------------------------------
# create sensitivity analysis
#-------------------------------------------------------------------------------

#------------------------------#
#    Static Simulation Params  #
#------------------------------#

#set sampling interval and lifespan
t <- sampling(mass_prey)
interval <- sampling(mass_prey, metric = "interval")

#number of individuals in arena
n_prey <- 5

#number of arenas
REPS <- 5

#number of generations
GENS <- 10


# Define your parameter combinations here
param_grid <- expand.grid(
  mass_prey = seq(2000, 10000, 2000),
  patch_width = seq(10, 100, 10),
  biomass = seq(1, 15, 3),
  calories = seq(1000, 4000, 500) 
)

# Prepare a list to store all results
all_results <- list()

for(i in seq_len(nrow(param_grid))) {
tic(paste("Scenario", i))
  params <- param_grid[i, ]
  
  cat("\n--- Running scenario", i, "of", nrow(param_grid), "---\n")
  cat(sprintf("Width: %d, Biomass: %.1f, Calories: %d\n", params$patch_width, params$biomass, params$calories))
  
  # Construct a prefix for saving files so each scenario saves to unique files
  save_prefix <- paste0("~/ballisticmovement/ballistic-movement/sim_results/sensitivity",
                        i, "_width_", params$patch_width, "_biomass_", params$biomass, "_calories_", params$calories)
  
  # Run the simulation with your exact code inside
  res <- run_simulation(
    mass_prey = params$mass_prey,
    GENS = GENS,
    REPS = REPS,
    n_prey = n_prey,
    patch_width = params$patch_width,
    biomass = params$biomass,
    calories = params$calories,
    save_prefix = save_prefix
  )
  
  all_results[[i]] <- res
  
  toc(log = TRUE)
}



all_results_df <- bind_rows(all_results)


# Summarize mean and variance of 'lv' over generations per scenario
summary_results <- all_results_df %>%
  group_by(patch_width, biomass, calories) %>%
  summarise(
    mean_off = mean(offspring, na.rm = TRUE),
    var_off = var(offspring, na.rm = TRUE),
    n_gens = n()
  ) %>%
  ungroup()

# Quick plot: effect of patch width on mean lv, faceted by biomass and calories
ggplot(summary_results, aes(x = patch_width, y = mean_off)) +
  geom_point() +
  geom_line() +
  facet_grid(biomass ~ calories) +
  labs(title = "Mean offspring by patch width, biomass, and calories",
       x = "Patch Width",
       y = "Mean offspring") +
  theme_minimal()

# Optional: Simple linear model to test sensitivity
lm_fit <- lm(mean_off ~ patch_width + biomass + calories, data = summary_results)
summary(lm_fit)

# More advanced: Random Forest for non-linear effects (requires 'randomForest' package)
# library(randomForest)
# rf_fit <- randomForest(mean_lv ~ width + biomass + calories, data = summary_results)
# print(rf_fit)



