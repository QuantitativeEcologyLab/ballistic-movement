
#preamble ----

#set the working directory
setwd("~/ballisticmovement/ballistic-movement")

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

#--------------------------------------------------------------------------
# create function for simulation ------------------------------------------
#--------------------------------------------------------------------------

run_sens <- function(mass_prey,
                     GENS,
                     REPS,
                     n_prey,
                     calories, 
                     width) {
  
  # Create or load FOOD raster based on parameters
  FOOD <- createFoodRaster(mass_prey, 
                           width = width,
                           calories = calories)
  
  # Other initialization 
  prey_res <- list()
  prey_details <- list()
  
  t <- sampling(mass_prey)
  
  PREY_tau_p <- numeric(0)
  PREY_tau_v <- numeric(0)
  PREY_sig <- numeric(0)
  
  for(G in 1:GENS) {
    
    prey <- list()
    
    for(R in 1:REPS){
      
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
          prey_lv <- sqrt((prey_tau_v/prey_tau_p) * prey_sig)
          
          PREY_mods[[i]] <- ctmm(tau = c(prey_tau_p, prey_tau_v),
                                 mu = c(CENTRES[i,1], CENTRES[i,2]),
                                 sigma = prey_sig)
        }
      }
      
      if (G!=1){
        CENTRES <- rbvpois(n = n_prey,
                           a = prey.SIG(mass_prey)*0.75,
                           b = prey.SIG(mass_prey)*0.75,
                           c = 0)
        CENTRES <- scale(CENTRES, scale = FALSE)
        
        PREY_mods <- list()
        for(i in 1:n_prey){
          prey_tau_p <- sample(PREY_tau_p,1) + rnorm(1, 0, 10)
          prey_tau_p <- ctmm:::clamp(prey_tau_p, min = 0.1, max = Inf)
          prey_tau_v <- sample(PREY_tau_v,1) + rnorm(1, 0, 2)
          prey_tau_v <- ctmm:::clamp(prey_tau_v, min = 0.1, max = Inf)
          prey_sig <- sample(PREY_sig,1)
          prey_lv <- sqrt((prey_tau_v/prey_tau_p)*prey_sig)
          
          PREY_mods[[i]] <- ctmm(tau = c(prey_tau_p, prey_tau_v),
                                 mu = c(CENTRES[i,1], CENTRES[i,2]),
                                 sigma = prey_sig)
        }
      }
      #simulate prey movement
      PREY_tracks <- list()
      for(i in 1:n_prey){
        PREY_tracks[[i]] <- simulate(PREY_mods[[i]], t = t)
      }
      
      #extract ids of patches entered
      benefits_prey <- vector("list", n_prey)
      for(i in 1:n_prey){
        benefits_prey[[i]] <- grazing(PREY_tracks[[i]], FOOD, metric = "ids")
      }
      
      #extract number of changes between patches
      patches <- vector("list", n_prey)
      for(i in 1:n_prey){
        patches[[i]] <- attr(benefits_prey[[i]], "patches")
      }
      
      #extract speed from model
      speed <- numeric(n_prey)
      for(i in 1:n_prey){
        speed[[i]] <- speed_val(models = PREY_mods[[i]])
      }
      
      #assign net calories to each individual
      cal_list <- vector("list", n_prey)
      cal_net <- numeric(n_prey)
      costs <- numeric(n_prey)
      for(i in 1:n_prey){
        mass <- if(length(mass_prey) == 1) mass_prey else mass_prey[i]
        
        cal_list[[i]] <- cals_net(IDs = benefits_prey[[i]], 
                                  habitat = FOOD, 
                                  mass = mass, 
                                  models = PREY_mods[[i]],
                                  speed = speed[[i]],
                                  t = t)
        # Defensive: check result is valid
        if (is.list(cal_list[[i]]) &&
            all(c("cal_net") %in% names(cal_list[[i]]))) {
          
          cal_net[i] <- cal_list[[i]]$cal_net
          costs[i] <- cal_list[[i]]$costs
          
        } else {
          cal_net[i] <- NA
          warning(sprintf("Invalid result from cals_net for individual %d", i))
        }
      }
      
      #compute prey offspring
      results <- prey.fitness(mass = mass_prey,
                              cal_net = cal_net)
      
      offspring_prey <- results$offspring
      mass_update_prey <- results$mass_update
      
      #get values
      prey_lvs <- vector()
      prey_TAU_V <- vector()
      prey_TAU_P <- vector()
      prey_SIGMA <- vector()
      for(i in 1:n_prey){
        prey_TAU_V[i] <- PREY_mods[[i]]$tau["velocity"]
        prey_TAU_P[i] <- PREY_mods[[i]]$tau["position"]
        prey_SIGMA[i] <- ctmm:::area.covm(PREY_mods[[i]]$sigma)
        prey_lvs[i] <- sqrt((prey_TAU_V[i]/prey_TAU_P[i])*prey_SIGMA[i])
      }
      
      #summarise
      prey[[R]] <- data.frame(generation = G,
                              tau_p = prey_TAU_P,
                              tau_v = prey_TAU_V,
                              sig = prey_SIGMA,
                              lv = prey_lvs,
                              patches = unlist(patches),
                              cal_net = cal_net,
                              costs = costs,
                              speed = unlist(speed),
                              offspring = unlist(offspring_prey),
                              mass = mass_prey,
                              mass_update = unlist(mass_update_prey),
                              patch_width = width,
                              calories = calories)
    }
    
    prey <- do.call(rbind, prey)
    
    prey_res[[G]] <- data.frame(generation = G, 
                                lv = mean(prey$lv),
                                var = var(prey$lv))
    
    prey_details[[G]] <- prey
    
    #Set up the parameters for the next generation based on
    #Fitness of current generation
    PREY_tau_p <- vector()
    PREY_tau_v <- vector()
    PREY_sig <- vector()
    for(i in 1:nrow(prey)){
      if(prey[i,"offspring"] >0){
        PREY_tau_p <- c(PREY_tau_p,
                        rep(prey[i,"tau_p"], prey[i,"offspring"]))
        
        PREY_tau_v <- c(PREY_tau_v,
                        rep(prey[i,"tau_v"], prey[i,"offspring"]))
        
        PREY_sig <- c(PREY_sig,
                      rep(prey[i,"sig"], prey[i,"offspring"]))
        
      } #Closes the if statement
    }
    # If no offspring, save results and stop simulation
    if(length(PREY_tau_p) == 0 || length(PREY_tau_v) == 0 || length(PREY_sig) == 0){
      warning(sprintf("Simulation stopped early at generation %d due to extinction (no offspring)", G))

      break
    }
    
    #progress report
    print(G)
  }
  
  return(prey_details)
  
}

#--------------------------------------------------------------------------
# example of function use -------------------------------------------------
#--------------------------------------------------------------------------

#number of individuals in arena
n_prey <- 5

#number of arenas
REPS <- 1

#number of generations
GENS <- 5

masses <- seq(1000, 10000, 1000)
calories <- seq(50, 500, 50)

# Define your parameter combinations here
param_grid <- expand.grid(
  mass = masses,
  calories = calories
)

# Prepare a list to store all results
all_results <- list()

for(i in 1:nrow(param_grid)) {
tic(paste("Scenario", i))
  
  params <- param_grid[i,]
  
  width <- round(sqrt(prey.SIG(params$mass))/15)
  
  cat("\n--- Running scenario", i, "of", nrow(param_grid), "---\n")
  cat(sprintf("Mass: %d Calories: %d\n", params$mass, params$calories))
  
  
  # Run the simulation with your exact code inside
  res <- tryCatch({
    run_sens(
      mass_prey = params$mass,
      calories = params$calories,
      GENS = GENS,
      REPS = REPS,
      n_prey = n_prey,
      width = width
    )
  }, error = function(e) {
    warning(sprintf("Scenario %d failed: %s", i, e$message))
    NULL  # store NULL so indexing is preserved
  })
  
  all_results[[i]] <- res
  
  save(all_results, file = "~/ballisticmovement/ballistic-movement/sens_sim_results/mass_calorie_results_width15_v6.Rda")
  
  toc(log = TRUE)
}

all_results_df <- bind_rows(all_results)

save(all_results_df, file = "~/ballisticmovement/ballistic-movement/sens_sim_results/mass_calorie_results_width15_v6_df.Rda")











