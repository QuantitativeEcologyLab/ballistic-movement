
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
                     fctr) {
  
  # Create or load FOOD raster based on parameters
  FOOD <- createFoodRaster(mass_prey, 
                           patch_width = round(sqrt(pred.SIG(mass_prey))/10),
                           fctr)
  
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
        PREY_tracks[[i]] <- simulate(PREY_mods[[i]], t = t, plot = FALSE)
      }
      
      #extract ids of patches entered
      benefits_prey <- vector("list", n_prey)
      for(i in 1:n_prey){
        benefits_prey[[i]] <- grazing(PREY_tracks[[i]], FOOD, mass = mass_prey)
      }
      
      #extract number of changes between patches
      patches <- vector("list", n_prey)
      for(i in 1:n_prey){
        patches[[i]] <- attr(benefits_prey[[i]], "patches")
      }
      
      #extract unscaled gross kJ
      kJ_gross_unsc <- vector("list", n_prey)
      for(i in 1:n_prey){
        kJ_gross_unsc[[i]] <- attr(benefits_prey[[i]], "kJ_gross_unsc")
      }
      
      #extract speed from model
      speed <- numeric(n_prey)
      for(i in 1:n_prey){
        speed[[i]] <- speed_val(models = PREY_mods[[i]])
      }
      
      #assign net calories to each individual
      kJ_list <- vector("list", n_prey)
      kJ_net <- numeric(n_prey)
      kJ_gross <- numeric(n_prey)
      cost <- numeric(n_prey)
      for(i in 1:n_prey){
        mass <- if(length(mass_prey) == 1) mass_prey else mass_prey[i]
        
        kJ_list[[i]] <- net_kJ_val(kJ_gross = benefits_prey[[i]], 
                                   habitat = FOOD, 
                                   mass = mass, 
                                   t = t,
                                   speed = speed[[i]])
        # Defensive: check result is valid
        if (is.list(kJ_list[[i]]) &&
            all(c("kJ_net") %in% names(kJ_list[[i]]))) {
          
          kJ_net[i] <- kJ_list[[i]]$kJ_net
          kJ_gross[i] <- kJ_list[[i]]$kJ_gross
          cost[i] <- kJ_list[[i]]$cost
          
        } else {
          kJ_net[i] <- NA
          warning(sprintf("Invalid result from cals_net for individual %d", i))
        }
      }
      
      #compute prey offspring
      results <- prey.fitness(mass = mass_prey,
                              kJ_net = kJ_net)
      
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
                              kJ_gross_unsc = unlist(kJ_gross_unsc),
                              kJ_net = kJ_net,
                              kJ_gross = kJ_gross,
                              cost = cost,
                              speed = unlist(speed),
                              offspring = unlist(offspring_prey),
                              mass = mass_prey,
                              mass_update = unlist(mass_update_prey))
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
      if(prey[i,"offspring"] > 0){
        PREY_tau_p <- c(PREY_tau_p,
                        rep(prey[i,"tau_p"], prey[i,"offspring"]))
        
        PREY_tau_v <- c(PREY_tau_v,
                        rep(prey[i,"tau_v"], prey[i,"offspring"]))
        
        PREY_sig <- c(PREY_sig,
                      rep(prey[i,"sig"], prey[i,"offspring"]))
        
      } #Closes the if statement
      
      # If no offspring, save results and stop simulation
      if(length(PREY_tau_p) == 0 || length(PREY_tau_v) == 0 || length(PREY_sig) == 0){
        warning(sprintf("Simulation stopped early at generation %d due to extinction (no offspring)", G))
        break
      }
    }
    
    cat("Completed generation", G, "\n")
  }
  return(prey_details)
}


#--------------------------------------------------------------------------
# example of function use -------------------------------------------------
#--------------------------------------------------------------------------

#number of individuals in arena
n_prey <- 2

#number of arenas
REPS <- 1

#number of generations
GENS <- 5

masses <- seq(1000, 20000, 1000)
fctrs <- seq(10, 100, 5)

# Define your parameter combinations here
param_grid <- expand.grid(
  mass = masses,
  fctr = fctrs
)

# Prepare a list to store all results
all_results <- list()

for(i in 1:nrow(param_grid)) {
tic(paste("Scenario", i))
  
  params <- param_grid[i,]
  
  cat("\n--- Running scenario", i, "of", nrow(param_grid), "---\n")
  cat(sprintf("Mass: %d Factor: %d\n", params$mass, params$fctr))
  
  
  # Run the simulation with your exact code inside
  res <- run_sens(
    mass_prey = params$mass,
    fctr = params$fctr,
    GENS = GENS,
    REPS = REPS,
    n_prey = n_prey
  )
  
  all_results[[i]] <- res
  
  saveRDS(all_results, file = "~/ballisticmovement/ballistic-movement/sens_sim_results/mass_fctr_sens_results.rds")
  
  toc(log = TRUE)
}

all_results_df <- bind_rows(all_results)

save(all_results_df, file = "~/ballisticmovement/ballistic-movement/sens_sim_results/sensitivity_analysis_mass_df.Rda")











