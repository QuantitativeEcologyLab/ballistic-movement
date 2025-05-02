#----------------------------------------------------------------------
# Preamble
#----------------------------------------------------------------------

# Set the working directory
setwd("C:/Users/lterp/Desktop/OneDrive - UBC/_URAproject/ballistic_movement_models/ballistic_movement_models")

# Set the random seed
set.seed(1)

# Import necessary packages
library(extraDistr)
library(ctmm)
library(terra)

# Source the functions (ensure 'functions.R' is available in the working directory)
source("functions.R")

#----------------------------------------------------------------------
# Set up the global parameters for the simulation
#----------------------------------------------------------------------

# Predator mass (g)
mass_pred <- 5000

# Prey mass (g)
mass_prey <- 1000

# "Lifespan" and sampling interval for the simulations
t <- sampling2.0(mass_prey, crossings = 20, risk_factor = 0.05)

# Energetic value of a patch (based on prey BMR)
CALS <- ((10^(0.774 + 0.727*log10(mass_prey)))^1.22)/150

# Number of preds & prey in each "arena"
n_prey <- floor(1/mass_pred^(-0.25))
n_pred <- 1

# Number of "arenas"
REPS <- 2

# Number of generations
GENS <- 2

# Build the raster of food patches for prey to feed on
FOOD <- patches(mass_pred,
                width = 20,
                pred = T, type = "uniform")

#----------------------------------------------------------------------
# Run the simulation
#----------------------------------------------------------------------

# Initialize lists for storing results
prey_res <- list()
prey_details <- list()
pred_res <- list()
pred_details <- list()

# Loop over the number of generations
for(G in 1:GENS){
  
  # Empty lists for storing the results of the current generation
  prey <- list()
  pred <- list()
  
  cat("Starting Generation", G, "\n")  # Debugging output
  
  for(R in 1:REPS){
    
    cat("Running Arena", R, "\n")  # Debugging output
    
    # Generate the prey movement models
    # If the first gen, generate movement parameters from the mass functions
    if(G == 1){
      # Generate the HR centres of the prey
      CENTRES <- rbvpois(n = n_prey,
                         a = pred.SIG(mass_pred)*.75,
                         b = pred.SIG(mass_pred)*.75,
                         c = 0)
      CENTRES <- scale(CENTRES, scale = FALSE)
      
      PREY_mods <- list()
      for(i in 1:n_prey){
        # Prey movement parameters
        prey_tau_p <- prey.tau_p(mass_prey, variance = TRUE)
        prey_tau_v <- sample(1:(prey_tau_p-1), 1)
        prey_sig <- prey.SIG(mass_prey)
        prey_lv <- sqrt((prey_tau_v/prey_tau_p)*prey_sig)
        
        PREY_mods[[i]] <- ctmm(tau = c(prey_tau_p,prey_tau_v),
                               mu = c(CENTRES[i,1],CENTRES[i,2]),
                               sigma = prey_sig)
      }
      
      PRED_mods <- list()
      for(i in 1:n_pred){
        # Predator movement parameters
        pred_tau_p <- pred.tau_p(mass_pred, variance = TRUE)
        pred_tau_v <- sample(1:(pred_tau_p-1), 1)
        pred_sig <- pred.SIG(mass_pred)
        pred_lv <- sqrt((pred_tau_v/pred_tau_p)*pred_sig)
        
        PRED_mods[[i]] <- ctmm(tau = c(pred_tau_p,
                                       pred_tau_v),
                               mu = c(0,0),
                               sigma = pred_sig)
      }
    }  # End of first generation
    
    # Simulate the prey movement using lapply (no parallelization)
    prey.tracks <- lapply(PREY_mods,
                          FUN = simulate,
                          t = t)
    
    # Debugging: Check if prey.tracks are populated
    cat("Simulated Prey Tracks\n")
    if (length(prey.tracks) == 0) {
      cat("No prey tracks generated! \n")
    }
    
    # Simulate the predator movement
    pred.tracks <- list()
    for(i in 1:n_pred){
      pred.tracks[[i]] <- simulate(PRED_mods[[i]], t = t)
    }
    
    # Calculate prey benefits
    benefits_prey <- vector()
    for(i in 1:n_prey){
      benefits_prey[i] <- grazing(prey.tracks[[i]], FOOD)
    }
    
    # Debugging: Check if benefits_prey are populated
    cat("Calculated Prey Benefits\n")
    if (length(benefits_prey) == 0) {
      cat("No benefits for prey! \n")
    }
    
    # Count the encounters (only setup for single predator/arena)
    encounters <- encounter(prey.tracks = prey.tracks,
                            pred.tracks = pred.tracks,
                            range = sqrt(pred.SIG(mass_pred))*.05) # Perceptual range scaled to HR size
    
    # Calculate prey fitness
    offspring_prey <- prey.fitness(benefits = benefits_prey,
                                   costs = encounters,
                                   mass_prey,
                                   crossings = 30,
                                   models = PREY_mods,
                                   calories = CALS)
    
    # Calculate predator fitness
    offspring_pred <- pred.fitness(encounters = encounters,
                                   mass = mass_pred,
                                   models = PRED_mods,
                                   time = t)
    
    # Get the values of the prey movement model parameters
    prey_lvs <- vector()
    prey_TAU_V <- vector() 
    prey_TAU_P <- vector()
    prey_SIGMA <- vector()
    prey_SPEED <- vector()
    for(i in 1:n_prey){
      prey_TAU_V[i] <- PREY_mods[[i]]$tau["velocity"]
      prey_TAU_P[i] <- PREY_mods[[i]]$tau["position"]
      prey_SIGMA[i] <- ctmm:::area.covm(PREY_mods[[i]]$sigma)
      prey_SPEED[i] <- if(nrow(summary(PREY_mods[[i]], units = FALSE)$CI)==4){
        summary(PREY_mods[[i]], units = FALSE)$CI[4,2]} else{Inf}
      prey_lvs[i] <- sqrt((prey_TAU_V[i]/prey_TAU_P[i])*prey_SIGMA[i])
    }
    
    # Summarize the results of the prey 
    prey[[R]] <- data.frame(generation = G,
                            tau_p = prey_TAU_P,
                            tau_v = prey_TAU_V,
                            sig = prey_SIGMA,
                            speed = prey_SPEED,
                            lv = prey_lvs,
                            patches = benefits_prey,
                            encounter = encounters,
                            offspring = offspring_prey)
    
    # Debugging: Check if prey data is generated
    cat("Prey Data Summary:\n")
    print(head(prey[[R]]))
    
    # Get the values of the predator movement model parameters
    pred_lvs <- vector()
    pred_TAU_V <- vector() 
    pred_TAU_P <- vector()
    pred_SIGMA <- vector()
    pred_SPEED <- vector()
    for(i in 1:n_pred){
      pred_TAU_V[i] <- PRED_mods[[i]]$tau["velocity"]
      pred_TAU_P[i] <- PRED_mods[[i]]$tau["position"]
      pred_SIGMA[i] <- ctmm:::area.covm(PRED_mods[[i]]$sigma)
      pred_SPEED[i] <- if(nrow(summary(PRED_mods[[i]], units = FALSE)$CI)==4){
        summary(PRED_mods[[i]], units = FALSE)$CI[4,2]} else{Inf}
      pred_lvs[i] <- sqrt((pred_TAU_V[i]/pred_TAU_P[i])*pred_SIGMA[i])
    }
    
    # Summarize the results of the predator
    pred[[R]] <- data.frame(generation = G,
                            tau_p = pred_TAU_P,
                            tau_v = pred_TAU_V,
                            sig = pred_SIGMA,
                            speed = pred_SPEED,
                            lv = pred_lvs,
                            encounter = sum(encounters),
                            offspring = offspring_pred)
    
  }
  
  # Compile the results from the generation
  prey <- do.call(rbind, prey)
  pred <- do.call(rbind, pred)
  
  # Save the results for the generation
  prey_res[[G]] <- data.frame(generation = G,
                              lv = mean(prey$lv),
                              var = var(prey$lv),
                              pred_lv = mean(pred$lv))
  
  prey_details[[G]] <- prey
  
  pred_res[[G]] <- data.frame(generation = G,
                              lv = mean(pred$lv),
                              var = var(pred$lv))
  
  pred_details[[G]] <- pred
  
  # Debugging: Print the results
  cat("Generation", G, "Results:\n")
  print(prey_res[[G]])
  print(pred_res[[G]])
  
}

