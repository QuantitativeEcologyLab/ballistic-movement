# Preamble ----

# Set the working directory
setwd("~/hdrive/GitHub/ballistic-movement")
# Set the random seed
set.seed(123)

# Import necessary packages
library(extraDistr)
library(parallel)
library(ctmm)
library(terra)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(patchwork)
library(tictoc)

# Source the functions (ensure 'functions.R' is available in the working directory)
source("scripts/functions.R")

Ncores <- 20

#----------------------------------------------------------------------
# run the simulation----
#----------------------------------------------------------------------

#calculate what masses are used for the analysis of prey individuals
# masses <- seq(500, 200000, length.out = 20)

# masses to test
#  500  11000  21500  32000  42500  53000  63500  74000  84500  95000 105500 116000 126500 137000 147500 158000 168500 179000 189500 200000

# masses to test first
# 500, 21500, 42500, 63500, 84500, 105500 126500 147500 168500 189500 200000

# Prey mass (g)
mass_prey <- 84500

#set sampling interval and lifespan
t <- sampling(mass_prey, x = 40.5)

#number of individuals in arena
n_prey <- 20

#number of arenas
REPS <- 10

#number of generations
GENS <- 1000

#updated food raster function
FOOD <- createFoodRaster(mass_prey, k = 240000, calories = 8)

#lists for storing results
prey_res <- list()
prey_details <- list()

#start the simulation ----
for(G in 1:GENS) {
  tic(paste("Generation", G))
  
  prey <- list()
  # pred <- list()
  
  for(R in 1:REPS){
    
    # generate movement models
    # if the first generation, generate movement parameters from mass
    if (G == 1) {
      CENTRES <- rbvpois(n = n_prey,
                         a = prey.SIG(mass_prey)*0.75,
                         b = prey.SIG(mass_prey)*0.75,
                         c = 0)
      
      CENTRES <- scale(CENTRES, scale = FALSE)
      
      PREY_mods <- list()
      for(i in 1:n_prey){
        prey_tau_p <- prey.tau_p(mass_prey, variance = TRUE)
        prey_tau_v <- prey.tau_v(mass_prey, variance = TRUE)
        prey_sig <- prey.SIG(mass_prey)
        prey_lv <- sqrt((prey_tau_v/prey_tau_p) * prey_sig)
        
        # create ctmm model
        PREY_mods[[i]] <- ctmm(tau = c(prey_tau_p, prey_tau_v),
                               mu = c(CENTRES[i,1], CENTRES[i,2]),
                               sigma = prey_sig)
      } # closes loop over n_prey
    } # closes generation 1
    
  # if any other generation, draw movement parameters from the offspring pool
    if (G!=1){
      CENTRES <- rbvpois(n = n_prey,
                         a = prey.SIG(mass_prey)*0.75,
                         b = prey.SIG(mass_prey)*0.75,
                         c = 0)
      CENTRES <- scale(CENTRES, scale = FALSE)
      
      PREY_mods <- list()
      for(i in 1:n_prey){
        prey_tau_p <- sample(PREY_tau_p,1)
        prey_tau_v <- sample(PREY_tau_v,1) + rnorm(1, 0, 2) # add 'mutation' based variance
        prey_tau_v <- ctmm:::clamp(prey_tau_v, min = 0.1, max = Inf) # clamp the minimum to 0
        prey_sig <- sample(PREY_sig,1)
        prey_lv <- sqrt((prey_tau_v/prey_tau_p)*prey_sig)
        
        # create ctmm model
        PREY_mods[[i]] <- ctmm(tau = c(prey_tau_p, prey_tau_v),
                               mu = c(CENTRES[i,1], CENTRES[i,2]),
                               sigma = prey_sig)
      } # closes loop over n_prey
      
    } # close if not generation 1
    
    # simulate prey movement
    
    # PREY_tracks <- list()
    # for(i in 1:n_prey){
    #   PREY_tracks[[i]] <- simulate(PREY_mods[[i]], t = t)
    # }
    
    # parallelised to speed up run times
    PREY_tracks <- mclapply(PREY_mods,
                            FUN = simulate,
                            t = t,
                            mc.cores = Ncores)
    
    #extract ids of patches entered
    # benefits_prey <- vector("list", n_prey)
    # for(i in 1:n_prey){
    #   benefits_prey[[i]] <- grazing(PREY_tracks[[i]], FOOD)
    # }
    
    benefits_prey <- mclapply(PREY_tracks, 
                              function(track) grazing(track, FOOD), 
                              mc.cores = Ncores)
    
    #extract number of changes between patches
    patches <- vector("list", n_prey)
    for(i in 1:n_prey){
      patches[[i]] <- attr(benefits_prey[[i]], "patches")
    }
    
    #extract prey speed from model
    prey_speed <- numeric(n_prey)
    for(i in 1:n_prey){
      prey_speed[[i]] <- get.speed(models = PREY_mods[[i]])
    }
    
    #assign net calories to each individual
    prey_cal_list <- vector("list", n_prey)
    prey_cal_net <- numeric(n_prey)
    prey_costs <- numeric(n_prey)
    for(i in 1:n_prey){
      mass <- if(length(mass_prey) == 1) mass_prey else mass_prey[i]

      prey_cal_list[[i]] <- prey.cals.net(IDs = benefits_prey[[i]],
                                          mass = mass,
                                          speed = prey_speed[[i]],
                                          t = t)

        prey_cal_net[i] <- prey_cal_list[[i]]$cal_net
        prey_costs[i] <- prey_cal_list[[i]]$costs
    }
    
    # compute prey offspring
    # output is a list of two variables
    prey_results <- prey.fitness(mass = mass_prey,
                                 cal_net = prey_cal_net)
    
    offspring_prey <- prey_results$offspring #assign offspring 
    mass_update_prey <- prey_results$mass_update #assign mass update
    
    # get prey values
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
    
    #summarise  # save(prey_res, file = 'sim_results/july16/5000000g_longlife_prey_res.Rda')
    # save(prey_details, file = 'sim_results/july16/5000000g_longlife_prey_details.Rda')
    prey[[R]] <- data.frame(generation = G,
                            tau_p = prey_TAU_P,
                            tau_v = prey_TAU_V,
                            sig = prey_SIGMA,
                            lv = prey_lvs,
                            patches = unlist(patches),
                            cal_net = prey_cal_net,
                            costs = prey_costs,
                            speed = unlist(prey_speed),
                            offspring = unlist(offspring_prey),
                            mass = mass_prey,
                            mass_update = unlist(mass_update_prey))
    
  } # closes loop over number of arenas
  
  prey <- do.call(rbind, prey)
  # pred <- do.call(rbind, pred)
  
  # save the results
  # prey
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
  } #closes loop over the number of prey
  
  # If no offspring, save results and stop simulation
    if(length(PREY_tau_p) == 0 || length(PREY_tau_v) == 0 || length(PREY_sig) == 0){
    warning(sprintf("Simulation stopped early at generation %d due to extinction (no offspring)", G))
    
    save(prey_res, file = 'prey_results/84500g_prey_res.Rda')
    save(prey_details, file = 'prey_results/84500g_prey_details.Rda')
    
    break
    }
  
  #save results
  save(prey_res, file = 'prey_results/84500g_prey_res.Rda')
  save(prey_details, file = 'prey_results/84500g_prey_details.Rda')
  
  toc(log = TRUE)
}

