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

Ncores <- 10

#----------------------------------------------------------------------
# run the simulation----
#----------------------------------------------------------------------

# predator mass 
# mass_pred <- 15000

# Prey mass (g)
mass_prey <- 30000

#set sampling interval and lifespan
t <- sampling(mass_prey, x = 40.5)

#number of individuals in arena
n_prey <- 10
# n_pred <- 1

#number of arenas
REPS <- 5

#number of generations
GENS <- 100

#updated food raster function
FOOD <- createFoodRaster(mass_prey, k = 240000, calories = 1.1)
plot(FOOD)

#lists for storing results
prey_res <- list()
prey_details <- list()

# pred_res <- list()
# pred_details <- list()

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
      
      # PRED_mods <- list()
      # for(i in 1:n_pred){
      #   pred_tau_p <- pred.tau_p(mass_pred, variance = TRUE)
      #   pred_tau_v <- pred.tau_v(mass_pred, variance = TRUE)
      #   pred_sig <- pred.SIG(mass_pred)
      #   pred_lv <- sqrt((pred_tau_v / pred_tau_p) * pred_sig)
      # 
      #   PRED_mods[[i]] <- ctmm(tau = c(pred_tau_p,
      #                                  pred_tau_v),
      #                          mu = c(0,0),
      #                          sigma = pred_sig)
      # } # closes loop over n_pred
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
        prey_tau_p <- sample(PREY_tau_p,1) + rnorm(1, 0, 10)
        prey_tau_p <- ctmm:::clamp(prey_tau_p, min = 0.1, max = Inf) 
        prey_tau_v <- sample(PREY_tau_v,1) + rnorm(1, 0, 2) # add 'mutation' based variance
        prey_tau_v <- ctmm:::clamp(prey_tau_v, min = 0.1, max = Inf) # clamp the minimum to 0
        prey_sig <- sample(PREY_sig,1)
        prey_lv <- sqrt((prey_tau_v/prey_tau_p)*prey_sig)
        
        # create ctmm model
        PREY_mods[[i]] <- ctmm(tau = c(prey_tau_p, prey_tau_v),
                               mu = c(CENTRES[i,1], CENTRES[i,2]),
                               sigma = prey_sig)
      } # closes loop over n_prey
      
      # PRED_mods <- list()
      # for(i in 1:n_pred){
      #   pred_tau_p <- sample(PRED_tau_p,1) + rnorm(1, 0, 10) # add some 'mutation' based variance
      #   pred_tau_p <- ctmm:::clamp(pred_tau_p, min = 0.1, max = Inf) # clamp the minimum to 0
      #   pred_tau_v <- sample(PRED_tau_v,1) + rnorm(1, 0, 2) # add some 'mutation' based variance
      #   pred_tau_v <- ctmm:::clamp(pred_tau_v, min = 0.1, max = Inf) # clamp the minimum to 0
      #   pred_sig <- sample(PRED_sig,1)
      #   pred_lv <- sqrt((pred_tau_v/pred_tau_p)*pred_sig)
      # 
      #   PRED_mods[[i]] <- ctmm(tau = c(pred_tau_p,
      #                                  pred_tau_v),
      #                          mu = c(0,0),
      #                          sigma = pred_sig)
      # } # closes loop over n_pred
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
    
    # simualte predator movement
    
    # PRED_tracks <- list()
    # for(i in 1:n_pred){
    #   PRED_tracks[[i]] <- simulate(PRED_mods[[i]], t = t)
    # }
    
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
    
    #extract predator speed from model
    # pred_speed <- numeric(n_pred)
    # for(i in 1:n_pred){
    #   pred_speed[[i]] <- speed_val(models = PRED_mods[[i]])
    # }
    
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
    
    # prey_cal_net <- numeric(n_prey)
    # for(i in 1:n_prey){
    #   prey_cal_net[[i]] <- prey_cals_net_nocost(IDs = benefits_prey[[i]])
    # }
    
    # count the encounters (only setup for a single predator/arena)
    # encounters <- encounter(prey.tracks = PREY_tracks,
    #                         pred.tracks = PRED_tracks,
    #                         range = sqrt(pred.SIG(mass_pred))*0.05) # perceptual range scaled to HR size
    
    # assign net calories to each predator 
    # pred_cal_list <- vector("list", n_pred)
    # pred_cal_net <- numeric(n_pred)
    # pred_costs <- numeric(n_pred)
    # for(i in 1:n_pred){
    #   mass <- if(length(mass_pred) == 1) mass_pred else mass_pred[i]
    # 
    #   pred_cal_list[[i]] <- pred_cals_net(encounters = encounters,
    #                                  mass = mass_pred,
    #                                  t = t,
    #                                  speed = pred_speed)
    # 
    #   pred_cal_net[i] <- pred_cal_list[[i]]$cal_net
    #   pred_costs[i] <- pred_cal_list[[i]]$costs
    # }
    
    # compute prey offspring
    # output is a list of two variables
    prey_results <- prey.fitness(mass = mass_prey,
                                 cal_net = prey_cal_net)
    
    offspring_prey <- prey_results$offspring #assign offspring 
    mass_update_prey <- prey_results$mass_update #assign mass update
    
    # compute predator offspring
    # pred_results <- pred.fitness(mass = mass_pred,
    #                              pred_cal_net = pred_cal_net)
    # 
    # offspring_pred <- pred_results$offspring
    # mass_update_pred <- pred_results$mass_update
    
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
    
    # get predator values
     # pred_lvs <- vector()
     # pred_TAU_V <- vector()
     # pred_TAU_P <- vector()
     # pred_SIGMA <- vector()
     # for(i in 1:n_pred){
     #   pred_TAU_V[i] <- PRED_mods[[i]]$tau["velocity"]
     #   pred_TAU_P[i] <- PRED_mods[[i]]$tau["position"]
     #   pred_SIGMA[i] <- ctmm:::area.covm(PRED_mods[[i]]$sigma)
     #   pred_lvs[i] <- sqrt((pred_TAU_V[i]/pred_TAU_P[i])*pred_SIGMA[i])
     # }
     # 
     # pred[[R]] <- data.frame(generation = G,
     #                         tau_p = pred_TAU_P,
     #                         tau_v = pred_TAU_V,
     #                         sig = pred_SIGMA,
     #                         lv = pred_lvs,
     #                         pred_cal_net = pred_cal_net,
     #                         pred_costs = pred_costs,
     #                         speed = unlist(pred_speed),
     #                         offspring = unlist(offspring_pred),
     #                         mass = mass_pred,
     #                         mass_update = unlist(mass_update_pred))
    
  } # closes loop over number of arenas
  
  prey <- do.call(rbind, prey)
  # pred <- do.call(rbind, pred)
  
  # save the results
  # prey
  prey_res[[G]] <- data.frame(generation = G, 
                              lv = mean(prey$lv),
                              var = var(prey$lv))
  
  prey_details[[G]] <- prey
  
  # predator
  # pred_res[[G]] <- data.frame(generation = G,
  #                             lv = mean(pred$lv),
  #                             var = var(pred$lv))
  # 
  # pred_details[[G]] <- pred
  
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
  
  #Fitness of current generation
  # PRED_tau_p <- vector()
  # PRED_tau_v <- vector()
  # PRED_sig <- vector()
  # for(i in 1:nrow(pred)){
  #   if(pred[i,"offspring"] >0){
  #     PRED_tau_p <- c(PRED_tau_p,
  #                     rep(pred[i,"tau_p"], pred[i,"offspring"]))
  # 
  #     PRED_tau_v <- c(PRED_tau_v,
  #                     rep(pred[i,"tau_v"], pred[i,"offspring"]))
  # 
  #     PRED_sig <- c(PRED_sig,
  #                   rep(pred[i,"sig"], pred[i,"offspring"]))
  # 
  #   } # closes the if statement
  # } #closes loop over the number of pred
  
  # If no offspring, save results and stop simulation
    if(length(PREY_tau_p) == 0 || length(PREY_tau_v) == 0 || length(PREY_sig) == 0){
    warning(sprintf("Simulation stopped early at generation %d due to extinction (no offspring)", G))
    
    # save(prey_res, file = 'sim_results/july16/5000000g_longlife_prey_res.Rda')
    # save(prey_details, file = 'sim_results/july16/5000000g_longlife_prey_details.Rda')
    
    break
    }
  
  #save results
  # save(prey_res, file = 'sim_results/july16/5000000g_longlife_prey_res.Rda')
  # save(prey_details, file = 'sim_results/july16/5000000g_longlife_prey_details.Rda')
  
  # save predator results
  # save(pred_res, file = 'sim_results/constant_resources_withpred/June23_15000gpredator_pred_res.Rda')
  # save(pred_details, file = 'sim_results/constant_resources_withpred/June23_15000gpredator_pred_details.Rda')
  
  toc(log = TRUE)
}

