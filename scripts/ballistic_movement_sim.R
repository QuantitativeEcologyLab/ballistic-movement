# Preamble ----

# Set the working directory
setwd("~/H/GitHub/ballistic-movement")
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
mass_prey <- 5000000

#set sampling interval and lifespan
t <- sampling(mass_prey)

#number of individuals in arena
n_prey <- 10
# n_pred <- 1

#number of arenas
REPS <- 1

#number of generations
GENS <- 100

#build food raster
FOOD <- createFoodRaster(mass_prey, calories = 0.015, width = round(sqrt(prey.SIG(mass_prey)))/10, heterogeneity = FALSE)

#lists for storing results
prey_res <- list()
prey_details <- list()

# pred_res <- list()
# pred_details <- list()

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
        prey_tau_v <- sample(PREY_tau_v,1) + rnorm(1, 0, 2)
        prey_tau_v <- ctmm:::clamp(prey_tau_v, min = 0.1, max = Inf)
        prey_sig <- sample(PREY_sig,1)
        prey_lv <- sqrt((prey_tau_v/prey_tau_p)*prey_sig)
        
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
      prey_speed[[i]] <- speed_val(models = PREY_mods[[i]])
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

      prey_cal_list[[i]] <- prey_cals_net(IDs = benefits_prey[[i]],
                                          mass = mass,
                                          speed = prey_speed[[i]],
                                          t = t)

        prey_cal_net[i] <- prey_cal_list[[i]]$cal_net
        prey_costs[i] <- prey_cal_list[[i]]$costs
    }
    
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
    prey_results <- prey.fitness(mass = mass_prey,
                                 cal_net = prey_cal_net)
    
    offspring_prey <- prey_results$offspring
    mass_update_prey <- prey_results$mass_update
    
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
    
    #summarise
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
    
    #save(prey_res, file = 'sim_results/calories_sensitivity/5000000g_11cal_prey_res.Rda')
    #save(prey_details, file = 'sim_results/calories_sensitivity/5000000g_11cal_prey_details.Rda')  
    
    break
    }
  
  #save results
  #save(prey_res, file = 'sim_results/calories_sensitivity/5000000g_11cal_prey_res.Rda')
  #save(prey_details, file = 'sim_results/calories_sensitivity/5000000g_11cal_prey_details.Rda')    
  
  # save predator results
  # save(pred_res, file = 'sim_results/constant_resources_withpred/June23_15000gpredator_pred_res.Rda')
  # save(pred_details, file = 'sim_results/constant_resources_withpred/June23_15000gpredator_pred_details.Rda')
  
  toc(log = TRUE)
}


#----------------------------------------------------------------------
# make diagnostic figures----
#----------------------------------------------------------------------

# make data sets compatible
prey_res_df <- do.call(rbind, prey_res)
prey_details_df <- do.call(rbind, prey_details)

# pred_res_df <- do.call(rbind, pred_res)
# pred_details_df <- do.call(rbind, pred_details)


# relative change in lv ~ gen
PREY_LV <- prey_res_df$lv[1]
prey_res_df$rel.lv <- prey_res_df$lv/PREY_LV
prey_res_df$rel_var <- prey_res_df$var / (PREY_LV^2)
prey_res_df$rel_sd <- sqrt(prey_res_df$rel_var)

prey_summary <- prey_details_df %>%
  group_by(generation) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop")

prey_details_df_sum <- prey_details_df %>% 
  filter(generation >= max(generation) - 49) %>% 
  ungroup()

# rel lv ~ generation
rel.lv.gen <- ggplot(prey_res_df, aes(x = generation, y = rel.lv)) +
  geom_line(color = "deeppink4", linewidth = 1) +
  geom_hline(yintercept = 1, color = 'grey30', linetype = "dashed") +
  geom_ribbon(aes(ymin = rel.lv - rel_sd,
                  ymax = rel.lv + rel_sd),
              fill = "deeppink4", alpha = 0.3) +
  labs(y = "relative change in lv",x = "generation") +
  theme_minimal()
# print(rel.lv.gen)

# tau_v ~ generation
tauv.gen <- ggplot(prey_summary, aes(x = generation, y = tau_v)) +
  stat_summary(fun = mean, geom = "line", col = "deeppink4", linewidth = 1)+
  labs(x = "generation", y = "tau_v") +
  theme_minimal()
# print(tauv.gen)

# tau_p ~ generation
taup.gen <- ggplot(prey_summary, aes(x = generation, y = tau_p)) +
  stat_summary(fun = mean, geom = "line", col = "deeppink4", linewidth = 1)+
  labs(x = "generation", y = "tau_p") +
  theme_minimal()
# print(taup.gen)

# sig ~ generation
sig.gen <- ggplot(prey_summary, aes(x = generation, y = sig)) +
  stat_summary(fun = mean, geom = "line", col = "deeppink4", linewidth = 1) +
  labs(x = "generation", y = "sig") +
  theme_minimal()
# print(sig.gen)

# lv ~ generation
lv.gen <- ggplot(prey_summary, aes(x = generation, y = lv)) +
  stat_summary(fun = mean, geom = "line", col = "deeppink4", linewidth = 1) +
  labs(x = "generation", y = "lv") +
  theme_minimal()
# print(lv.gen)

# patches ~ generation
patches.gen <- ggplot(prey_summary, aes(x = generation, y = patches)) +
  stat_summary(fun = mean, geom = "line", col = "deeppink4", linewidth = 1) +
  labs(x = "generation", y = "patches visited") +
  theme_minimal()
# print(patches.gen)

# cost ~ generation
cost.gen <- ggplot(prey_summary, aes(x = generation, y = costs)) +
  stat_summary(fun = mean, geom = "line", col = "deeppink4", linewidth = 1) +
  labs(x = "generation", y = "metabolic costs (cal)") +
  theme_minimal()
# print(cost.gen)

# cal_net ~ generation
cal.gen <- ggplot(prey_summary, aes(x = generation, y = cal_net)) +
  stat_summary(fun = mean, geom = "line", col = "deeppink4", linewidth = 1) +
  labs(x = "generation", y = "net calories") +
  theme_minimal()
# print(cal.gen)

# speed ~ generation
speed.gen <- ggplot(prey_summary, aes(x = generation, y = speed)) +
  stat_summary(fun = mean, geom = "line", col = "deeppink4", linewidth = 1) +
  labs(x = "generation", y = "speed (m/s)") +
  theme_minimal()
# print(speed.gen)

# offspring ~ generation
offspring.gen <- ggplot(prey_summary, aes(x = generation, y = offspring)) +
  stat_summary(fun = mean, geom = "line", col = "deeppink4", linewidth = 1) +
  labs(x = "generation", y = "offspring") +
  theme_minimal()
# print(offspring.gen)

# mass ~ gen
mass.gen <- ggplot(prey_summary, aes(x = generation, y = mass_update)) +
  stat_summary(fun = mean, geom = "line", col = "deeppink4", linewidth = 1) +
  labs(x = "generation", y = "mass_updated (g)") +
  theme_minimal()
# print(mass.gen)

# BMR ~ generation
# BMR.gen <- ggplot(prey_summary, aes(x = generation, y = BMR)) +
#   stat_summary(fun = mean, geom = "line", col = "deeppink4", linewidth = 1) +
#   labs(x = "generation",y = "BMR (cal)") +
#   theme_minimal()
# print(BMR.gen)

# Movement cost ~ generation
# move.gen <- ggplot(prey_summary, aes(x = generation, y = Move)) +
#   stat_summary(fun = mean, geom = "line", col = "deeppink4", linewidth = 1) +
#   labs(x = "generation", y = "movement costs (cal)") +
#   theme_minimal()
# print(move.gen)

# number of patches ~ lv
patches.lv <- ggplot(prey_details_df, aes(x = lv, y = patches)) +
  geom_point(col = "deeppink4", alpha = 0.3) +
  labs(x = "lv", y = "patches visited") +
  theme_minimal()
# print(patches.lv)

# number of patches ~ speed
patches.speed <- ggplot(prey_details_df, aes(x = speed, y = patches)) +
  geom_point(col = "deeppink4", alpha = 0.3) +
  labs(x = "speed (m/s)", y = "patches visited") +
  theme_minimal()

# patches ~ offspring
patches.off <- ggplot(prey_details_df, aes(x = patches, y = offspring)) +
  geom_point(col = "deeppink4", alpha = 0.3) +
  labs(x = "patch visited", y = "offspring") + 
  theme_minimal()

# cal_net ~ mass_update
cal.mass <- ggplot(prey_details_df, aes(x = mass_update, y = cal_net)) +
  geom_point(col = "deeppink4", alpha = 0.3) +
  labs(x = "mass_update (g)", y = "net calories") +
  theme_minimal()

# cal_net ~ lv
cal.lv <- ggplot(prey_details_df, aes(x = lv, y = cal_net)) +
  geom_point(col = "deeppink4", alpha = 0.3) +
  labs(x = "lv", y = "net calories") +
  theme_minimal()

# cost ~ mass_update
cost.mass <- ggplot(prey_details_df, aes(x = costs, y = mass_update)) +
  geom_point(col = "deeppink4", alpha = 0.3) +
  labs(x = "movement costs (cal)", y = "mass_update (g)") +
  theme_minimal()
# print(cost.mass)

# cost ~ speed
cost.speed <- ggplot(prey_details_df, aes(x = speed, y = costs)) +
  geom_point(col = "deeppink4", alpha = 0.3) +
  labs(x = "speed (m/s)", y = "metabolic costs (cal)") + 
  theme_minimal()

# speed ~ cal_net 
speed.cal <- ggplot(prey_details_df, aes(x = speed, y = cal_net)) +
  geom_point(col = "deeppink4", alpha = 0.3) +
  labs(x = "speed", y = "cal_net") +
  theme_minimal()
# print(speed.cal)

# speed ~ lv
speed.lv <- ggplot(prey_details_df, aes(x = lv, y = speed)) +
  geom_point(col = "deeppink4", alpha = 0.3) +
  labs(x = "lv", y = "speed (m/s)") +
  theme_minimal()

# speed ~ mass_update
speed.mass <- ggplot(prey_details_df, aes(x = speed, y = mass_update)) +
  geom_point(col = "deeppink4", alpha = 0.3) +
  labs(x = "speed (m/s)", y = "mass_update (g)") +
  theme_minimal()
# print(speed.mass)

# offspring ~ speed
offspring.speed <- ggplot(prey_details_df, aes(x = speed, y = offspring)) +
  geom_point(col = "deeppink4", alpha = 0.3) +
  labs(x = "speed", y = "offspring") +
  theme_minimal()

# offspring ~ lv
offspring.lv <- ggplot(prey_details_df, aes(x = lv, y = offspring)) +
  geom_point(col = "deeppink4", alpha = 0.3) +
  labs(x = "lv", y = "offspring") +
  theme_minimal()

# offspring ~ cal_net
offspring.cal <- ggplot(prey_details_df, aes(x = cal_net, y = offspring)) +
  geom_point(col = "deeppink4", alpha = 0.3) +
  labs(x = "net calories", y = "offspring") +
  theme_minimal()

# offspring ~ mass_update
offspring.mass <- ggplot(prey_details_df, aes(x = mass_update, y = offspring)) +
  geom_point(col = "deeppink4", alpha = 0.3) +
  labs(x = "mass_update (g)", y = "offspring") +
  theme_minimal()

# Move/BMR ~ generation
# BMR.move <- ggplot(prey_details_df, aes(y = Move/BMR, x = generation)) +
#   geom_point(col = "deeppink4", alpha = 0.3) +
#   labs(y = "movement costs/BMR", x = "generation") +
#   theme_minimal()

# speed ~ tau_p
taup.speed <- ggplot(prey_details_df, aes(y = speed, x = tau_p)) +
  geom_point(col = "deeppink4", alpha = 0.3) +
  labs(y = "speed (m/s)", x = "tau_p (s)") +
  theme_minimal()

plots <- list(rel.lv.gen,
              taup.gen,
              tauv.gen,
              sig.gen,
              lv.gen,
              patches.gen,
              # cost.gen,
              offspring.gen,
              cal.gen,
              mass.gen,
              speed.gen,
              patches.lv,
              speed.lv,
              patches.speed,
              cal.lv,
              offspring.lv,
              # cost.speed,
              offspring.speed,
              cal.mass,
              speed.mass,
              speed.cal,
              # cost.mass,
              offspring.mass,
              offspring.cal,
              patches.off,
              taup.speed)

final.plot <- wrap_plots(plots[1:10], ncol = 4)
final <- final.plot + plot_annotation('Panel 1: 5000000g, 11 calories')
print(final)
ggsave("~/H/GitHub/ballistic-movement/sim_results/calories_sensitivity/5000000g_11cal_v2_panel1.PNG", plot = final, width = 15, height = 8, dpi = 800)

final.plot2 <- wrap_plots(plots[11:23], ncol = 4)
final2 <- final.plot2 + plot_annotation('Panel 2: 5000000g, 11 calories')
print(final2)
ggsave("~/H/GitHub/ballistic-movement/sim_results/calories_sensitivity/5000000g_11cal_v2_panel2.PNG", plot = final2, width = 15, height = 8, dpi = 800)
 

gen.plots <- grid.arrange(rel.lv.gen, taup.gen, tauv.gen, speed.gen, ncol = 2, nrow = 2)

ggsave(gen.plots, file = 'sim_results/sensitivity/30000g_finalpanel1.PNG', 
       width = 11, height = 5, units = "in")

other.plots <- grid.arrange(patches.lv, speed.lv, patches.speed, cal.lv, speed.cal, offspring.lv, ncol = 2, nrow = 3)

ggsave(other.plots, file = 'sim_results/sensitivity/30000g_finalpanel2.PNG',
       width = 11, height = 8, units = "in")




