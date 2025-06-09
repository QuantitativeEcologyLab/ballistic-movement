

# Preamble ----

# Set the working directory
setwd("~/H/GitHub/ballistic-movement")
# Set the random seed
set.seed(123)

# Import necessary packages
library(extraDistr)
library(ctmm)
library(terra)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(patchwork)
library(tictoc)

# Source the functions (ensure 'functions.R' is available in the working directory)
source("scripts/functions.R")

#----------------------------------------------------------------------
# run the simulation----
#----------------------------------------------------------------------

# Prey mass (g)
mass_prey <- 1000

#set sampling interval and lifespan
t <- sampling(mass_prey)

#number of individuals in arena
n_prey <- 10

#number of arenas
REPS <- 5

#number of generations
GENS <- 500

#build food raster
FOOD <- createFoodRaster(mass_prey, calories = 10, width = round(sqrt(prey.SIG(mass_prey)))/10)
#lists for storing results
prey_res <- list()
prey_details <- list()

for(G in 1:GENS) {
  tic(paste("Generation", G))
  
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
      benefits_prey[[i]] <- grazing(PREY_tracks[[i]], FOOD)
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
    
      save(prey_res, file = "~/H/GitHub/ballistic-movement/sim_results/June6_lv_Evo_1000g_10cal_prey_res.Rda")
      save(prey_details, file = "~/H/GitHub/ballistic-movement/sim_results/June6_lv_Evo_1000g_10cal_prey_details.Rda")
      
    break
    }
  
  save(prey_res, file = "~/H/GitHub/ballistic-movement/sim_results/June6_lv_Evo_1000g_10cal_prey_res.Rda")
  save(prey_details, file = "~/H/GitHub/ballistic-movement/sim_results/June6_lv_Evo_1000g_10cal_prey_details.Rda")
  
  #progress report
  print(G)
  
  toc(log = TRUE)
}


print(prey_res)
print(prey_details)

#----------------------------------------------------------------------
# make diagnostic figures----
#----------------------------------------------------------------------

load("H:/GitHub/ballistic-movement/sim_results/constant_resources_nopred_varycalories/June6_lv_Evo_1000g_10cal_prey_details.Rda")
load("H:/GitHub/ballistic-movement/sim_results/constant_resources_nopred_varycalories/June6_lv_Evo_1000g_10cal_prey_res.Rda")

#make data sets compatible
prey_res_df <- do.call(rbind, prey_res)
print(prey_res_df)
prey_details_df <- do.call(rbind, prey_details)


# relative change in lv ~ generation
PREY_LV <- prey_res_df$lv[1]
prey_res_df$rel.lv <- prey_res_df$lv/PREY_LV
prey_res_df$rel_var <- prey_res_df$var / (PREY_LV^2)
prey_res_df$rel_sd <- sqrt(prey_res_df$rel_var)

# rel lv ~
rel.lv.gen <- ggplot(prey_res_df, aes(x = generation, y = rel.lv)) +
  geom_line(color = "deeppink4", linewidth = 1) +
  geom_hline(yintercept = 1, color = 'grey30', linetype = "dashed") +
  geom_ribbon(aes(ymin = rel.lv - rel_sd,
                  ymax = rel.lv + rel_sd),
              fill = "deeppink4", alpha = 0.3) +
  labs(y = "relative change in lv",x = "generation") +
  theme_minimal()

# tau_v ~ generation
tauv.gen <- ggplot(prey_details_df, aes(x = generation, y = tau_v)) +
  stat_summary(fun = mean, geom = "line", col = "deeppink4", linewidth = 1)+
  labs(y = "tau_v", x = "generation") +
  theme_minimal()

# tau_p ~ generation
taup.gen <- ggplot(prey_details_df, aes(x = generation, y = tau_p)) +
  stat_summary(fun = mean, geom = "line", col = "deeppink4", linewidth = 1)+
  labs(y = "tau_p", x = "generation") +
  theme_minimal()

# sig ~ generation
sig.gen <- ggplot(prey_details_df, aes(x = generation, y = sig)) +
  stat_summary(fun = mean, geom = "line", col = "deeppink4", linewidth = 1) +
  labs(y = "sig", x = "generation") +
  theme_minimal()

# lv ~ mass_update
lv.mass <- ggplot(prey_details_df, aes(x = generation, y = lv)) +
  stat_summary(fun = mean, geom = "line", col = "deeppink4", linewidth = 1) +
  labs(y = "lv", x = "mass_update (g)") +
  theme_minimal()

# patches ~ generation
patches.gen <- ggplot(prey_details_df, aes(x = generation, y = patches)) +
  stat_summary(fun = mean, geom = "line", col = "deeppink4", linewidth = 1) +
  labs(y = "patches visited", x = "generation") +
  theme_minimal()

# number of patches ~ lv
patches.lv <- ggplot(prey_details_df, aes(x = lv, y = patches)) +
  stat_summary(fun = mean, geom = "line", col = "deeppink4", linewidth = 1) +
  labs(y = "patches visited", x = "lv") +
  theme_minimal()

# number of patches ~ speed
patches.speed <- ggplot(prey_details_df, aes(x = speed, y = patches)) +
  stat_summary(fun = mean, geom = "line", col = "deeppink4", linewidth = 1) +
  labs(y = "patches visited", x = "speed (m/s)") +
  theme_minimal()

# patches ~ offspring
patches.off <- ggplot(prey_details_df, aes(x = offspring, y = patches)) +
  stat_summary(fun = mean, geom = "line", col = "deeppink4", linewidth = 1) +
  labs(x = "offspring", y = "patches visited") + 
  theme_minimal()

# cal_net ~ generation
cal.gen <- ggplot(prey_details_df, aes(x = generation, y = cal_net)) +
  stat_summary(fun = mean, geom = "line", col = "deeppink4", linewidth = 1) +
  labs(x = "generation", y = "net calories") +
  theme_minimal()

# cal_net ~ mass_update
cal.mass <- ggplot(prey_details_df, aes(x = mass_update, y = cal_net)) +
  stat_summary(fun = mean, geom = "line", col = "deeppink4", linewidth = 1) +
  labs(x = "mass_update (g)", y = "net calories") +
  theme_minimal()

# cal_net ~ lv
cal.lv <- ggplot(prey_details_df, aes(x = lv, y = cal_net)) +
  stat_summary(fun = mean, geom = "line", col = "deeppink4", linewidth = 1) +
  labs(y = "net calories", x = "lv") +
  theme_minimal()

# cost ~ generation
cost.gen <- ggplot(prey_details_df, aes(x = generation, y = costs)) +
  stat_summary(fun = mean, geom = "line", col = "deeppink4", linewidth = 1) +
  labs(y = "metabolic costs (cal)", y = "generation") +
  theme_minimal()

# cost ~ mass_update
cost.mass <- ggplot(prey_details_df, aes(x = mass_update, y = costs)) +
  stat_summary(fun = mean, geom = "line", col = "deeppink4", linewidth = 1) +
  labs(y = "metabolic costs (cal)", y = "mass_update") +
  theme_minimal()

# cost ~ speed
cost.speed <- ggplot(prey_details_df, aes(x = speed, y = costs)) +
  stat_summary(fun = mean, geom = "line", col = "deeppink4", linewidth = 1) +
  labs(y = "metabolic costs (cal)", y = "speed (m/s)") + 
  theme_minimal()

# speed ~ generation
speed.gen <- ggplot(prey_details_df, aes(x = generation, y = speed)) +
  stat_summary(fun = mean, geom = "line", col = "deeppink4", linewidth = 1) +
  labs(y = "speed (m/s)", x = "generation") +
  theme_minimal()

# speed ~ cal_net 
speed.cal <- ggplot(prey_details_df, aes(x = cal_net, y = speed)) +
  stat_summary(fun = mean, geom = "line", col = "deeppink4", linewidth = 1) +
  labs(x = "net calories", y = "speed (m/s)") +
  theme_minimal()

# speed ~ lv
speed.lv <- ggplot(prey_details_df, aes(x = lv, y = speed)) +
  stat_summary(fun = mean, geom = "line", col = "deeppink4", linewidth = 1) +
  labs(y = "speed (m/s)", x = "lv") +
  theme_minimal()

# speed ~ mass_update
speed.mass <- ggplot(prey_details_df, aes(x = mass_update, y = speed)) +
  stat_summary(fun = mean, geom = "line", col = "deeppink4", linewidth = 1) +
  labs(y = "speed (m/s)", x = "mass_update (g)") +
  theme_minimal()

# offspring ~ generation
offspring.gen <- ggplot(prey_details_df, aes(x = generation, y = offspring)) +
  stat_summary(fun = mean, geom = "line", col = "deeppink4", linewidth = 1) +
  labs(x = "generation", y = "offspring") +
  theme_minimal()

# offspring ~ speed
offspring.speed <- ggplot(prey_details_df, aes(x = speed, y = offspring)) +
  stat_summary(fun = mean, geom = "line", col = "deeppink4", linewidth = 1) +
  labs(y = "offspring", x = "speed") +
  theme_minimal()

# offspring ~ lv
offspring.lv <- ggplot(prey_details_df, aes(x = lv, y = offspring)) +
  stat_summary(fun = mean, geom = "line", col = "deeppink4", linewidth = 1) +
  labs(y = "offspring", x = "lv") +
  theme_minimal()

# offspring ~ cal_net
offspring.cal <- ggplot(prey_details_df, aes(x = cal_net, y = offspring)) +
  stat_summary(fun = mean, geom = "line", col = "deeppink4", linewidth = 1) +
  labs(y = "offspring", x = "net calories") +
  theme_minimal()

# offspring ~ mass_update
offspring.mass <- ggplot(prey_details_df, aes(x = mass_update, y = offspring)) +
  stat_summary(fun = mean, geom = "line", col = "deeppink4", linewidth = 1) +
  labs(y = "offspring", x = "mass_update (g)") +
  theme_minimal()

# mass ~ gen
mass.gen <- ggplot(prey_details_df, aes(x = generation, y = mass_update)) +
  stat_summary(fun = mean, geom = "line", col = "deeppink4", linewidth = 1) +
  labs(y = "mass_updated (g)", x = "generation") +
  theme_minimal()


plots <- list(rel.lv.gen,
              taup.gen,
              tauv.gen,
              sig.gen,
              patches.gen,
              cost.gen,
              offspring.gen,
              cal.gen,
              mass.gen,
              speed.gen,
              patches.lv,
              speed.lv,
              patches.speed,
              cal.lv,
              offspring.lv,
              cost.speed,
              offspring.speed,
              cal.mass,
              lv.mass,
              speed.mass,
              speed.cal,
              offspring.mass,
              cost.mass,
              offspring.cal,
              patches.off)

final.plot <- wrap_plots(plots, ncol = 5)
final <- final.plot + plot_annotation('1000g, 10 calories')
print(final)
ggsave("H:/GitHub/ballistic-movement/sim_results/constant_resources_nopred_varycalories/figures/1000g_10cal_finalplot.PNG", plot = final, width = 15, height = 8, dpi = 800)


