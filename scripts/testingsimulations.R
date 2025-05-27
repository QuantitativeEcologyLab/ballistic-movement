
# Preamble ----

# Set the working directory
setwd("~/ballisticmovement/ballistic-movement/scripts")
# Set the random seed
set.seed(1)

# Import necessary packages
library(extraDistr)
library(ctmm)
library(terra)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(patchwork)

# Source the functions (ensure 'functions.R' is available in the working directory)
source("functions.R")

#----------------------------------------------------------------------
# set up global parameters----
#----------------------------------------------------------------------

#predator mass (g)
mass_pred <- 20000

# Prey mass (g)
mass_prey <- 5000

#set sampling interval and lifespan
t <- sampling(mass_prey)
interval <- sampling(mass_prey, metric = "interval")

#number of individuals in arena
n_prey <- 5

#number of arenas
REPS <- 2

#number of generations
GENS <- 500

#build food raster
FOOD <- createFoodRaster(mass_prey, patch_width = 35)
plot(FOOD, col = "steelblue")
grid(nx = ncol(FOOD), ny = nrow(FOOD), col = "black", lty = "dotted")

#lists for storing results
prey_res <- list()
prey_details <- list()

#----------------------------------------------------------------------
# run the simulation----
#----------------------------------------------------------------------

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
      patches[[i]] <- grazing(PREY_tracks[[i]], FOOD, metric = "patches")
    }
    
    #extract speed from model
    speed <- numeric(n_prey)
    for(i in 1:n_prey){
      speed[[i]] <- speed_val(models = PREY_mods[[i]])
    }
    
    #assign net calories to each individual
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
      # Defensive: check result is valid
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
                            cal = cal_net,
                            cal_max = cal_max,
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
      
      save(prey_res, file = '~/ballisticmovement/ballistic-movement/sim_results/lv_Evo_5000g_prey_res2.Rda')
      save(prey_details, file = '~/ballisticmovement/ballistic-movement/sim_results/lv_Evo_5000g_prey_details2.Rda')
      
      break
    }
  }
  
  #Save prey results
  save(prey_res, file = '~/ballisticmovement/ballistic-movement/sim_results/lv_Evo_5000g_prey_res2.Rda')
  save(prey_details, file = '~/ballisticmovement/ballistic-movement/sim_results/lv_Evo_5000g_prey_details2.Rda')
  
  #progress report
  print(G)
}

print(prey_res)
print(prey_details)

#----------------------------------------------------------------------
# make diagnostic figures----
#----------------------------------------------------------------------

#make data sets compatible
prey_res_df <- do.call(rbind, prey_res)
print(prey_res_df)
prey_details_df <- do.call(rbind, prey_details)
print(prey_details_df)

# relative change in lv ~ generation
PREY_LV <- prey_res_df$lv[1]
prey_res_df$rel.lv <- prey_res_df$lv/PREY_LV
prey_res_df$rel_var <- prey_res_df$var / (PREY_LV^2)
prey_res_df$rel_sd <- sqrt(prey_res_df$rel_var)

rel.lv.gen <- ggplot(prey_res_df, aes(x = generation, y = rel.lv)) +
  geom_line(color = "red", linewidth = 1) +
  geom_hline(yintercept = 1, color = 'grey30', linetype = "dashed") +
  geom_ribbon(aes(ymin = rel.lv - rel_sd,
                  ymax = rel.lv + rel_sd),
              fill = "red", alpha = 0.3) +
  labs(
    y = "relative change in lv",
    x = "generation") +
  theme_minimal()

print(rel.lv.gen)

# calories ~ lv
cal.lv <- ggplot(prey_details_df, aes(x = lv, y = cal)) +
  geom_line(color = "red", linewidth = 1) +
  labs(
    y = "net calories",
    x = "lv") +
  theme_minimal()

print(cal.lv)

# offspring ~ lv
offspring.lv <- ggplot(prey_details_df, aes(x = lv, y = offspring)) +
  geom_line(color = "red", linewidth = 1) +
  labs(
    y = "offspring",
    x = "lv") +
  theme_minimal()

print(offspring.lv)

# mass ~ gen
mass.gen <- ggplot(prey_details_df, aes(x = generation, y = mass_update)) +
  geom_line(color = "red", linewidth = 1) +
  labs(
    y = "mass_updated",
    x = "generation") +
  theme_minimal()

print(mass.gen)

# number of patches ~ lv
patches.lv <- ggplot(prey_details_df, aes(x = lv, y = patches)) +
  geom_line(color = "red", linewidth = 1) +
  labs(
    y = "number of patches",
    x = "lv") +
  theme_minimal()

print(patches.lv)

# speed ~ lv
speed.lv <- ggplot(prey_details_df, aes(x = lv, y = speed)) +
  geom_line(color = "red", linewidth = 1) +
  labs(
    y = "speed",
    x = "lv") +
  theme_minimal()

print(speed.lv)

# offspring ~ calories
offspring.cal <- ggplot(prey_details_df, aes(x = cal, y = offspring)) +
  geom_line(color = "red", linewidth = 1) +
  labs(
    y = "offspring",
    x = "calories") +
  theme_minimal()

print(offspring.cal)

# offspring ~ mass
offspring.mass <- ggplot(prey_details_df, aes(x = mass_update, y = offspring)) +
  geom_line(color = "red", linewidth = 1) +
  labs(
    y = "offspring",
    x = "mass_updated") +
  theme_minimal()

print(offspring.mass)

# speed ~ mass
speed.mass <- ggplot(prey_details_df, aes(x = mass_update, y = speed)) +
  geom_line(color = "red", linewidth = 1) +
  labs(
    y = "speed",
    x = "mass_updated") +
  theme_minimal()

print(speed.mass)

# offspring ~ speed
offspring.speed <- ggplot(prey_details_df, aes(x = speed, y = offspring)) +
  geom_line(color = "red", linewidth = 1) +
  labs(
    y = "offspring",
    x = "speed") +
  theme_minimal()

print(offspring.speed)

# number of patches ~ speed
patches.speed <- ggplot(prey_details_df, aes(x = speed, y = patches)) +
  geom_line(color = "red", linewidth = 1) +
  labs(
    y = "patches",
    x = "speed") +
  theme_minimal()

print(patches.speed)

# speed ~ generation
speed.gen <- ggplot(prey_details_df, aes(x = generation, y = speed)) +
  stat_summary(fun = mean, geom = "line", col = "red", linewidth = 1)+
  labs(
    y = "speed",
    x = "generation") +
  theme_minimal()

print(speed.gen)

# tau_v ~ generation
tauv.gen <- ggplot(prey_details_df, aes(x = generation, y = tau_v)) +
  stat_summary(fun = mean, geom = "line", col = "red", linewidth = 1)+
  labs(
    y = "tau_v",
    x = "generation") +
  theme_minimal()

print(tauv.gen)

# tau_p ~ generation
taup.gen <- ggplot(prey_details_df, aes(x = generation, y = tau_p)) +
  stat_summary(fun = mean, geom = "line", col = "red", linewidth = 1)+
  labs(
    y = "tau_p",
    x = "generation") +
  theme_minimal()

print(taup.gen)

# movement paths
track_df <- do.call(rbind, lapply(1:length(PREY_tracks), function(i) {
  data.frame(x = PREY_tracks[[i]]$x,
             y = PREY_tracks[[i]]$y,
             t = PREY_tracks[[i]]$t,
             id = as.factor(i))
}))

food_df <- as.data.frame(FOOD, xy = TRUE)
colnames(food_df) <- c("x", "y", "value")

move <- ggplot() +
  geom_raster(data = food_df, aes(x = x, y = y, fill = value)) +
  geom_path(data = track_df, aes(x = x, y = y, group = id, color = id), inherit.aes = FALSE) +
  theme_minimal() +
  labs(title = "Prey Movement over Food Landscape", x = "X", y = "Y") +
  theme(legend.position = "none")
print(move)

#final plot
plots <- list(rel.lv.gen,
              cal.lv, 
              offspring.lv, 
              patches.lv, 
              patches.speed,
              speed.lv, 
              speed.mass,
              speed.gen,
              offspring.cal, 
              offspring.mass,
              offspring.speed, 
              mass.gen,
              tauv.gen, 
              taup.gen)

final.plot <- wrap_plots(plots, ncol = 3)
print(final.plot)

#prey lv and fitness
# Prepare filtered data (offspring != 0)
DATA <- prey_details_df[which(prey_details_df$offspring != 0), ]

# Round lv for aggregation
prey_details_df$lv2 <- round(prey_details_df$lv, 2)
# Aggregate offspring mean by lv2
AGG <- aggregate(offspring ~ lv2, data = prey_details_df, FUN = mean)

# Check if offspring varies enough to fit model
if(length(unique(AGG$offspring)) > 1) {
  FIT <- nls(offspring ~ a * (lv2 + c) * exp(-b * (lv2 + c)),
             start = list(a = 6, b = 0.1, c = -2),
             data = AGG)
  
  ricker <- function(x) {
    coef(FIT)[1] * (x + coef(FIT)[3]) * exp(-coef(FIT)[2] * (x + coef(FIT)[3]))
  }
  
  x <- seq(min(prey_details_df$lv2), max(prey_details_df$lv2), 0.01)
  y <- ricker(x)
  
  plot(offspring ~ lv,
       data = DATA,
       pch = 20,
       col = adjustcolor("grey70", alpha = 0.4),
       xlab = "Ballistic lengthscale (m)",
       ylab = "Prey fitness (offspring)")
  
  points(offspring ~ lv2,
         data = AGG,
         pch = 20,
         col = adjustcolor("blue", alpha = 0.8))
  
  # Add fitted curve only if FIT exists
  lines(x, y, col = "steelblue", lwd = 2)
  
} else {
  # No variation: just plot raw and aggregated data
  plot(offspring ~ lv,
       data = DATA,
       pch = 20,
       col = adjustcolor("grey70", alpha = 0.4),
       xlab = "Ballistic lengthscale (m)",
       ylab = "Prey fitness (offspring)")
  
  points(offspring ~ lv2,
         data = AGG,
         pch = 20,
         col = adjustcolor("blue", alpha = 0.8))
  
  message("Skipped nonlinear fit: insufficient variation in offspring")
}

# Add mean and variance shading in both cases
mu <- mean(DATA$lv, na.rm = TRUE)
sig <- var(DATA$lv, na.rm = TRUE)

abline(v = mu, col = '#046C9A', lty = 'dashed')
rect(xleft = mu - sig, xright = mu + sig,
     ybottom = par("usr")[3], ytop = par("usr")[4], 
     border = NA,
     col = adjustcolor("#046C9A", alpha = 0.3))

