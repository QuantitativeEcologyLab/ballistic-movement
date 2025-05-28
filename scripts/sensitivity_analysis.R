
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
  
  #generate sampling interval 
  t <- sampling(mass_prey)
  interval <- sampling(mass_prey, metric = "interval")
  
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
    save(prey_res, file = paste0(save_prefix, "_res.Rda")) #delete this line
    save(prey_details, file = paste0(save_prefix, "_details.Rda"))
    
    cat("Completed generation", G, "\n")
  }
  return(prey_details)
}


#-------------------------------------------------------------------------------
# create sensitivity analysis----
#-------------------------------------------------------------------------------

#------------------------------#
#    Static Simulation Params  #
#------------------------------#

#number of individuals in arena
n_prey <- 5

#number of arenas
REPS <- 3

#number of generations
GENS <- 5


# Define your parameter combinations here
param_grid <- expand.grid(
  mass_prey = seq(2000, 10000, 2000),
  patch_width = seq(10, 60, 10),
  biomass = seq(1, 15, 3),
  calories = seq(1000, 4000, 500) 
)

# Prepare a list to store all results
all_results <- list()

for(i in seq_len(nrow(param_grid))) {
tic(paste("Scenario", i))
  params <- param_grid[i, ]
  
  cat("\n--- Running scenario", i, "of", nrow(param_grid), "---\n")
  cat(sprintf("Width: %d, Biomass: %d, Calories: %d, Mass: %d\n", params$patch_width, params$biomass, params$calories, params$mass_prey))
  
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

save(all_results, file = "~/ballisticmovement/ballistic-movement/sim_results/results_df.Rda")

all_results_df <- bind_rows(all_results)


# --- Summary Tables ----

# 1. Main summary table: mean offspring by mass, calories, biomass, patch_width
summary_table <- all_results_df %>%
  group_by(mass, calories, biomass, patch_width) %>%
  summarise(
    mean_offspring = mean(offspring, na.rm = TRUE),
    sd_offspring = sd(offspring, na.rm = TRUE),
    n = n(),
    prop_extinct = mean(offspring == 0, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_table)

# 2. Overview: mean offspring by mass (collapsed across others)
mass_effects <- all_results_df %>%
  group_by(mass) %>%
  summarise(
    mean_offspring = mean(offspring, na.rm = TRUE),
    sd_offspring = sd(offspring, na.rm = TRUE),
    prop_extinct = mean(offspring == 0, na.rm = TRUE),
    .groups = "drop"
  )

print(mass_effects)

# 3. Plot: Mean offspring by mass (with error bars)
mass.off <- ggplot(mass_effects, aes(x = mass, y = mean_offspring)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_offspring - sd_offspring,
                    ymax = mean_offspring + sd_offspring), width = 100) +
  theme_minimal() +
  labs(title = "Mean Offspring by Mass",
       y = "Mean Offspring", x = "Initial Mass")
print(mass.off)

# 4. Heatmap (3D slicing): Calories vs Patch Width per Mass bin
heat <- ggplot(all_results_df, aes(x = factor(patch_width), y = factor(calories), fill = offspring)) +
  geom_tile() +
  facet_wrap(~ mass, ncol = 2) +
  scale_fill_viridis_c() +
  labs(title = "Offspring by Patch Width and Calories per Mass",
       x = "Patch Width", y = "Calories") +
  theme_minimal()

# 5. Optional: 4D interaction plot using facets (collapsed by biomass)
interaction <- ggplot(all_results_df, aes(x = biomass, y = offspring)) +
  geom_jitter(alpha = 0.3) +
  geom_smooth(method = "loess", se = FALSE) +
  facet_grid(mass ~ calories, labeller = label_both) +
  theme_minimal() +
  labs(title = "Offspring by Biomass across Mass & Calories",
       y = "Offspring", x = "Biomass")
print(interaction)

# 6. run linear model to data
lm_model <- lm(offspring ~ patch_width + biomass + calories + mass, data = all_results_df)
summary(lm_model)

lm <- ggplot(all_results_df, aes(x = calories, y = offspring, color = as.factor(mass))) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm") +
  facet_grid(biomass ~ patch_width) +
  theme_minimal()

print(lm)

sens.plots <- list(mass.off,
              heat,
              interaction,
              lm)

final.sens.plot <- wrap_plots(sens.plots, ncol = 2)
print(final.sens.plot)

# Random Forest for non-linear effects (requires 'randomForest' package)
library(randomForest)
rf_fit <- randomForest(offspring ~ patch_width + biomass + calories, data = all_results_df)
print(rf_fit)

# more analysis of sensitivity

# Load required packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(car)
library(patchwork)  # for combining plots
library(MASS)       # for stepwise selection

# Basic checks on multicollinearity and model quality
lm_model <- lm(offspring ~ patch_width + biomass + calories + mass, data = all_results_df)
summary(lm_model)
car::vif(lm_model)  # check variance inflation factors

# Interaction model
lm_int <- lm(offspring ~ (patch_width + biomass + calories) * mass, data = all_results_df)
summary(lm_int)

# Stepwise model selection to reduce complexity
step_model <- MASS::stepAIC(lm_int, direction = "both")
summary(step_model)

# Plot response surfaces: calories vs patch_width for various masses
mass_vals <- unique(all_results_df$mass)
mass_vals <- sort(mass_vals)

plots <- lapply(mass_vals, function(m) {
  df_sub <- all_results_df %>% filter(mass == m)
  ggplot(df_sub, aes(x = calories, y = patch_width, fill = offspring)) +
    geom_tile() +
    scale_fill_viridis_c() +
    labs(title = paste("Mass:", m), x = "Calories", y = "Patch Width", fill = "Offspring") +
    theme_minimal()
})

# Combine all plots into one grid for comparison
wrap_plots(plots)

# Aggregate offspring by scenario for robustness metric
robust_summary <- all_results_df %>%
  group_by(patch_width, biomass, calories) %>%
  summarise(
    avg_offspring = mean(offspring),
    min_offspring = min(offspring),
    max_offspring = max(offspring),
    prop_success = mean(offspring > 0),
    .groups = 'drop'
  ) %>%
  arrange(desc(prop_success), desc(avg_offspring))

# Plot robustness heatmap
robust_plot <- ggplot(robust_summary, aes(x = calories, y = patch_width, fill = prop_success)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Proportion of Successful Simulations", x = "Calories", y = "Patch Width", fill = "Prop. Success") +
  theme_minimal()

robust_plot

# Optional: Export the top parameter sets
best_sets <- robust_summary %>% filter(prop_success == 1)
write.csv(best_sets, "best_parameter_sets.csv", row.names = FALSE)




