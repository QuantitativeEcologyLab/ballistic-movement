#----------------------------------------------------------------------
# Preamble
#----------------------------------------------------------------------

# Set the working directory
setwd("~/ballistic_movement_models")
# Set the random seed
set.seed(1)

# Import necessary packages
library(extraDistr)
library(ctmm)
library(terra)
library(ggplot2)
library(dplyr)

# Source the functions (ensure 'functions.R' is available in the working directory)
source("functions.R")

#----------------------------------------------------------------------
# Set up the global parameters for the simulation
#----------------------------------------------------------------------

# Predator mass (g)
mass_pred <- 5000

# Prey mass (g)
mass_prey <- prey.mass(mass_pred)

#----------------------------------------------------------------------
# testing prey.mod function
#----------------------------------------------------------------------

#generate movement model

mod <- prey.mod(mass_prey, variance = FALSE)

#define the sampling duration/schedule

t <- seq(0,5 %#% 'month', 1 %#% 'hr')

#simulate model

SIM <- simulate(mod,t = t)

#fit model to it

GUESS <- ctmm.guess(SIM, interactive = FALSE)

FIT <- ctmm.select(SIM, GUESS, cores = -1)

summary(FIT)

# home range estimate

HR <- akde(SIM, FIT)

#visualize the results

plot(SIM, UD = HR)

#wahoo prey.mod creates movement data with ctmm

#----------------------------------------------------------------------
# testing movement parameters based on mass function(s)
#----------------------------------------------------------------------

#prey movement parameters
  prey_tau_p <- prey.tau_p(mass_prey, variance = FALSE)
  prey_tau_v <- prey.tau_v(mass_prey, variance = FALSE)
  prey_sig <- prey.SIG(mass_prey)
  prey_lv <- sqrt((prey_tau_v/prey_tau_p)*prey_sig)
  
#generate model
  mod2 <- ctmm(tau = c(prey_tau_p,prey_tau_v),
                         mu = c(0,0),
                         sigma = prey_sig)

#define the sampling duration/schedule

t2 <- seq(0,5 %#% 'month', 1 %#% 'hr')

#simulate model

SIM2 <- simulate(mod2,t = t2)

#fit model to it

GUESS2 <- ctmm.guess(SIM2, interactive = FALSE)

FIT2 <- ctmm.select(SIM2, GUESS2, cores = -1)

summary(FIT2)

HR2 <- akde(SIM2, FIT2)

#visualize the results

plot(SIM2, UD = HR2)

#----------------------------------------------------------------------
# layer the food raster
#----------------------------------------------------------------------

#prey movement parameters
prey_tau_p <- prey.tau_p(mass_prey, variance = FALSE)
prey_tau_v <- prey.tau_v(mass_prey, variance = FALSE)
prey_sig <- prey.SIG(mass_prey)
prey_lv <- sqrt((prey_tau_v/prey_tau_p)*prey_sig)

#generate raster
FOOD <- patches(mass = mass_prey, width = round(sqrt(prey.SIG(mass_prey))/10), pred = FALSE, type = "random")
plot(FOOD, main = "Resource Landscape")

#generate model
mod3 <- ctmm(tau = c(prey_tau_p, prey_tau_v),
             mu = c(0,0),
             sigma = prey_sig)

#simulate movement
t3 <- seq(0, 5 %#% "month", 1 %#% "hour")
SIM3 <- simulate(mod3, t = t3)

#extract raster values along path (not needed?)
track_df <- data.frame(x = SIM3$x, y = SIM3$y)
track_prey <- vect(track_df, geom = c("x", "y"))
res_vals <- extract(FOOD, track_prey)[, 2]
total_resources <- sum(res_vals, na.rm = TRUE)
cat("total resource encountered by prey:", total_resources, "\n")

#use grazing function
patch_visits <- grazing(track_df, FOOD, metric = "patches")
print(patch_visits)
time_between <- grazing(track_df, FOOD, metric = "time")
print(time_between)

#fit and visualize

GUESS3 <- ctmm.guess(SIM3, interactive = FALSE)
FIT3 <- ctmm.select(SIM3, GUESS3, cores = -1)

HR3 <- akde(SIM3, FIT3)

plot(SIM3, UD = HR3)



#----------------------------------------------------------------------
# calculate fitness 
#----------------------------------------------------------------------

#prey movement parameters
prey_tau_p <- prey.tau_p(mass_prey, variance = FALSE)
prey_tau_v <- prey.tau_v(mass_prey, variance = FALSE)
prey_sig <- prey.SIG(mass_prey)
prey_lv <- sqrt((prey_tau_v/prey_tau_p)*prey_sig)

#generate raster
FOOD2 <- patches(mass = mass_prey, width = 20, pred = FALSE, type = "uniform")
plot(FOOD2, main = "Resource Landscape")

#generate model
mod4 <- ctmm(tau = c(prey_tau_p, prey_tau_v),
             mu = c(0,0),
             sigma = prey_sig)

#simulate movement
t4 <- seq(0, 5 %#% "month", 1 %#% "hour")
SIM4 <- simulate(mod4, t = t4)

#extract raster values along path
track_df2 <- data.frame(x = SIM4$x, y = SIM4$y)
track_prey2 <- vect(track_df2, geom = c("x", "y"))
res_vals2 <- extract(FOOD, track_prey2)[, 2]
total_resources2 <- sum(res_vals2, na.rm = TRUE)
cat("total resource encountered by prey:", total_resources, "\n")

#use grazing function
patch_visits2 <- grazing(track_df2, FOOD, metric = "patches")
print(patch_visits2)
time_between2 <- grazing(track_df2, FOOD, metric = "time")
print(time_between2)

#fit and visualize

GUESS4 <- ctmm.guess(SIM4, interactive = FALSE)
FIT4 <- ctmm.select(SIM4, GUESS4, cores = -1)

HR4 <- akde(SIM4, FIT4)

plot(SIM4, UD = HR4)

#calculate fitness and lifespan
prey_offspring <- prey.fitness.deb(benefits = patch_visits2, 
                                   mass = mass_prey, 
                                   costs = 0,
                                   models = mod4,
                                   crossings = 20,
                                   calories = 10,
                                   risk_factor = 0,
                                   constant = 1,
                                   DEB = FALSE,
                                   metric = "offspring"
                                   )
print(prey_offspring)

prey_offspring.deb <- prey.fitness.deb(benefits = patch_visits2, 
                                   mass = mass_prey, 
                                   costs = 0,
                                   models = mod4,
                                   crossings = 20,
                                   calories = 10,
                                   risk_factor = 0,
                                   constant = 1,
                                   DEB = TRUE)
print(prey_offspring.deb)

#----------------------------------------------------------------------
# trying loops
#----------------------------------------------------------------------

#set sampling interval and lifespan
t <- sampling(mass_prey, crossings = 20)

#energetic value of a patch
CALS <- ((10^(0.774 + 0.727*log10(mass_prey)))^1.22)/150

#number of individuals in arena
n_prey <- 10

#number of arenas
REPS <- 5

#number of generations
GENS <- 10

#build food raster
FOOD <- patches(mass_prey, width = 20, pred = FALSE, type = "uniform")

#lists for storing results
prey_res <- list()
prey_details <- list()

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
    
    benefits_prey <- vector()
    for(i in 1:n_prey){
      benefits_prey[i] <- grazing(PREY_tracks[[i]], FOOD)
    }
    
    f <- vector()
    for(i in 1:n_prey){
      f[i] <- calculate_f(benefits_prey)
    }
    
    offspring_prey <- prey.fitness.debkiss(mass = mass_prey,
                                           f = f,
                                           costs = 0,
                                           crossings = 20,
                                           calories = 10,
                                           DEBkiss = TRUE,
                                           models = PREY_mods,
                                           metric = "offspring")
    
    #get values
    prey_lvs <- vector()
    prey_TAU_V <- vector()
    prey_TAU_P <- vector()
    prey_SIGMA <- vector()
    prey_SPEED <-vector()
    for(i in 1:n_prey){
      prey_TAU_V[i] <- PREY_mods[[i]]$tau["velocity"]
      prey_TAU_P[i] <- PREY_mods[[i]]$tau["position"]
      prey_SIGMA[i] <- ctmm:::area.covm(PREY_mods[[i]]$sigma)
      prey_SPEED[i] <- if(nrow(summary(PREY_mods[[i]], units = FALSE)$CI)==4){summary(PREY_mods[[i]], units = FALSE)$CI[4,2]} else{Inf}
      prey_lvs[i] <- sqrt((prey_TAU_V[i]/prey_TAU_P[i])*prey_SIGMA[i])
    }
    
    #summarise
    prey[[R]] <- data.frame(generation = G,
                            tau_p = prey_TAU_P,
                            tau_v = prey_TAU_V,
                            sig = prey_SIGMA,
                            speed = prey_SPEED,
                            lv = prey_lvs,
                            patches = benefits_prey,
                            encounters = 0, 
                            offspring = offspring_prey)
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
  
  print(G)
}

print(prey_res)

#model change in lv over generations
res_df <- do.call(rbind, prey_res)

res_df$rel_change_lv <- res_df$lv / res_df$lv[1]

plot(res_df$generation, res_df$rel_change_lv, type = "b", pch = 19,
     xlab = "gen", ylab = "rel change in lv")
abline(h = 1, lty = 2, col = "gray")


##ggplot version
res_df <- do.call(rbind, lapply(seq_along(prey_res), function(i){
  df <- prey_res[[i]]
  df$sim <- i
  df
}))

summary_df <- res_df %>%
  group_by(generation) %>%
  summarize(
    lv_mean = mean(lv),
    lv_sd = sd(lv),
    .groups = "drop"
  ) %>%
  mutate(
    lv_se = lv_sd / sqrt(n_prey),
    lv_ci95 = 1.96 * lv_se,
    rel_change_lv = lv_mean / lv_mean[1],
    lower = (lv_mean - lv_sd) / lv_mean[1],
    upper = (lv_mean + lv_sd) / lv_mean[1]
  )


ggplot(summary_df, aes(x = generation, y = rel_change_lv)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "red", alpha = 0.2) +
  geom_line(color = "red", size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  labs(
    x = "Generation",
    y = "Relative change in lv"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line()
  )

#plot the movement paths
# Flatten list of tracks into data frame
track_df <- do.call(rbind, lapply(1:length(PREY_tracks), function(i) {
  data.frame(x = PREY_tracks[[i]]$x,
             y = PREY_tracks[[i]]$y,
             t = PREY_tracks[[i]]$t,
             id = as.factor(i))
}))

ggplot(track_df, aes(x = x, y = y, group = id, color = id)) +
  geom_path() +
  theme_minimal() +
  labs(title = "Prey Movement Tracks", x = "X", y = "Y") +
  theme(legend.position = "none")




