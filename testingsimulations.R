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
library(deSolve)

# Source the functions (ensure 'functions.R' is available in the working directory)
source("functions.R")

#----------------------------------------------------------------------
# Set up the global parameters for the simulation
#----------------------------------------------------------------------

# Predator mass (g)
mass_pred <- 5000

# Prey mass (g)
mass_prey <- 1000


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
prey_offspring <- prey.fitness.deb(benefits = total_resources2, mass = mass_prey, costs = 0, 
                 models = mod4, 
                 crossings = 20, 
                 calories = 10, 
                 alpha = 0.25,
                 beta = 0.75,
                 kap = 0.5,
                 metric = "offspring")
print(prey_offspring)
prey_lifespan <- prey.fitness.deb(benefits = total_resources2, mass = mass_prey, costs = 0, 
                                  models = mod4, 
                                  crossings = 20, 
                                  calories = 10, 
                                  alpha = 0.25,
                                  beta = 0.75,
                                  kap = 0.5,
                                  metric = "lifespan")
print(prey_lifespan)

#----------------------------------------------------------------------
# trying loops
#----------------------------------------------------------------------

#number of generations
GENS <- 1

#number of 'arenas' 
REPS <- 5

#number of prey
n_prey <- 2

#sampling interval
t <- sampling2.0(mass_prey, crossings = 30)

#Empty lists for storing the results of the current generation
prey <- vector("list", GENS)

for(G in 1:GENS) {
  prey[[G]] <- vector("list", REPS)
  
  for(R in 1:REPS){
    
    #Generate the prey movement models
    #If the first gen, generate movement parameters from the mass functions
    if(G == 1){
      #Generate the HR centres of the prey
      CENTRES <- rbvpois(n = n_prey,
                         a = pred.SIG(mass_pred)*.75,
                         b = pred.SIG(mass_pred)*.75,
                         c = 0)
      CENTRES <- scale(CENTRES, scale = FALSE)
      
      PREY_mods <- vector("list", n_prey)
      for(i in 1:n_prey){
        # Prey movement parameters
        prey_tau_p <- prey.tau_p(mass_prey, variance = TRUE)
        prey_tau_v <- prey.tau_v(mass_prey, variance = TRUE)
        prey_sig <- prey.SIG(mass_prey)
        prey_lv <- sqrt((prey_tau_v/prey_tau_p)*prey_sig)
        
        PREY_mods[[i]] <- ctmm(tau = c(prey_tau_p,prey_tau_v),
                               mu = c(CENTRES[i,1],CENTRES[i,2]),
                               sigma = prey_sig)
      } #Closes loop over n_prey
      
      #simulate prey movement
      prey_tracks <- vector("list", n_prey)
      for(i in 1:n_prey){
        prey_tracks[[i]] <- simulate(PREY_mods[[i]],t = t)
      } #closes loop over prey_tracks
      
      #calculate prey benefits
      benefits_prey <- vector()
      for(i in 1:n_prey){
        benefits_prey[i] <- grazing(prey_tracks[[i]], FOOD)
      }
      
      #Calculate prey fitness
      offspring_prey <- prey.fitness.deb(benefits = total_resources2, 
                                         mass = mass_prey, 
                                         costs = 0, 
                                         models = PREY_mods, 
                                         crossings = 30, 
                                         calories = 10, 
                                         alpha = 0.25,
                                         beta = 0.75,
                                         kap = 0.5,
                                         metric = "offspring")
      
      #Get the values of the prey movement model parameters
      prey_lvs <- numeric(n_prey)
      prey_TAU_V <- numeric(n_prey)
      prey_TAU_P <- numeric(n_prey)
      prey_SIGMA <- numeric(n_prey)
      prey_SPEED <- numeric(n_prey)
      for(i in 1:n_prey){
        prey_TAU_V[i] <- PREY_mods[[i]]$tau["velocity"]
        prey_TAU_P[i] <- PREY_mods[[i]]$tau["position"]
        prey_SIGMA[i] <- ctmm:::area.covm(PREY_mods[[i]]$sigma)
        prey_SPEED[i] <- if(nrow(summary(PREY_mods[[i]], units = FALSE)$CI)==4){summary(PREY_mods[[i]], units = FALSE)$CI[4,2]} else{Inf}
        prey_lvs[i] <- sqrt((prey_TAU_V[i]/prey_TAU_P[i])*prey_SIGMA[i])
      }
      
      #Summarise the results of the prey 
      prey[[G]][[R]] <- data.frame(generation = G,
                              replicate = R,
                              tau_p = prey_TAU_P,
                              tau_v = prey_TAU_V,
                              sig = prey_SIGMA,
                              speed = prey_SPEED,
                              lv = prey_lvs,
                              patches = benefits_prey,
                              offspring = offspring_prey)
    }
  }
}

results <- do.call(
  rbind,
  lapply(prey, function(gen) {
    do.call(rbind, lapply(gen, function(rep) {
      if (is.data.frame(rep)) rep else NULL
    }))
  })
)

head(results)
print(summary(results))


