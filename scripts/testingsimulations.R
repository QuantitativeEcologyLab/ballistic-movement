#----------------------------------------------------------------------
# Preamble
#----------------------------------------------------------------------

# Set the working directory
setwd("~scripts")
# Set the random seed
set.seed(1)

# Import necessary packages
library(extraDistr)
library(ctmm)
library(terra)
library(ggplot2)
library(dplyr)

# Source the functions (ensure 'functions.R' is available in the working directory)
source("scripts/functions.R")

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
print(SIM3)

#use grazing function
patch_visits <- grazing(SIM3, FOOD, metric = "patches")
print(patch_visits)
time_between <- grazing(SIM3, FOOD, metric = "time")
print(time_between)

#fit and visualize

GUESS3 <- ctmm.guess(SIM3, interactive = FALSE)
FIT3 <- ctmm.select(SIM3, GUESS3, cores = -1)

HR3 <- akde(SIM3, FIT3)

plot(SIM3, UD = HR3)



#----------------------------------------------------------------------
# troubleshooting 
#----------------------------------------------------------------------

#patches and grazing

mass3 <- 500

habitat <- patches(mass = mass3, type = "random", calories = 10, width = 10)

#generate model
mod3 <- ctmm(tau = c(prey_tau_p,prey_tau_v),
             mu = c(0,0),
             sigma = prey_sig)

#define the sampling duration/schedule

t3 <- seq(0,5 %#% 'month', 1 %#% 'hr')

#simulate model

track <- simulate(mod3,t = t3)

patch_count <- grazing(track, habitat, metric = "patches")
calories_consumed <- grazing(track, habitat, metric = "calories")
time <- grazing(track, habitat, metric = "time")

print(patch_count)
print(calories_consumed)
print(time)

plot(habitat)
lines(track$x, track$y, col = "gray", lwd = 2)

#test reproduction and benefits 
preymasstest <- 50000

testt <- sampling2.0(preymasstest, metric = "t")

CALS <- ((10^(0.774 + 0.727*log10(preymasstest)))^1.22)/150

foodtest <- patches(preymasstest, width = 20, pred = FALSE, type = "uniform", calories = 20)

taup <- prey.tau_p(preymasstest)
tauv <- prey.tau_v(preymasstest)
sig <- prey.SIG(preymasstest)
lv <- sqrt((tauv/taup)*sig)

testmod <- ctmm(tau = c(taup, tauv),
                      mu = c(0,0),
                      sigma = sig)

testtrack <- replicate(100, ctmm::simulate(testmod, t = testt), simplify = FALSE)
print(testtrack)

plot(foodtest)
lines(testtrack$x, testtrack$y, col = "red", lwd = 2)

maxcal <- est.max.benefit(testtrack, foodtest, metric = "calories")
print(maxcal)

graze <- grazing(track, foodtest, metric = "calories")

f <- calculate_f(benefits = graze, saturation = maxcal$saturation)
print(f)

fitness <- prey.fitness.debkiss(mass = preymasstest,
                                f = f)
print(fitness)

#----------------------------------------------------------------------
# trying loops
#----------------------------------------------------------------------

# Predator mass (g)
mass_pred <- 5000

# Prey mass (g)
mass_prey <- prey.mass(mass_pred)

#set sampling interval and lifespan
t <- sampling2.0(mass_prey)

#energetic value of a patch
CALS <- ((10^(0.774 + 0.727*log10(mass_prey)))^1.22)/150

#number of individuals in arena
n_prey <- 5

#number of arenas
REPS <- 2

#number of generations
GENS <- 200

#build food raster
FOOD <- patches(mass_prey, width = 20, pred = FALSE, type = "uniform", calories = CALS)

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
      benefits_prey[i] <- grazing(PREY_tracks[[i]], FOOD, metric = "calories")
    }
    
    maxben <- as.list(est.max.benefit(tracks = PREY_tracks, habitat = FOOD, metric = "calories"))
    
    f <- vector()
    for(i in 1:n_prey){
      f[i] <- calculate_f(benefits_prey[[i]], saturation = maxben$saturation)
    }
    
    offspring_prey <- prey.fitness.debkiss(mass = mass_prey,
                                           f = f)
    
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

prey_details[[G]]$offspring
sum(prey_details[[G]]$offspring)  # total offspring in gen G
mean(prey_details[[G]]$offspring) # average offspring per individual in gen G

#model change in lv over generations
res_df <- do.call(rbind, prey_res)

res_df$rel_change_lv <- res_df$lv / res_df$lv[1]

plot(res_df$generation, res_df$rel_change_lv, type = "b", pch = 19,
     xlab = "gen", ylab = "rel change in lv")
abline(h = 1, lty = 2, col = "gray")

##ggplot version w/ variance ribbon
lv0 <- res_df$lv[1]
res_df$rel_var <- res_df$var / (lv0^2)
res_df$rel_sd <- sqrt(res_df$rel_var)

ggplot(res_df, aes(x = generation, y = rel_change_lv)) +
  geom_ribbon(aes(ymin = rel_change_lv - rel_sd,
                  ymax = rel_change_lv + rel_sd),
              fill = "red", alpha = 0.3) + 
  geom_line(color = "red", linewidth = 1) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "gray")+
  labs(x = "Generation", y = "Relative change in lv") +
  theme_minimal()

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




