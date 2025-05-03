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





