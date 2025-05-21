# This scripts generates the functions necessary for carrying out the 
# simulation study aimed at exploring the evolution of ballistic motion


#Written by Michael Noonan and Lynndsay Terpsma

#Last updated: May 20th 2025


#----------------------------------------------------------------------
# Package import

library(ctmm)
library(raster)
library(terra)

#----------------------------------------------------------------------
# Calculate the euclidean distance between two points
#----------------------------------------------------------------------

SLD <- function(x_1, y_1, x_2, y_2){
  sqrt((x_1 - x_2)^2 + (y_1 - y_2)^2)
}

#----------------------------------------------------------------------
# Generate prey movement model based on prey mass (in g)
#----------------------------------------------------------------------

# Model comes from Noonan et al. 2020  https://doi.org/10.1111/cobi.13495

prey.mod <- function(mass, mu = c(0,0), variance = FALSE){
  #Calculate
  HR <- 0.5078955 + 1.372162*log10(mass)
  #Back transform
  HR <- 10^(HR)
  #Add variance if desired
  if(variance == TRUE){HR <- rchisq(n = length(mass), df = HR)}
  
  #Convert from 95% HR to var[position]
  SIG <- HR/(-2*log(0.05)*pi)
  
  
  #Calculate tau_p based on correlation between tau_p and 95% HR
  tau_p <- 1.115028 + 0.576379*log10(HR) + rnorm(n = length(mass),
                                                 mean = 0,
                                                 sd = sqrt(0.3945408))
  #Back transform
  tau_p <- 10^(tau_p)
  
  #Calculate tau_v based on correlation between tau_v and 95% HR
  tau_v <- 0.7840590 + 0.2396508*log10(HR) + rnorm(n = length(mass),
                                                   mean = 0,
                                                   sd = sqrt(1.136595))
  #Back transform
  tau_v <- 10^(tau_v)
  
  mod <- ctmm(tau = c(tau_p,tau_v),
              mu = mu,
              sigma = SIG)
  
  #Return
  return(mod)
}

#----------------------------------------------------------------------
# Generate prey var[position] based on mass (in g)
#----------------------------------------------------------------------

# Model comes from Noonan et al. 2020  https://doi.org/10.1111/cobi.13495

prey.SIG <- function(mass, variance = FALSE) {
  #Calculate
  HR <- 0.5078955 + 1.372162*log10(mass)
  #Back transform
  HR <- 10^(HR)
  #Convert from 95% HR to var[position]
  SIG <- HR/(-2*log(0.05)*pi)
  #Add variance if desired
  if(variance == TRUE){SIG <- rchisq(n = length(mass), df = SIG)}
  #Return
  return(SIG)
}

#----------------------------------------------------------------------
# Generate prey E[tau_p] based on mass (in g)
#----------------------------------------------------------------------

# Model comes from Noonan et al. 2020  https://doi.org/10.1111/cobi.13495

prey.tau_p <- function(mass, variance = FALSE) {
  #Calculate
  tau_p <- 1.2994292 + 0.8129125*log10(mass)
  #Back transform
  tau_p <- 10^(tau_p)
  #Add variance if desired
  if(variance == TRUE){tau_p <- rchisq(n = length(mass), df = tau_p)}
  #Return
  return(tau_p)
}

#----------------------------------------------------------------------
# Generate prey E[tau_v] based on mass (in g)
#----------------------------------------------------------------------

# Model comes from Noonan et al. 2020  https://doi.org/10.1111/cobi.13495

prey.tau_v <- function(mass, variance = FALSE) {
  #Calculate
  tau_v <- -1.365200 + 0.787177*log10(mass)
  #Back transform
  tau_v <- 10^(tau_v)
  #Add variance if desired
  if(variance == TRUE){tau_v <- rchisq(n = length(mass), df = tau_v)}
  #Return
  return(tau_v)
}

#----------------------------------------------------------------------
# Generate predator var[position] based on mass (in g)
#----------------------------------------------------------------------

# Model comes from Noonan et al. 2020  https://doi.org/10.1111/cobi.13495

pred.SIG <- function(mass, variance = FALSE) {
  #Calculate
  HR <- 1.089972 + 1.478050*log10(mass)
  #Back transform
  HR <- 10^(HR)
  #Convert from 95% HR to var[position]
  SIG <- HR/(-2*log(0.05)*pi)
  #Add variance if desired
  if(variance == TRUE){SIG <- rchisq(n = length(mass), df = SIG)}
  #Return
  return(SIG)
}


#----------------------------------------------------------------------
# Generate predator E[tau_p] based on mass (in g)
#----------------------------------------------------------------------

# Model comes from Noonan et al. 2020  https://doi.org/10.1111/cobi.13495

pred.tau_p <- function(mass, variance = FALSE) {
  #Calculate
  tau_p <- 1.612761 + 0.766461*log10(mass)
  #Back transform
  tau_p <- 10^(tau_p)
  #Add variance if desired
  if(variance == TRUE){tau_p <- rchisq(n = length(mass), df = tau_p)}
  #Return
  return(tau_p)
}

#----------------------------------------------------------------------
# Generate predator E[tau_v] based on mass (in g)
#----------------------------------------------------------------------

# Model comes from Noonan et al. 2020  https://doi.org/10.1111/cobi.13495

pred.tau_v <- function(mass, variance = FALSE) {
  #Calculate
  tau_v <- -0.1005302 + 0.7403169*log10(mass)
  #Back transform
  tau_v <- 10^(tau_v)
  #Add variance if desired
  if(variance == TRUE){tau_v <- rchisq(n = length(mass), df = tau_v)}
  #Return
  return(tau_v)
}
#----------------------------------------------------------------------
# Generate E[mass_prey] based on mass_pred (in g)
#----------------------------------------------------------------------

# Model comes from Tucker & Rogers 2014  https://doi.org/10.1371/journal.pone.0106402

prey.mass <- function(mass, variance = FALSE) {
  #Convert to kg
  mass <- mass * (1 %#% "g")
  #Calculate
  prey_mass <- -0.87 + 1.26*log10(mass)
  #Back transform
  prey_mass <- 10^(prey_mass)
  #Convert to g
  prey_mass <-prey_mass / (1 %#% "g")
  #Add variance if desired
  if(variance == TRUE){tau_v <- rchisq(n = length(mass), df = prey_mass)}
  #Return
  return(prey_mass)
}

#----------------------------------------------------------------------
# Generate raster of food patches based on mass_prey (g)
#----------------------------------------------------------------------

createFoodRaster <- function(mass, width = 20, pred = FALSE, calories = 10,
                    type = c("uniform", "random")) {
  
  type <- match.arg(type)
  
  # var[position]
  if(pred){SIG <- pred.SIG(mass)} else{
    SIG <- prey.SIG(mass)}
  
  # Range of raster based on 99.9% HR area
  EXT <- round(sqrt((-2*log(0.0001)*pi)* SIG))
  
  #number of patches based on fixed patch width
  #N <- EXT/n
  N <- ceiling(2*EXT/width)
  
  #Build the raster
  FOOD <- rast(ncol = N, nrow = N, 
               xmin = -EXT, xmax = EXT, 
               ymin = -EXT, ymax = EXT)
  
  # assign values based on type of habitat
  if (type == "uniform") {
    cell_values <- 1
  } else if (type == "random")
  {cell_values <- runif(ncell(FOOD), min = 0.1, max = 1.5)
  }  
  
  cell_values <- cell_values * calories
  
  values(FOOD) <- cell_values
  
  #Return the raster of food patches
  return(FOOD)
}

#----------------------------------------------------------------------
# Count the number of patches visited (assumes immediate renewal)
#----------------------------------------------------------------------

grazing <- function(track, habitat, metric = "ids") {
  
  #convert track to data frame
  coords <- data.frame(x = track$x, y = track$y)
  
  #Patch identities
  IDs <- cellFromXY(habitat, coords)
  
  # Count the number of times it moved to a new food patch
  PATCHES <- sum(diff(IDs) != 0)
  
  #Mean time between patches 
  TIME <- mean(rle(c(FALSE, diff(IDs) != 0))$lengths)
  
  if(metric == "patches"){return(PATCHES)}
  if(metric == "time"){return(TIME)}
  if(metric == "ids") {return(IDs)}
    stop("Invalid metric. Use 'patches', 'time' or 'ids'.")
}

#----------------------------------------------------------------------
# calculate speed
#----------------------------------------------------------------------

speed_val <- function(models){
  # Extract movement speeds from the models
  model_summary <- summary(models, units = FALSE)
  
  # Ensure model_summary has $CI before accessing
  if(!is.null(model_summary$CI) && nrow(model_summary$CI) == 5){
    SPEED <- model_summary$CI[4, "est"]
  } else {
    SPEED <- Inf
  }
  return(SPEED)
}

#----------------------------------------------------------------------
# net benefits attempt
#----------------------------------------------------------------------

cals_net <- function(IDs, habitat, mass, models, speed, interval){

  #gain in calories
  patch_values <- values(habitat)[IDs]
  gain_total <- sum(patch_values, na.rm = TRUE)
  
  # Basal metabolic rate (in kj/day) from Nagy 1987 https://doi.org/10.2307/1942620 
  BMR <- 0.774 + 0.727*log10(mass)
  #Back transform 
  BMR <- 10^BMR
  #convert to calories/sec
  BMR <- (BMR * 239.005736) / 86400
  
  #lifespan calculations
  lifespan <- (4.88*mass^0.153) * 31536000 # years to seconds
  time_total <- lifespan * 0.001 # 1/1000 of a lifespan
  
  BMR_cost <- BMR * time_total
  
  # Metabolic cost of movement in watts/kg from Taylor et al. 1982 https://doi.org/10.1242/jeb.97.1.1 
  E = 10.7*(mass/1000)^(-0.316)*speed + 6.03*(mass/1000)^(-0.303)
  #Convert to kJ/s
  E <- (E * (mass/1000))/1000
  #Convert to calories/s
  E <- E * 239.005736
  
  num_movements <- sum(diff(IDs) != 0)
  cost_move <- num_movements*E*interval 
  
  cost_total <- BMR_cost + cost_move
  
  cal_net <- gain_total - cost_total
  
  return(cal_net)
}

#----------------------------------------------------------------------
# define "lifespan" and sampling interval
#----------------------------------------------------------------------

#sampling function with lifespan scaled to body mass

sampling <- function(mass, metric = "t") {
  
  #calculate lifespan in seconds
  #de Magalhaes et al (2008) 
  lifespan <- (4.88*mass^0.153) * 31536000 # years to seconds
  lifespan_int <- lifespan * 0.001 # 1/1000 of a lifespan
  
  #sampling interval (tau_v) in seconds
  interval <- max(1, round(prey.tau_v(mass)))
  
  #lifespan and sampling interval for simulations
  t <- seq(0,
           lifespan_int,
           interval)
  
  #return vector of sampling times
  if(metric == "t"){return(t)}
  if(metric == "lifespan"){return(lifespan)}
  if(metric == "interval"){return(interval)}
  stop("Invalid metric. Use 't' or 'offspring' or 'lifespan'.")
}


#----------------------------------------------------------------------
# Prey fitness function 
#----------------------------------------------------------------------

#calculate fitness 
prey.fitness <- function(mass, 
                         cal_net,
                         costs = NULL) 
{
  
  #DEBkiss model: Jager T (2024). DEBkiss. A simple framework for animal energy budgets. 
  #Version 3.1. Leanpub: https://leanpub.com/debkiss_book.
  
  # Standardize mass input
  if (length(mass) == 1) mass <- rep(mass, n_prey)
  
  #update weight
  #let weight_gain (g) = calorie surplus (kcal) / energy cost of tissue (kcal/g)
  #estimate value for cost from Waterlow et al. (1981)
  weight.gain <- cal_net / 20
  mass.update <- mass+weight.gain
  
  #reproductive buffer, can use litter mass as a proxy, therefore
  #m_litter = 0.637 * mass^0.778 from Huijsmans et al. (2024)
  W_R <- 0.637*(mass.update^0.778) 

  #birth weight via allometric scaling in mammals (Blueweiss et al. 1978)
  #wet weight $\approx$ 0.75 total weight, therefore dry mass $\approx$ 0.25 (Fusch et al. 1999)
  W_B0 <- 0.25*(0.097*mass.update^(0.92))
  
  #total offspring scaled to functional response
  offspring <- floor(W_R/W_B0) 
 
  #set offspring to 0 is cal_net <= 0
  offspring[cal_net <= 0] <- 0
  
  #If predator encounters are being considered,
  #individuals that encountered a predator are killed and don't reproduce.
  if(!is.null(costs)){
    offspring[costs] <- 0
  }
  
  #clamp minimum offspring to 0
  offspring <- ctmm:::clamp(offspring, min = 0, max = Inf) #Clamp the minimum to 0
  
  return(list(offspring = offspring, mass_update = mass.update))
}

#----------------------------------------------------------------------
# Identify Encounter Events
#----------------------------------------------------------------------

encounter <- function(prey.tracks, pred.tracks, range = 50){
  distances <- list()
  encounters <- vector()
  for(i in 1:length(prey.tracks)){
    #Pairwise separation distances over time
    distances[[i]] <- SLD(pred.tracks[[1]]$x,pred.tracks[[1]]$y,
                          prey.tracks[[i]]$x, prey.tracks[[i]]$y)
    
    #Did it encounter a predator
    encounters[i] <- any(distances[[i]]<range)
  }
  return(encounters)
}

#----------------------------------------------------------------------
# Predator fitness function
#----------------------------------------------------------------------

pred.fitness <- function(encounters, mass, costs = NULL, models, time = t, calories = 10, constant = 1){
  
  # Extract movement speeds from the models
  SPEED <- vector()
  for(i in 1:length(models)){SPEED[i] <- if(nrow(summary(models[[i]], units = FALSE)$CI)==4){summary(models[[i]], units = FALSE)$CI[4,2]} else{Inf}}
  
  # Basal metabolic rate (in kj/day) from Nagy 1987 https://doi.org/10.2307/1942620 
  BMR <- 0.774 + 0.727*log10(mass)
  
  #Back transform 
  BMR <- 10^BMR
  
  # total lifespan in days (based on number of range crossings)
  #lifespan <- round(pred.tau_p(mass)*crossings) /60/60/24
  lifespan <- tail(t, n=1) /60/60/24
  
  # Metabolic cost of movement in watts/kg from Taylor et al. 1982 https://doi.org/10.1242/jeb.97.1.1 
  E = 10.7*(mass/1000)^(-0.316)*SPEED + 6.03*(mass/1000)^(-0.303)
  
  #Convert to kJ/s
  E <- (E * (mass/1000))/1000
  
  # Maximum running speed in km/hr from Hirt et al. 2017 https://doi.org/10.1038/s41559-017-0241-4
  v_max <- 25.5 * (mass/1000)^(0.26) * (1 - exp(-22*(mass/1000)^(-0.66)))
  
  #Convert to m/s
  v_max <- v_max/3.6
  
  #Total energetic cost in kj as a function of BMR and movement speed
  COST <- BMR * lifespan + E*tail(t, n=1)*constant
  
  #Energy intake in kj based on Gorecki 1965 http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.915.5227&rep=rep1&type=pdf
  mass_prey <- prey.mass(mass)
  intake <- 1.5 * mass_prey * 4.184 * sum(encounters)
  
  # Excess energy
  excess <- intake - COST
  excess[is.infinite(excess)] <- NA
  
  # Define number of prey offspring based on their excess energy and metabolic rate
  offspring <- floor(excess/BMR)
  offspring[is.na(offspring)] <- 0
  offspring <- ctmm:::clamp(offspring, min = 0, max = Inf) #Clamp the minimum to 0
  
  return(offspring)
}

#----------------------------------------------------------------------
# 
#----------------------------------------------------------------------
