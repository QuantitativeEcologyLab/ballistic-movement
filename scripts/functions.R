# This scripts generates the functions necessary for carrying out the 
# simulation study aimed at exploring the evolution of ballistic motion


#Written by Michael Noonan and Lynndsay Terpsma

#Last updated: June 4th 2025


#----------------------------------------------------------------------
# Package import

library(ctmm) #for generating movement models
library(raster) #for creating the food raster
library(terra) #for creating the food raster

#----------------------------------------------------------------------
# create figure theme----
#----------------------------------------------------------------------

theme.qel <- function(legend = TRUE){
  theme <- theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_text(size=6, family = "sans", face = "bold"),
          axis.title.x = element_text(size=6, family = "sans", face = "bold"),
          axis.text.y = element_text(size=4, family = "sans"),
          axis.text.x  = element_text(size=4, family = "sans"),
          plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
          plot.background = element_rect(fill = "transparent", color = NA),
          plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
  if(legend){
    theme <- theme +
      theme(
        legend.position = "right",
        legend.text = element_text(size = 4, family = "sans"),
        legend.title = element_text(size = 5, family = "sans", face = "bold"),
        legend.key.size = unit(0.2, "cm"),
        legend.spacing.y = unit(0.1, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent"))
  } else {
     theme <- theme + 
       theme(legend.position = "none",
             panel.grid = element_blank())
   }

  return(theme)
  
}

#----------------------------------------------------------------------
# Calculate the euclidean distance between two points----
#----------------------------------------------------------------------

#used to help calculate encounters
SLD <- function(x_1, y_1, x_2, y_2){
  sqrt((x_1 - x_2)^2 + (y_1 - y_2)^2)
}

#----------------------------------------------------------------------
# re-parameterize rgamma() as a function of mean and variance----
#----------------------------------------------------------------------

#used to add variance in movement parameters and food raster
rgamma2 <- function(mu, sigma2, N = n()) {
  # mean = k * theta
  # sigma^2 = k * theta^2
  rgamma(n = N,
         shape = mu^2 / sigma2, # (k * theta)^2 / (k * theta^2)
         scale = sigma2 / mu) # (k * theta^2) / (k * theta)
}

#----------------------------------------------------------------------
# Generate prey movement model based on prey mass (in g)----
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
# Generate prey var[position] based on mass (in g)----
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
  if(variance == TRUE){
    sigma2 <-SIG * 10
    SIG <- rgamma2(SIG, sigma2, N = length(mass))
  }
  #Return
  return(SIG)
}

#----------------------------------------------------------------------
# Generate prey E[tau_p] based on mass (in g)----
#----------------------------------------------------------------------

# Model comes from Noonan et al. 2020  https://doi.org/10.1111/cobi.13495

prey.tau_p <- function(mass, variance = FALSE) {
  #Calculate
  tau_p <- 1.2994292 + 0.8129125*log10(mass)
  #Back transform
  tau_p <- 10^(tau_p)
  #Add variance if desired
  if(variance == TRUE){
    sigma2 <- tau_p * 10
    tau_p <- rgamma2(tau_p, sigma2, N = length(mass))}
  #Return
  return(tau_p)
}

#----------------------------------------------------------------------
# Generate prey E[tau_v] based on mass (in g)----
#----------------------------------------------------------------------

# Model comes from Noonan et al. 2020  https://doi.org/10.1111/cobi.13495

prey.tau_v <- function(mass, variance = FALSE) {
  #Calculate
  tau_v <- -1.365200 + 0.787177*log10(mass)
  #Back transform
  tau_v <- 10^(tau_v)
  #Add variance if desired
  if(variance == TRUE){
    sigma2 <- tau_v * 10
    tau_v <- rgamma2(tau_v, sigma2, N = length(mass))}
  #Return
  return(tau_v)
}

#----------------------------------------------------------------------
# Generate predator var[position] based on mass (in g)----
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
  if(variance == TRUE){
    sigma2 <-SIG * 10
    SIG <- rgamma2(SIG, sigma2, N = length(mass))}
  #Return
  return(SIG)
}

#----------------------------------------------------------------------
# Generate predator E[tau_p] based on mass (in g)----
#----------------------------------------------------------------------

# Model comes from Noonan et al. 2020  https://doi.org/10.1111/cobi.13495

pred.tau_p <- function(mass, variance = FALSE) {
  #Calculate
  tau_p <- 1.612761 + 0.766461*log10(mass)
  #Back transform
  tau_p <- 10^(tau_p)
  #Add variance if desired
  if(variance == TRUE){
    sigma2 <- tau_p * 10
    tau_p <- rgamma2(tau_p, sigma2, N = length(mass))
  }
  #Return
  return(tau_p)
}

#----------------------------------------------------------------------
# Generate predator E[tau_v] based on mass (in g)----
#----------------------------------------------------------------------

# Model comes from Noonan et al. 2020  https://doi.org/10.1111/cobi.13495

pred.tau_v <- function(mass, variance = FALSE) {
  #Calculate
  tau_v <- -0.1005302 + 0.7403169*log10(mass)
  #Back transform
  tau_v <- 10^(tau_v)
  #Add variance if desired
  if(variance == TRUE){
    sigma2 <- tau_v * 10
    tau_v <- rgamma2(tau_v, sigma2, N = length(mass))}
  #Return
  return(tau_v)
}
#----------------------------------------------------------------------
# Generate E[mass_prey] based on mass_pred (in g)----
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
# Generate raster of food patches based on mass_prey (g)----
#----------------------------------------------------------------------

# food raster function utilizing patches per 95% HR area (new)

createFoodRaster <- function(mass, k, pred = FALSE, 
                             calories = 0.015) {
  
  #var[position]
  if(pred){SIG <- pred.SIG(mass)} else{
    SIG <- prey.SIG(mass)}
  
  #range of raster based on 99.9% HR area
  EXT <- round(sqrt((-2*log(0.0001)*pi)* SIG))
  
  # 95% HR radius
  HR <- round(sqrt((-2*log(0.05)*pi)*SIG))

  # 95% HR area
  HR_area <- pi * HR^2
  # area of each patch based on set number of patches in 95% HR
  patch_area <- HR_area / k #where k is the number of patches in the 95% HR
  # back calculate the width of each patch
  width <- sqrt(patch_area)
  
  # assign number of cells based on EXT and width
  N <- ceiling(2 * EXT / width)
  
  #create raster with terra
  food_raster <- rast(ncol = N, nrow = N,
                      xmin = -EXT, xmax = EXT,
                      ymin = -EXT, ymax = EXT)
  
  #assign caloric values to cells
  values(food_raster) <- calories # calories is defined as calories per patch
  
  #return calorie raster
  return(food_raster)
}

#..............................................................................
# old food raster function, using patch width scaled to prey.SIG
# use for seed2STEM heterogeneity investigation

makefood <- function(mass, width, pred = FALSE, 
                     calories = 0.015, # calories per unit area
                     heterogeneity = FALSE) {
  # width = round(sqrt(prey.SIG(mass_prey))/10
  
  #var[position]
  if(pred){SIG <- pred.SIG(mass)} else{
    SIG <- prey.SIG(mass)}
  
  #range of raster based on 99.9% HR area
  EXT <- round(sqrt((-2*log(0.0001)*pi)* SIG))
  
  #calculate the number of cells based on the EXT and the cell width
  N <- EXT / width
  
  #calculate patch area from width
  patch_area <- width^2
  
  #create raster with terra
  food_raster <- rast(ncol = N, nrow = N,
                      xmin = -EXT, xmax = EXT,
                      ymin = -EXT, ymax = EXT)
  
  cal_per_patch <- calories * patch_area # convert cal_per_m2 to cal_per_patch
  var <- 10
  sigma2 <- cal_per_patch * var # creates sigma value, increase var to increase the variance
  
  #assign caloric values to raster
  if (heterogeneity) {
    values(food_raster) <- rgamma2(mu = cal_per_patch, sigma2 = sigma2, N = ncell(food_raster))
  } else {
    values(food_raster) <- cal_per_patch
  }
  #return calorie raster
  return(food_raster)
}

#----------------------------------------------------------------------
# Count the number of patches visited (assumes immediate renewal)----
#----------------------------------------------------------------------

grazing <- function(track, habitat) {
  
  #convert track to data frame
  coords <- data.frame(x = track$x, y = track$y)
  
  #patch identities
  IDs <- cellFromXY(habitat, coords)
  
  #count the number of times it moved to a new food patch (total movements)
  PATCHES <- sum(diff(IDs) != 0)
  
  #mean time between patches 
  TIME <- mean(rle(c(FALSE, diff(IDs) != 0))$lengths)
  
  #find indices (cell identity) for when it moved to a new food patch
  entry_ids <- c(1, which(diff(IDs) != 0) + 1)
  
  #get IDs for those patches
  entry_IDs <- IDs[entry_ids]
  
  #get values from the raster
  patch_values <- habitat[entry_IDs]
  
  #attributes
  attr(IDs, "patch_values") <- patch_values
  attr(IDs, "patches") <- PATCHES
  attr(IDs, "time") <- TIME
 
  return(IDs) 
}

#----------------------------------------------------------------------
# extract speed----
#----------------------------------------------------------------------

get.speed <- function(models){
  #extract movement speeds from the models
  model_summary <- summary(models, units = FALSE)
  
  #ensure model_summary has $CI before accessing
  if(!is.null(model_summary$CI) && nrow(model_summary$CI) == 5){
    SPEED <- model_summary$CI[4, 2]
  } else {
    SPEED <- Inf
  }
  
  #return speed
  return(SPEED)
}

#----------------------------------------------------------------------
# define "lifespan" and sampling interval----
#----------------------------------------------------------------------

#sampling function with lifespan scaled to body mass

sampling <- function(mass, x = 1) {
  
  #calculate lifespan in seconds from de Magalhaes et al (2008) https://doi.org/10.1093/gerona/62.2.149
  lifespan <- (4.88*mass^0.153) * 31536000 # years to seconds
  time_total <- lifespan * 0.001 # 1/1000 of a lifespan
  
  #sampling interval (tau_v) in seconds
  #increasing x decreases interval, making sampling more frequent
  interval <- max(1, round(prey.tau_v(mass))) / x
  
  #lifespan and sampling interval for simulations
  t <- seq(0,
           time_total,
           interval)
  
  # assign attributes
  attr(t, "interval") <- interval
  attr(t, "lifespan") <- lifespan
  attr(t, "time_total") <- time_total
  
  return(t)
}

#----------------------------------------------------------------------
# net calories from grazing----
#----------------------------------------------------------------------

prey.cals.net <- function(IDs, mass, speed, t){
  
  time_total <- attr(t, "time_total")
  
  #extract calorie values from which the movement track overlaps
  patch_values <- attr(IDs, "patch_values")
  cal_gross <- sum(patch_values, na.rm = TRUE)
  
  #metabolic rate (kj/day) from Nagy 1987 https://doi.org/10.2307/1942620
  # BMR <- 0.774 + 0.727 * log10(mass)
  # #back transform
  # BMR <- 10^BMR
  # #convert to cal/s
  # BMR <- (BMR * 239.005736) / 86400
  
  #calculate total BMR cost over sample period 
  # BMR_cost <- BMR * time_total
  
  #calculate movement cost (watts/kg) from Taylor et al. 1982 https://doi.org/10.1242/jeb.97.1.1
  E <- 10.7 * (mass / 1000)^(-0.316) * speed + 6.03 * (mass / 1000)^(-0.303)  #convert to kJ/s
  E <- (E * mass/1000)/1000
  #convert to kcal/s
  E <- E * 0.239005736
  
  #calculate total movement costs
  #cal/s to cal 
  move_cost <- E * time_total
  
  #calculate total energetic costs
  # cost_total <- BMR_cost * 0.2 + move_cost
  
  #assign net calories
  cal_net <- cal_gross - move_cost
  
  #return cal_net and cal_max
  return(cal_net)
}

#.........................................................................
# test function with no costs of movement of metabolism
#.........................................................................

prey_cals_net_nocost <- function(IDs){
  
  #extract calorie values from which the movement track overlaps
  patch_values <- attr(IDs, "patch_values")
  cal_gross <- sum(patch_values, na.rm = TRUE)
  
  #assign net calories
  cal_net <- cal_gross
  
  #return cal_net and cal_max
  return(cal_net)
}

#----------------------------------------------------------------------
# Prey fitness function----
#----------------------------------------------------------------------

#calculate fitness 
prey.fitness <- function(mass, 
                         cal_net,
                         costs = NULL) 
{
  #standardize mass input
  # if (length(mass) == 1) mass <- rep(mass, n_prey)
  
  #update weight
  cal_net[cal_net < 0] <- 0 #prevent negative
  growth_cal <- cal_net*0.8 #allocation to soma
  repro_cal <- cal_net*0.2 #allocation to reproduction
  
  weight.gain <- growth_cal / 5
  mass.update <- mass + weight.gain
  
  #using mass allocated to reproduction to determine W_R
  W_R <- repro_cal / 5
  
  #birth weight via allometric scaling in mammals from Blueweiss et al. 1978 https://doi.org/10.1007/BF00344996
  #wet weight $\approx$ 0.75 total weight
  ##therefore dry mass $\approx$ 0.25 from Fusch et al. 1999 https://doi.org/10.1203/00006450-199910000-00018
  W_B0 <- 0.25*(0.097*mass.update^(0.92))
  
  #total offspring based on updated mass
  offspring <- floor(W_R/W_B0) 
  
  #set offspring to 0 is cal_net <= 0
  offspring[cal_net <= 0] <- 0
  
  #If predator encounters are being considered,
  #individuals that encountered a predator are killed and don't reproduce.
  if(!is.null(costs)){
    offspring[costs] <- 0
  }
  
  #clamp minimum offspring to 0
  offspring <- ctmm:::clamp(offspring, min = 0, max = Inf) #clamp the minimum to 0
  
  return(offspring)
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
# Predator calorie intake
#----------------------------------------------------------------------

pred_cals_net <- function(encounters, mass, t, speed){
  
  time_total <- attr(t, "time_total")
  
  # metabolic rate (kj/day) from Nagy 1987 https://doi.org/10.2307/1942620
  # BMR <- 0.774 + 0.727 * log10(mass)
  # #back transform
  # BMR <- 10^BMR
  # #convert to cal/s
  # BMR <- (BMR * 239.005736) / 86400
  
  #calculate movement cost (watts/kg) from Taylor et al. 1982 https://doi.org/10.1242/jeb.97.1.1
  E <- 10.7 * (mass / 1000)^(-0.316) * (speed) + 6.03 * (mass / 1000)^(-0.303)
  #convert to kJ/s
  E <- (E * mass/1000)/1000
  #convert to cal/s
  E <- E * 0.239005736
  
  # # Maximum running speed in km/hr from Hirt et al. 2017 https://doi.org/10.1038/s41559-017-0241-4
  # v_max <- 25.5 * (mass/1000)^(0.26) * (1 - exp(-22*(mass/1000)^(-0.66)))
  # 
  # #Convert to m/s
  # v_max <- v_max/3.6
  
  #Total energetic cost in kj as a function of BMR and movement speed
  move_cost <- E * time_total
  
  #Energy intake in cal based on Gorecki 1965 http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.915.5227&rep=rep1&type=pdf
  mass_prey <- prey.mass(mass)
  intake <- 1500 * mass_prey * sum(encounters)
  
  pred_cal_net <- intake - move_cost
  
  return(list(cal_net = pred_cal_net, costs = move_cost))
}

#----------------------------------------------------------------------
# Predator fitness ---- 
#----------------------------------------------------------------------

#calculate fitness 
pred.fitness <- function(mass, 
                         pred_cal_net) 
{
  #standardize mass input
  if (length(mass) == 1) mass <- rep(mass, n_pred)
  
  #update weight
  pred_cal_net[pred_cal_net < 0] <- 0 #prevent negative
  growth_cal <- pred_cal_net*0.8 #allocation to soma
  repro_cal <- pred_cal_net*0.2 #allocation to reproduction
  
  weight.gain <- growth_cal / 5
  mass.update <- mass + weight.gain
  
  #using mass allocated to reproduction to determine W_R
  W_R <- repro_cal / 5
  
  #birth weight via allometric scaling in mammals from Blueweiss et al. 1978 https://doi.org/10.1007/BF00344996
  #wet weight $\approx$ 0.75 total weight
  ##therefore dry mass $\approx$ 0.25 from Fusch et al. 1999 https://doi.org/10.1203/00006450-199910000-00018
  W_B0 <- 0.25*(0.097*mass.update^(0.92))
  
  #total offspring based on updated mass
  offspring <- floor(W_R/W_B0) 
  
  #set offspring to 0 is cal_net <= 0
  offspring[pred_cal_net <= 0] <- 0
  
  #clamp minimum offspring to 0
  offspring <- ctmm:::clamp(offspring, min = 0, max = Inf) #clamp the minimum to 0
  
  return(list(offspring = offspring, mass_update = mass.update))
}

