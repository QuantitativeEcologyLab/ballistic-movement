# This script generates the functions necessary for carrying out the 
# prey simulation study aimed at exploring the evolution of ballistic motion


#Written by Michael Noonan and Lynndsay Terpsma

#Last updated: December 16, 2025

#.........................................................................
# Package import

library(ctmm) #for generating movement models
library(spatstat) #for creating point process landscapes
library(RANN) # for calculating nearest neighbors

#.........................................................................
# Calculate the euclidean distance between two points----
#.........................................................................

#used to help calculate encounters
SLD <- function(x_1, y_1, x_2, y_2){
  sqrt((x_1 - x_2)^2 + (y_1 - y_2)^2)
}

#.........................................................................
# re-parameterize rgamma() as a function of mean and variance----
#.........................................................................

#used to add variance in movement parameters and food raster
rgamma2 <- function(mu, sigma2, N = n()) {
  # mean = k * theta
  # sigma^2 = k * theta^2
  rgamma(n = N,
         shape = mu^2 / sigma2, # (k * theta)^2 / (k * theta^2)
         scale = sigma2 / mu) # (k * theta^2) / (k * theta)
}



#.........................................................................
# Generate prey movement model based on prey mass (in g)----
#.........................................................................

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

#.........................................................................
# Generate prey var[position] based on mass (in g)----
#.........................................................................

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

#.........................................................................
# Generate prey E[tau_p] based on mass (in g)----
#.........................................................................

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

#.........................................................................
# Generate prey E[tau_v] based on mass (in g)----
#.........................................................................

# Model comes from Noonan et al. 2020  https://doi.org/10.1111/cobi.13495

prey.tau_v <- function(mass, variance = FALSE) {
  #Calculate
  tau_v <- -1.365200 + 0.787177*log10(mass)
  #Back transform
  tau_v <- 10^(tau_v)
  #Add variance if desired
  if(variance == TRUE){
    sigma2 <- tau_v *10
    tau_v <- rgamma2(tau_v, sigma2, N = length(mass))}
  #Return
  return(tau_v)
}

#.........................................................................
# Generate point process of food patches based on mass_prey (g)----
#.........................................................................

makeHabitat <- function(mass, 
                        kappa, 
                        r, 
                        mu, 
                        tile_size = 500, 
                        cv = 0,
                        target_n = NULL,
                        cal = 1){
  sig <- prey.SIG(mass)
  EXT <- round(sqrt((-2*log10(0.001)*pi)*sig))
  win <- owin(c(-EXT, EXT), c(-EXT, EXT))
  
  win_x <- win$xrange
  win_y <- win$yrange
  
  xs <- seq(floor(win_x[1])-tile_size, ceiling(win_x[2])+tile_size, by = tile_size)
  ys <- seq(floor(win_y[1])-tile_size, ceiling(win_y[2])+tile_size, by = tile_size)
  
  pts_list <- list()
  idx <- 1
  
  for(i in seq_len(length(xs) - 1)) {
    for(j in seq_len(length(ys) - 1)) {
      tile_win <- owin(c(xs[i], xs[i+1]), c(ys[j], ys[j+1]))
      
      local_pp <- rMatClust(kappa=kappa, r=r, mu=mu, win=tile_win)
      
      pts_list[[idx]] <- local_pp
      idx <- idx+1
    }
  }
  pp_all <- do.call(superimpose, pts_list)
  pp_all[win]
  
  if(!is.null(target_n)){
    N <- npoints(pp_all)
    p <- target_n/N*1.01 #thin to slightly over desired points
    pp_all <- rthin(pp_all, p)
    
    if(npoints(pp_all) < target_n){
      stop("target n greater than simulated points")
    } 
    
    #refine to exact desired
    if(npoints(pp_all) > target_n){
      pp_all <- pp_all[sample.int(npoints(pp_all), target_n)]
    }
  }
  
  # control CoV via gamma
  if (cv > 0){
    vals <- rgamma2(cal, cv, pp_all$n)
  } else {
    vals <- rep(cal, pp_all$n)
  }
  
  marks(pp_all) <- vals
  
  return(pp_all)
}

#.........................................................................
# food raster function utilizing patches per 95% HR area 
#.........................................................................

createFoodRaster <- function(mass, 
                             k = 240000, # number of patches in the 95% HR area
                             pred = FALSE, 
                             calories = 1,
                             var = 0) {
  
  #var[position]
  if(pred){SIG <- pred.SIG(mass)} else{
    SIG <- prey.SIG(mass)}
  
  #range of raster based on 99.9% HR area
  EXT <- round(sqrt((-2*log(0.0001)*pi)* SIG))
  
  #95% HR radius
  HR <- round(sqrt((-2*log(0.05)*pi)*SIG))

  #95% HR area
  HR_area <- pi * HR^2
  #area of each patch based on set number of patches in 95% HR
  patch_area <- HR_area / k #where k is the number of patches in the 95% HR
  #back calculate the width of each patch
  width <- sqrt(patch_area)
  
  #assign number of cells based on EXT and width
  N <- ceiling(2 * EXT / width)
  
  #create raster with terra
  food_raster <- rast(ncol = N, nrow = N,
                      xmin = -EXT, xmax = EXT,
                      ymin = -EXT, ymax = EXT)
  
  #assign caloric values to raster
  #heterogeneous landscapes created with a gamma distribution based on the mean and variance
  if (var == 0) {
    values(food_raster) <- calories 
  } else {
    #create variance as a function of the number of calories
    #this helps keep is more consistent across the calorie spectrum (and thus the mass spectrum)
    sigma2 <- calories * var # creates sigma value, increase var to increase the variance
    values(food_raster) <- rgamma2(mu = calories, sigma2 = sigma2, N = ncell(food_raster))
  }
  
  #return calorie raster
  return(food_raster)
}

#.........................................................................
# Count the number of patches visited (assumes immediate renewal)----
#.........................................................................

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

#.........................................................................
# grazing function optimized for point process habitat ----
#.........................................................................

grazing_point <- function(mass, track, habitat){
  # SIG <- prey.SIG(mass)
  # HR <- round(sqrt((-2*log(0.05)*pi)*SIG))
  # #95% HR area
  # HR_area <- pi * HR^2
  # #area of each patch based on set number of patches in 95% HR
  # patch_area <- HR_area / 500 #where k is the number of patches in the 95% HR
  # #back calculate the width of each patch
  # pr <- sqrt(patch_area)
  
  pr <- sqrt(prey.SIG(mass))*0.05
  
  hab_df <- data.frame(habitat$x, habitat$y)
  hab_mat <- as.matrix(hab_df)
  
  track_df <- data.frame(x = track$x, y = track$y)
  track_mat <- as.matrix(track_df)
  
  nn <- nn2(
    data = hab_mat,
    query = track_mat,
    radius = pr,
    searchtype = "radius"
  )
  
  nearest_id <- nn$nn.idx[,1]
  
  within_range <- nearest_id != 0 
  
  marks_vec <- marks(habitat)
  
  IDs <- ifelse(within_range, nearest_id, NA)
  
  PATCHES <- sum(diff(IDs) != 0, na.rm = TRUE)
  
  entry_ids <- c(1, which(diff(IDs) != 0) + 1)
  entry_IDs <- IDs[entry_ids]
  
  patch_values <- marks_vec[entry_IDs]

  attr(IDs, "patch_values") <- patch_values
  attr(IDs, "patches") <- PATCHES
  
  return(IDs)
}

#.........................................................................
# extract speed----
#.........................................................................

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

#.........................................................................
# define "lifespan" and sampling interval----
#.........................................................................

#sampling function with lifespan scaled to body mass

sampling <- function(mass, x = 40.5) {
  
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

#.........................................................................
# net calories from grazing----
#.........................................................................

prey.cals.net <- function(IDs, mass, speed, t){
  
  time_total <- attr(t, "time_total")
  
  #extract calorie values from which the movement track overlaps
  patch_values <- attr(IDs, "patch_values")
  cal_gross <- sum(patch_values, na.rm = TRUE)
  
  # #metabolic rate (kj/day) from Nagy 1987 https://doi.org/10.2307/1942620
  BMR <- 0.774 + 0.727 * log10(mass)
  #back transform
  BMR <- 10^BMR
  #convert to kcal/s
  BMR <- (BMR * 0.239005736) / 86400

  #calculate total BMR cost over sample period
  cost_BMR <- BMR * time_total
  
  #calculate movement cost (watts/kg) from Taylor et al. 1982 https://doi.org/10.1242/jeb.97.1.1
  E <- 10.7 * (mass / 1000)^(-0.316) * speed + 6.03 * (mass / 1000)^(-0.303)  #convert to kJ/s
  E <- (E * mass/1000)/1000
  #convert to kcal/s
  E <- E * 0.239005736
  
  #calculate total movement costs
  #kcal/s to kcal 
  cost_move <- E * time_total
  
  #calculate total energetic costs
  cost_total <- cost_BMR + cost_move
  
  #assign net calories
  cal_net <- cal_gross - cost_total
  
  #return cal_net and cal_max
  return(list(cal_net = cal_net, costs = cost_total))
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

#.........................................................................
# Prey fitness function----
#.........................................................................

#calculate fitness 
prey.fitness <- function(mass, 
                         cal_net,
                         costs = NULL) 
{
  #standardize mass input
  if (length(mass) == 1) mass <- rep(mass, n_prey)
  
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
  ##therefore dry mass $\approx$ 0.25 from Young et al. 2021 https://doi.org/10.1136/archdischild-2020-321112
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
  
  return(list(offspring = offspring, mass_update = mass.update))
}

