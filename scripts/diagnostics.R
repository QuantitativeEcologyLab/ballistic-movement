# preamble

# Set the working directory
setwd("~/hdrive/GitHub/ballistic-movement")

# Import necessary packages
library(ggplot2)
library(dplyr)
library(gridExtra)
library(patchwork)

# Source the functions (ensure 'functions.R' is available in the working directory)
source("scripts/functions.R")

#----------------------------------------------------------------------
# make diagnostic figures----
#----------------------------------------------------------------------

#load in your data
load('prey_results/105500g_prey_res.Rda')
load('prey_results/105500g_prey_details.Rda')

# make data sets compatible with ggplot2
prey_res_df <- do.call(rbind, prey_res)
prey_details_df <- do.call(rbind, prey_details)

# relative change in lv ~ gen
PREY_LV <- prey_res_df$lv[1]
prey_res_df$rel.lv <- prey_res_df$lv/PREY_LV
prey_res_df$rel_var <- prey_res_df$var / (PREY_LV^2)
prey_res_df$rel_sd <- sqrt(prey_res_df$rel_var)

prey_summary <- prey_details_df %>%
  group_by(generation) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop")

# prey_details_df_sum <- prey_details_df %>% 
#   filter(generation >= max(generation) - 49) %>% 
#   ungroup()

# rel lv ~ generation
rel.lv.gen <- ggplot(prey_res_df, aes(x = generation, y = rel.lv)) +
  geom_line(color = "deeppink4", linewidth = 1) +
  geom_hline(yintercept = 1, color = 'grey30', linetype = "dashed") +
  geom_ribbon(aes(ymin = rel.lv - rel_sd,
                  ymax = rel.lv + rel_sd),
              fill = "deeppink4", alpha = 0.3) +
  labs(y = "relative change in lv",x = "generation") +
  theme_minimal()
# print(rel.lv.gen)

# tau_v ~ generation
tauv.gen <- ggplot(prey_summary, aes(x = generation, y = tau_v)) +
  geom_line(col = "deeppink4", linewidth = 1)+
  labs(x = "generation", y = "tau_v") +
  theme_minimal()
# print(tauv.gen)

# tau_p ~ generation
taup.gen <- ggplot(prey_summary, aes(x = generation, y = tau_p)) +
  geom_line(col = "deeppink4", linewidth = 1)+
  labs(x = "generation", y = "tau_p") +
  theme_minimal()
# print(taup.gen)

# sig ~ generation
sig.gen <- ggplot(prey_summary, aes(x = generation, y = sig)) +
  geom_line(col = "deeppink4", linewidth = 1) +
  labs(x = "generation", y = "sig") +
  theme_minimal()
# print(sig.gen)

# lv ~ generation
lv.gen <- ggplot(prey_summary, aes(x = generation, y = lv)) +
  geom_line(col = "deeppink4", linewidth = 1) +
  labs(x = "generation", y = "lv") +
  theme_minimal()
# print(lv.gen)

# patches ~ generation
patches.gen <- ggplot(prey_summary, aes(x = generation, y = patches)) +
  geom_line(col = "deeppink4", linewidth = 1) +
  labs(x = "generation", y = "patches visited") +
  theme_minimal()
# print(patches.gen)

# cost ~ generation
cost.gen <- ggplot(prey_summary, aes(x = generation, y = costs)) +
  geom_line(col = "deeppink4", linewidth = 1) +
  labs(x = "generation", y = "metabolic costs (cal)") +
  theme_minimal()
# print(cost.gen)

# cal_net ~ generation
cal.gen <- ggplot(prey_summary, aes(x = generation, y = cal_net)) +
  geom_line(col = "deeppink4", linewidth = 1) +
  labs(x = "generation", y = "net calories") +
  theme_minimal()
# print(cal.gen)

# speed ~ generation
speed.gen <- ggplot(prey_summary, aes(x = generation, y = speed)) +
  geom_line(col = "deeppink4", linewidth = 1) +
  labs(x = "generation", y = "speed (m/s)") +
  theme_minimal()
# print(speed.gen)

# offspring ~ generation
offspring.gen <- ggplot(prey_summary, aes(x = generation, y = offspring)) +
  geom_line(col = "deeppink4", linewidth = 1) +
  labs(x = "generation", y = "offspring") +
  theme_minimal()
# print(offspring.gen)

# mass ~ gen
mass.gen <- ggplot(prey_summary, aes(x = generation, y = mass_update)) +
  geom_line(col = "deeppink4", linewidth = 1) +
  labs(x = "generation", y = "mass_updated (g)") +
  theme_minimal()
# print(mass.gen)

# BMR ~ generation
# BMR.gen <- ggplot(prey_summary, aes(x = generation, y = BMR)) +
#   stat_summary(fun = mean, geom = "line", col = "deeppink4", linewidth = 1) +
#   labs(x = "generation",y = "BMR (cal)") +
#   theme_minimal()
# print(BMR.gen)

# Movement cost ~ generation
# move.gen <- ggplot(prey_summary, aes(x = generation, y = Move)) +
#   stat_summary(fun = mean, geom = "line", col = "deeppink4", linewidth = 1) +
#   labs(x = "generation", y = "movement costs (cal)") +
#   theme_minimal()
# print(move.gen)

# number of patches ~ lv
patches.lv <- ggplot(prey_details_df, aes(x = lv, y = patches)) +
  geom_point(col = "deeppink4", alpha = 0.3) +
  labs(x = "lv", y = "patches visited") +
  theme_minimal()
# print(patches.lv)

# number of patches ~ speed
patches.speed <- ggplot(prey_details_df, aes(x = speed, y = patches)) +
  geom_point(col = "deeppink4", alpha = 0.3) +
  labs(x = "speed (m/s)", y = "patches visited") +
  theme_minimal()

# patches ~ offspring
patches.off <- ggplot(prey_details_df, aes(x = patches, y = offspring)) +
  geom_point(col = "deeppink4", alpha = 0.3) +
  labs(x = "patch visited", y = "offspring") + 
  theme_minimal()

# cal_net ~ mass_update
cal.mass <- ggplot(prey_details_df, aes(x = mass_update, y = cal_net)) +
  geom_point(col = "deeppink4", alpha = 0.3) +
  labs(x = "mass_update (g)", y = "net calories") +
  theme_minimal()

# cal_net ~ lv
cal.lv <- ggplot(prey_details_df, aes(x = lv, y = cal_net)) +
  geom_point(col = "deeppink4", alpha = 0.3) +
  labs(x = "lv", y = "net calories") +
  theme_minimal()

# cost ~ mass_update
cost.mass <- ggplot(prey_details_df, aes(x = costs, y = mass_update)) +
  geom_point(col = "deeppink4", alpha = 0.3) +
  labs(x = "movement costs (cal)", y = "mass_update (g)") +
  theme_minimal()
# print(cost.mass)

# cost ~ speed
cost.speed <- ggplot(prey_details_df, aes(x = speed, y = costs)) +
  geom_point(col = "deeppink4", alpha = 0.3) +
  labs(x = "speed (m/s)", y = "metabolic costs (cal)") + 
  theme_minimal()

# speed ~ cal_net 
speed.cal <- ggplot(prey_details_df, aes(x = speed, y = cal_net)) +
  geom_point(col = "deeppink4", alpha = 0.3) +
  labs(x = "speed", y = "net calories") +
  theme_minimal()
# print(speed.cal)

# speed ~ lv
speed.lv <- ggplot(prey_details_df, aes(x = lv, y = speed)) +
  geom_point(col = "deeppink4", alpha = 0.3) +
  labs(x = "lv", y = "speed (m/s)") +
  theme_minimal()

# speed ~ mass_update
speed.mass <- ggplot(prey_details_df, aes(x = speed, y = mass_update)) +
  geom_point(col = "deeppink4", alpha = 0.3) +
  labs(x = "speed (m/s)", y = "mass_update (g)") +
  theme_minimal()
# print(speed.mass)

# offspring ~ speed
offspring.speed <- ggplot(prey_details_df, aes(x = speed, y = offspring)) +
  geom_point(col = "deeppink4", alpha = 0.3) +
  labs(x = "speed", y = "offspring") +
  theme_minimal()

# offspring ~ lv
offspring.lv <- ggplot(prey_details_df, aes(x = lv, y = offspring)) +
  geom_point(col = "deeppink4", alpha = 0.3) +
  labs(x = "lv", y = "offspring") +
  theme_minimal()

# offspring ~ cal_net
offspring.cal <- ggplot(prey_details_df, aes(x = cal_net, y = offspring)) +
  geom_point(col = "deeppink4", alpha = 0.3) +
  labs(x = "net calories", y = "offspring") +
  theme_minimal()

# offspring ~ mass_update
offspring.mass <- ggplot(prey_details_df, aes(x = mass_update, y = offspring)) +
  geom_point(col = "deeppink4", alpha = 0.3) +
  labs(x = "mass_update (g)", y = "offspring") +
  theme_minimal()

# Move/BMR ~ generation
# BMR.move <- ggplot(prey_details_df, aes(y = Move/BMR, x = generation)) +
#   geom_point(col = "deeppink4", alpha = 0.3) +
#   labs(y = "movement costs/BMR", x = "generation") +
#   theme_minimal()

# speed ~ tau_p
taup.speed <- ggplot(prey_details_df, aes(y = speed, x = tau_p)) +
  geom_point(col = "deeppink4", alpha = 0.3) +
  labs(y = "speed (m/s)", x = "tau_p (s)") +
  theme_minimal()

#create panel 1: figures of traits over generation time
p1 <- grid.arrange(rel.lv.gen, taup.gen, tauv.gen, sig.gen, patches.gen, speed.gen, cost.gen, offspring.gen, cal.gen, mass.gen, 
                   ncol = 2, top = "105500 grams")
ggsave(p1, file = "prey_results/figures/overview/105500g_panel1.png")

#create panel 2: figure comparing traits to one another
p2 <- grid.arrange(patches.lv, speed.lv, patches.speed, cal.lv, speed.cal, cost.speed, 
                   cal.mass, speed.mass, cost.mass, taup.speed, offspring.lv, offspring.speed,
                   offspring.mass, offspring.cal, patches.off, 
                   ncol = 3, top = "105500 grams")
ggsave(p2, file = "prey_results/figures/overview/105500g_panel2.png")

#plot raster with track
#create data frame from food raster
df_raster <- as.data.frame(FOOD, xy = TRUE)
colnames(df_raster) <- c("x", "y", "calories")

#collect movement tracks from all individuals
prey_list <- lapply(seq_along(PREY_tracks), function(i){
  df <- as.data.frame(PREY_tracks[[i]])
  colnames(df)[2:3] <- c("x", "y")
  return(df)
})

#get x and y lines to show the patches in the raster
xlines <- unique(df_raster$x)
ylines <- unique(df_raster$y)

#custome function used to draw circles with ggplot
circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

#calculate the 95% (HR) and 99.9% (EXT) home range area
HR <- round(sqrt((-2*log(0.05)*pi)*prey.SIG(mass)))
EXT <- round(sqrt((-2*log(0.0001)*pi)* prey.SIG(mass)))

#create circles to add to ggplot
HR_area <- circleFun(diameter = 2*HR)
EXT_area <- circleFun(diameter = 2*EXT)

#plot it all together
ggplot() +
  geom_raster(data = df_raster, aes(x = x, y = y, fill = calories)) +
  geom_vline(xintercept = xlines, color = "white", alpha = 0.5) +
  geom_hline(yintercept = ylines, color = "white", alpha = 0.5) +
  scale_fill_viridis_c() +
  geom_path(data = PREY_tracks[[1]], aes(x = x, y = y), color = "black", linewidth = 0.7, alpha = 0.8) +
  geom_path(dat = HR_area, aes(x,y), color = "#467378") +
  geom_path(dat = EXT_area, aes(x,y), color = "#68855C") + 
  coord_equal() +
  theme_minimal() 

#save the figure
#ggsave(file = 'sim_results/july16/figures/30000g_1000thlifespan_movepath.png', width = 8, height = 8, dpi = 900, bg = "white")

#........................................................................
# reduced diagnostics ----
#------------------------------------------------------------------------

#this section creates plots using the same data, filtered by generations to 
#investigate where changes in the lineages occur

lastgens <- prey_details_df %>%
  filter(generation >= max(generation) - 199)

p1 <-ggplot(lastgens, aes(x = generation, y = lv)) +
  stat_summary(fun = mean, geom = "line", col = "steelblue") +
  labs(title = "mean l_v over generations") +
  theme_minimal()

#tau_v vs generation
p2 <- ggplot(lastgens, aes(x = generation, y = tau_v)) +
  stat_summary(fun = mean, geom = "line", col = "steelblue") +
  labs(title = "mean tau_v over generations") +
  theme_minimal()

#tau_p vs generation
p3 <- ggplot(lastgens, aes(x = generation, y = tau_p)) +
  stat_summary(fun = mean, geom = "line", col = "steelblue") +
  labs(title = "mean tau_p over generations") +
  theme_minimal()

#lv vs patches
p4 <- ggplot(lastgens, aes(x = lv, y = patches)) +
  geom_point(col = "steelblue", alpha = 0.3) +
  labs(title = "lv vs patches") +
  theme_minimal()

#lv vs speed
p5 <- ggplot(lastgens, aes(x = lv, y = speed)) +
  geom_point(col = "steelblue", alpha = 0.3) +
  labs(title = "lv vs speed", y = "speed (m/s)") +
  theme_minimal()

#speed vs patches
p6 <- ggplot(lastgens, aes(x = speed, y = patches)) +
  geom_point(col = "steelblue", alpha = 0.3) +
  labs(title = "patches vs speed", x = "speed (m/s)") +
  theme_minimal()

#speed vs calories
p7 <- ggplot(lastgens, aes(x = speed, y = cal_net)) +
  geom_point(col = "steelblue", alpha = 0.3) +
  labs(title = "calories gained vs speed", x = "speed (m/s)") +
  theme_minimal()

#calories vs lv
p8 <- ggplot(lastgens, aes(x = lv, y = cal_net)) +
  geom_point(col = "steelblue", alpha = 0.3) +
  labs(title = "calories gained vs l_v") +
  theme_minimal()

#offspring vs lv
p9 <- ggplot(lastgens, aes(x = lv, y = offspring)) +
  geom_point(col = "steelblue", alpha = 0.3) +
  labs(title = "offspring vs l_v") +
  theme_minimal()

p10 <- ggplot(lastgens, aes(x = speed, y = offspring)) +
  geom_point(col = "steelblue", alpha = 0.3) +
  labs(title = "offspring vs speed", x = "speed (m/s)") +
  theme_minimal()

#create the final figure
FIG <- grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10, ncol = 3, nrow = 4, top = "100g, generations 500 to 1000")

#save the figure
#ggsave(FIG, file = '~/hdrive/GitHub/ballistic-movement/sim_results/sensitivity/100g_case_study/500_to_1000_gens.png', width = 12, height = 8, dpi = 900, bg = "white")


#-------------------------------------------------------------------------------
# comparing tau_p and sig ----
#-------------------------------------------------------------------------------

#tau_p able to evolve
load('sim_results/july16/5000000g_newsetup_prey_res.Rda')
load('sim_results/july16/5000000g_newsetup_prey_details.Rda')

# make data sets compatible with ggplot2
prey_res_df <- do.call(rbind, prey_res)
prey_details_df <- do.call(rbind, prey_details)

taup.sig <- ggplot(prey_details_df, aes(x = generation, y = (tau_p/sig))) +
  geom_point() +
  theme.qel()
print(taup.sig)

#held tau_p
load('sim_results/sensitivity/movement_params/5000000g_heldtaup2_prey_res.Rda')
load('sim_results/sensitivity/movement_params/5000000g_heldtaup2_prey_details.Rda')

# make data sets compatible with ggplot2
prey_res_df <- do.call(rbind, prey_res)
prey_details_df <- do.call(rbind, prey_details)

taup.sig2 <- ggplot(prey_details_df, aes(x = generation, y = (tau_p/sig))) +
  geom_point() +
  ylim(0, 0.2) +
  theme.qel()
print(taup.sig2)
