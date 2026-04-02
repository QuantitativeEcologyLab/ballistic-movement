# this script is designed to create diagnostic figure functions for the ballistic movement project
# this allows for quick visualization of many parameters and their relationships

#import necessary packages

library(ggplot2) # for creating figures
library(paletteer) # for continuous color palettes
library(ggthemes) # for custom themes
library(gridExtra) # grid combining figures
library(grid) # for combining figures
library(dplyr) # for data cleaning
library(patchwork) # for combined legend
library(viridis) #color palette

#----------------------------------------------------------------------
# simulation check function----
#----------------------------------------------------------------------

simple.diag <- function(details_df) {
  p1 <- 
    ggplot() +
    geom_point(data = details_df, aes(x = generation, y = lv), 
               alpha = 0.2, color = "grey70", size = 0.1) +
    stat_summary(data = details_df, aes(x = generation, y = lv), fun = "mean", geom = "line") +
    labs(y = expression(bold(l[v])), x = "Generation") +
    theme_bw()
  
  p5 <- 
    ggplot(details_df, aes(x = lv, y = cal_net, color = generation)) +
    geom_point(size = 0.1) +
    scale_color_viridis_c(option = "magma") +
    labs(x = expression(bold(l[v])), y = "Net calories") +
    theme_bw()
  
  p6 <-
    ggplot(details_df, aes(x = speed, y = cal_net, color = generation)) +
    geom_point(size = 0.1) +
    scale_color_viridis_c(option = "magma") +
    labs(x = "Speed (m/s)", y = "Net calories") +
    theme_bw()
  
  p7 <- 
    ggplot(details_df, aes(x = lv, y = speed, color = generation)) +
    geom_point(size = 0.1) +
    scale_color_viridis_c(option = "magma") +
    labs(x = expression(bold(l[v])), y = "Speed (m/s)") +
    theme_bw()
  
  p2 <-
    ggplot() +
    geom_point(data = details_df, aes(x = generation, y = patches), 
               alpha = 0.2, color = "grey70", size = 0.1) +
    stat_summary(data = details_df, aes(x = generation, y = patches), fun = "mean", geom = "line") +
    labs(y = "Number of Patches Encountered", x = "Generation") +
    theme_bw()
  
  p3 <-
    ggplot()+
    geom_point(data=details_df, aes(x = generation, y = offspring),
               alpha = 0.2, color = "grey70", size = 0.1) + 
    stat_summary(data = details_df, aes(x = generation, y = offspring), fun = "mean", geom = "line") +
    labs(x = "Generation", y = "Offspring per Individual") +
    theme_bw()
  
  
  final <- p1 + p2 + p3 + p5 + p6 + p7 + plot_layout(ncol = 2)

  return(final)
}

explore.gen <- function(details_df) {
  p1 <- 
    ggplot() +
    geom_point(data = details_df, aes(x = generation, y = lv), 
               alpha = 0.2, color = "grey70", size = 0.1) +
    stat_summary(data = details_df, aes(x = generation, y = lv), fun = "mean", geom = "line") +
    labs(y = expression(bold(l[v])), x = "Generation") +
    theme_bw()
  
  p2 <- 
    ggplot() +
    geom_point(data = details_df, aes(x = generation, y = tau_v), 
               alpha = 0.2, color = "grey70", size = 0.1) +
    stat_summary(data = details_df, aes(x = generation, y = tau_v), fun = "mean", geom = "line") +
    theme_bw()
  
  p3 <- 
    ggplot() +
    geom_point(data = details_df, aes(x = generation, y = tau_p), 
               alpha = 0.2, color = "grey70", size = 0.1) +
    stat_summary(data = details_df, aes(x = generation, y = tau_p), fun = "mean", geom = "line") +
    theme_bw()
  
  p4 <- 
    ggplot() +
    geom_point(data = details_df, aes(x = generation, y = speed),
               alpha = 0.2, color = "grey70", size = 0.1) +
    stat_summary(data = details_df, aes(x = generation, y = speed), fun = "mean", geom = "line") +
    theme_bw()
  
  p5 <-
    ggplot() + 
    geom_point(data = details_df, aes(x = generation, y = offspring)) + 
    stat_summary(data = details_df, aes(x = generation, y = speed), fun = "mean", geom = "line") +
    theme_bw()
  
  p6 <- 
    ggplot() +
    geom_point(data = details_df, aes(x = generation, y = patches), 
               alpha = 0.2, color = "grey70", size = 0.1) +
    stat_summary(data = details_df, aes(x = generation, y = patches), fun = "mean", geom = "line") +
    labs(y = "Number of Patches Encountered", x = "Generation") +
    theme_bw()
  
  final <- p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 3)
  
  return(final)
}

explore.var <- function(details_df) {
  p1 <- 
    ggplot(details_df, aes(x = lv, y = cal_net, color = generation)) +
    geom_point(size = 0.1) +
    scale_color_viridis_c(option = "magma") +
    labs(x = expression(bold(l[v])), y = "Net calories") +
    theme_bw()
  
  p2 <-
    ggplot(details_df, aes(x = speed, y = cal_net, color = generation)) +
    geom_point(size = 0.1) +
    scale_color_viridis_c(option = "magma") +
    labs(x = "Speed (m/s)", y = "Net calories") +
    theme_bw()
  
  p3 <- 
    ggplot(details_df, aes(x = lv, y = speed, color = generation)) +
    geom_point(size = 0.1) +
    scale_color_viridis_c(option = "magma") +
    labs(x = expression(bold(l[v])), y = "Speed (m/s)") +
    theme_bw()
  
  final <- p1 + p2 + p3 + plot_layout(guides = "collect", ncol = 3)
  
  return(final)
}
