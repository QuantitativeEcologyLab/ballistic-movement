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
  prey_summary <- details_df %>%
    group_by(generation) %>%
    summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop")
  
  p1 <- 
    ggplot() +
    ggtitle("A") +
    geom_point(data = details_df, aes(x = generation, y = lv), 
               alpha = 0.2, color = "grey70", size = 0.1) +
    geom_line(data = prey_summary, aes(x = generation, y = lv)) +
    labs(y = "Ballistic Length Scale", x = "Generation") +
    theme_bw()
  
  p5 <- 
    ggplot(details_df, aes(x = lv, y = cal_net, color = generation)) +
    ggtitle("B") +
    geom_point(size = 0.1) +
    scale_color_viridis_c(option = "magma") +
    labs(x = expression(bold(l[v])), y = "Net calories") +
    theme_bw()
  
  p6 <-
    ggplot(details_df, aes(x = speed, y = cal_net, color = generation)) +
    ggtitle("C") +
    geom_point(size = 0.1) +
    scale_color_viridis_c(option = "magma") +
    labs(x = "Speed (m/s)", y = "Net calories") +
    theme_bw()
  
  p7 <- 
    ggplot(details_df, aes(x = lv, y = speed, color = generation)) +
    ggtitle("D") +
    geom_point(size = 0.1) +
    scale_color_viridis_c(option = "magma") +
    labs(x = expression(bold(l[v])), y = "Speed (m/s)") +
    theme_bw()
  
  final <- p1 + p5 + p6 + p7 + plot_layout(ncol = 2)

  return(final)
}
