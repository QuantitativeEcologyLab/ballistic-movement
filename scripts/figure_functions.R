## import packages

# this script is designed to create diagnostic figure functions for the ballistic movement project

library(ggplot2) # for creating figures
library(rcartocolor) # for color palettes
library(ggthemes) # for custom themes

# for pred / prey colors
geyser <- carto_pal(7, "Geyser")

prey_col <- geyser[[1]]
pred_col <- geyser[[7]]

#----------------------------------------------------------------------
# create figure theme----
#----------------------------------------------------------------------

#theme is designed for saving with ggsave() with the following specs
# width = 6.86, height = 3.5, units = "in", dpi = 600

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
# prey diagnostics----
#----------------------------------------------------------------------

prey.gen.summary <- function(res_df, details_df) {
  
  # relative change in lv ~ gen
  prey_LV <- res_df$lv[1]
  res_df$rel.lv <- res_df$lv/prey_LV
  res_df$rel_var <- res_df$var / (prey_LV^2)
  res_df$rel_sd <- sqrt(res_df$rel_var)
  
  prey_summary <- details_df %>%
    group_by(generation) %>%
    summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop")
  
  # relative lv over generation
  p1 <- 
    ggplot(res_df, aes(x = generation, y = rel.lv)) +
    geom_line(color = prey_col, linewidth = 1) +
    geom_hline(yintercept = 1, color = "grey30", linetype = "dashed") +
    geom_ribbon(aes(ymin = rel.lv - rel_sd,
                    ymax = rel.lv + rel_sd),
                fill = prey_col, alpha = 0.3) +
    labs(y = "relative change in lv", x = "generation") +
    theme_few()
  
  # lv ~ generation
  p2 <- 
    ggplot(prey_summary, aes(x = generation, y = lv)) +
    geom_line(col = prey_col, linewidth = 1) +
    labs(x = "generation", y = "lv (m)") +
    theme_few()
  
  # tau_p ~ generation
  p3 <-
    ggplot(prey_summary, aes(x = generation, y = tau_p)) +
    geom_line(col = prey_col, linewidth = 1) +
    labs(x = "generation", y = "tau_p") +
    theme_few()
  
  # tau_v ~ generation
  p4 <-
    ggplot(prey_summary, aes(x = generation, y = tau_v)) +
    geom_line(col = prey_col, linewidth = 1) +
    labs(x = "generation", y = "tau_v") +
    theme_few()
  
  # sig ~ generation
  p5 <- 
    ggplot(prey_summary, aes(x = generation, y = sig)) +
    geom_line(col = prey_col, linewidth = 1) +
    labs(x = "generation", y = "sig") +
    theme_few()
  
  # encounters ~ generation
  p6 <- 
    ggplot(prey_summary, aes(x = generation, y = patches)) +
    geom_line(col = prey_col, linewidth = 1) +
    labs(x = "generation", y = "patches visited") +
    theme_few()
  
  # cal_net ~ generation
  p7 <- 
    ggplot(prey_summary, aes(x = generation, y = cal_net)) +
    geom_line(col = prey_col, linewidth = 1) +
    labs(x = "generation", y = "net calories") +
    theme_few()
  
  # speed ~ generation
  p8 <- 
    ggplot(prey_summary, aes(x = generation, y = speed)) +
    geom_line(col = prey_col, linewidth = 1) +
    labs(x = "generation", y = "speed (m/s)") +
    theme_few()
  
  # offspring ~ generation
  p9 <- 
    ggplot(prey_summary, aes(x = generation, y = offspring)) +
    geom_line(col = prey_col, linewidth = 1) +
    labs(x = "generation", y = "offspring") +
    theme_few()
  
  plot <- grid.arrange(p1, p2, p3,
                       p5, p4, p6, 
                       p7, p8, p9,
                       ncol = 3)
  
  return(plot)
}
  
# for comparing other variables
prey.var.summary <- function(details_df){
  # number of patches ~ lv
  p1 <- 
    ggplot(details_df, aes(x = lv, y = patches)) +
    geom_point(col = prey_col, alpha = 0.3) +
    labs(x = "lv (m)", y = "patches visited") +
    theme_few()
  
  # speed ~ lv
  p2 <- 
    ggplot(details_df, aes(x = lv, y = speed)) +
    geom_point(col = prey_col, alpha = 0.3) +
    labs(x = "lv (m)", y = "speed (m/s)") +
    theme_few()
  
  # offspring ~ lv
  p3 <- 
    ggplot(details_df, aes(x = lv, y = offspring)) +
    geom_point(col = prey_col, alpha = 0.3) +
    labs(x = "lv (m)", y = "offspring") +
    theme_few()
  
  # cal_net ~ lv
  p4 <- 
    ggplot(details_df, aes(x = lv, y = cal_net)) +
    geom_point(col = prey_col, alpha = 0.3) +
    labs(x = "lv (m)", y = "net calories") +
    theme_few()
  
  # speed ~ cal_net 
  p5 <- 
    ggplot(details_df, aes(x = speed, y = cal_net)) +
    geom_point(col = prey_col, alpha = 0.3) +
    labs(x = "speed (m/s)", y = "net calories") +
    theme_few()
  
  # number of patches ~ speed
  p6 <- 
    ggplot(details_df, aes(x = speed, y = patches)) +
    geom_point(col = prey_col, alpha = 0.3) +
    labs(x = "speed (m/s)", y = "patches visited") +
    theme_few()
  
  # offspring ~ speed
  p7 <- 
    ggplot(details_df, aes(x = speed, y = offspring)) +
    geom_point(col = prey_col, alpha = 0.3) +
    labs(x = "speed (m/s)", y = "offspring") +
    theme_few()
  
  # encounters ~ offspring
  p8 <- 
    ggplot(details_df, aes(x = patches, y = offspring)) +
    geom_point(col = prey_col, alpha = 0.3) +
    labs(x = "patches visited", y = "offspring") + 
    theme_few()
  
  # offspring ~ cal_net
  p9 <- 
    ggplot(details_df, aes(x = cal_net, y = offspring)) +
    geom_point(col = prey_col, alpha = 0.3) +
    labs(x = "net calories", y = "offspring") +
    theme_few()
  
  # speed ~ tau_p
  p10 <- 
    ggplot(details_df, aes(y = speed, x = tau_p)) +
    geom_point(col = prey_col, alpha = 0.3) +
    labs(y = "speed (m/s)", x = "tau_p (s)") +
    theme_few()
  
  plot <- grid.arrange(p1, p2,
                       p3, p4, 
                       p5, p6, 
                       p7, p8,
                       p9, p10,
                       ncol = 2)
  
  return(plot)
}


#----------------------------------------------------------------------
# predator diagnostics----
#----------------------------------------------------------------------

# for comparing changes over the number of generations
pred.gen.summary <- function(res_df, details_df){
  
  # relative change in lv ~ gen
  pred_LV <- res_df$lv[1]
  res_df$rel.lv <- res_df$lv/pred_LV
  res_df$rel_var <- res_df$var / (pred_LV^2)
  res_df$rel_sd <- sqrt(res_df$rel_var)
  
  pred_summary <- details_df %>%
    group_by(generation) %>%
    summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop")
  
  # relative lv over generation
  p1 <- 
    ggplot(res_df, aes(x = generation, y = rel.lv)) +
    geom_line(color = pred_col, linewidth = 1) +
    geom_hline(yintercept = 1, color = "grey30", linetype = "dashed") +
    geom_ribbon(aes(ymin = rel.lv - rel_sd,
                    ymax = rel.lv + rel_sd),
                fill = pred_col, alpha = 0.3) +
    labs(y = "relative change in lv", x = "generation") +
    theme_few()
  
  # lv ~ generation
  p2 <- 
    ggplot(pred_summary, aes(x = generation, y = lv)) +
    geom_line(col = pred_col, linewidth = 1) +
    labs(x = "generation", y = "lv (m)") +
    theme_few()
  
  # tau_p ~ generation
  p3 <-
    ggplot(pred_summary, aes(x = generation, y = tau_p)) +
    geom_line(col = pred_col, linewidth = 1) +
    labs(x = "generation", y = "tau_p") +
    theme_few()
  
  # tau_v ~ generation
  p4 <-
    ggplot(pred_summary, aes(x = generation, y = tau_v)) +
    geom_line(col = pred_col, linewidth = 1) +
    labs(x = "generation", y = "tau_v") +
    theme_few()
  
  # sig ~ generation
  p5 <- 
    ggplot(pred_summary, aes(x = generation, y = sig)) +
    geom_line(col = pred_col, linewidth = 1) +
    labs(x = "generation", y = "sig") +
    theme_few()

  # encounters ~ generation
  p6 <- 
    ggplot(pred_summary, aes(x = generation, y = encounters)) +
    geom_line(col = pred_col, linewidth = 1) +
    labs(x = "generation", y = "prey encountered") +
    theme_few()
  
  # cal_net ~ generation
  p7 <- 
    ggplot(pred_summary, aes(x = generation, y = cal_net)) +
    geom_line(col = pred_col, linewidth = 1) +
    labs(x = "generation", y = "net calories") +
    theme_few()
  
  # speed ~ generation
  p8 <- 
    ggplot(pred_summary, aes(x = generation, y = speed)) +
    geom_line(col = pred_col, linewidth = 1) +
    labs(x = "generation", y = "speed (m/s)") +
    theme_few()
  
  # offspring ~ generation
  p9 <- 
    ggplot(pred_summary, aes(x = generation, y = offspring)) +
    geom_line(col = pred_col, linewidth = 1) +
    labs(x = "generation", y = "offspring") +
    theme_few()
  
  plot <- grid.arrange(p1, p2, p3,
                       p5, p4, p6, 
                       p7, p8, p9,
                       ncol = 3)
  
  return(plot)
}

# for comparing other variables
pred.var.summary <- function(details_df){
  # number of encounters ~ lv
  p1 <- 
    ggplot(details_df, aes(x = lv, y = encounters)) +
    geom_point(col = pred_col, alpha = 0.3) +
    labs(x = "lv (m)", y = "prey encountered") +
    theme_few()
  
  # speed ~ lv
  p2 <- 
    ggplot(details_df, aes(x = lv, y = speed)) +
    geom_point(col = pred_col, alpha = 0.3) +
    labs(x = "lv (m)", y = "speed (m/s)") +
    theme_few()
  
  # offspring ~ lv
  p3 <- 
    ggplot(details_df, aes(x = lv, y = offspring)) +
    geom_point(col = pred_col, alpha = 0.3) +
    labs(x = "lv (m)", y = "offspring") +
    theme_few()
  
  # cal_net ~ lv
  p4 <- 
    ggplot(details_df, aes(x = lv, y = cal_net)) +
    geom_point(col = pred_col, alpha = 0.3) +
    labs(x = "lv (m)", y = "net calories") +
    theme_few()
  
  # speed ~ cal_net 
  p5 <- 
    ggplot(details_df, aes(x = speed, y = cal_net)) +
    geom_point(col = pred_col, alpha = 0.3) +
    labs(x = "speed (m/s)", y = "net calories") +
    theme_few()
  
  # number of encounters ~ speed
  p6 <- 
    ggplot(details_df, aes(x = speed, y = encounters)) +
    geom_point(col = pred_col, alpha = 0.3) +
    labs(x = "speed (m/s)", y = "prey encountered") +
    theme_few()
  
  # offspring ~ speed
  p7 <- 
    ggplot(details_df, aes(x = speed, y = offspring)) +
    geom_point(col = pred_col, alpha = 0.3) +
    labs(x = "speed (m/s)", y = "offspring") +
    theme_few()
  
  # encounters ~ offspring
  p8 <- 
    ggplot(details_df, aes(x = encounters, y = offspring)) +
    geom_point(col = pred_col, alpha = 0.3) +
    labs(x = "prey encountered", y = "offspring") + 
    theme_few()
  
  # offspring ~ cal_net
  p9 <- 
    ggplot(details_df, aes(x = cal_net, y = offspring)) +
    geom_point(col = pred_col, alpha = 0.3) +
    labs(x = "net calories", y = "offspring") +
    theme_few()
  
  # speed ~ tau_p
  p10 <- 
    ggplot(details_df, aes(y = speed, x = tau_p)) +
    geom_point(col = pred_col, alpha = 0.3) +
    labs(y = "speed (m/s)", x = "tau_p (s)") +
    theme_few()
  
  plot <- grid.arrange(p1, p2,
                       p3, p4, 
                       p5, p6, 
                       p7, p8,
                       p9, p10,
                       ncol = 2)
  
  return(plot)
}

#----------------------------------------------------------------------
# project update figures----
#----------------------------------------------------------------------

prey.proj.update <- function(res_df, details_df) {
  # relative change in lv ~ gen
  PREY_LV <- res_df$lv[1]
  res_df$rel.lv <- res_df$lv/PREY_LV
  res_df$rel_var <- res_df$var / (PREY_LV^2)
  res_df$rel_sd <- sqrt(res_df$rel_var)
  
  prey_summary <- details_df %>%
    group_by(generation) %>%
    summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop")
  
  p1 <- 
    ggplot(res_df, aes(x = generation, y = rel.lv)) +
    ggtitle("A") +
    geom_line() +
    geom_hline(yintercept = 1, color = 'grey30', linetype = "dashed") +
    geom_ribbon(aes(ymin = rel.lv - rel_sd,
                    ymax = rel.lv + rel_sd),
                alpha = 0.3) +
    labs(y = "Relative Change in lv", x = "Generation") +
    theme.qel()
  
  p2 <- 
    ggplot(prey_summary, aes(x = generation, y = tau_v)) +
    ggtitle("B") +
    geom_line() +
    labs(x = "Generation", y = expression(bold(tau[v]))) +
    theme.qel()
  
  p3 <- 
    ggplot() +
    ggtitle("C") +
    geom_point(data = details_df, aes(x = generation, y = tau_p), 
               alpha = 0.2, color = "grey70", size = 0.5) +
    geom_line(data = prey_summary, aes(x = generation, y = tau_p)) +
    labs(x = "Generation", y = expression(bold(tau[p]))) +
    theme.qel()
  
  p4 <- 
    ggplot(prey_summary, aes(x = generation, y = sig)) +
    ggtitle("D") +
    geom_line() +
    labs(x = "Generation", y = expression(bold(sigma))) +
    theme.qel()
  
  p5 <- 
    ggplot(details_df, aes(x = lv, y = cal_net, color = generation)) +
    ggtitle("E") +
    geom_point(size = 0.2) +
    labs(x = expression(bold(l[v])), y = "Net calories") +
    theme.qel()
  
  p6 <-
    ggplot(details_df, aes(x = speed, y = cal_net, color = generation)) +
    ggtitle("F") +
    geom_point(size = 0.2) +
    labs(x = "Speed (m/s)", y = "Net calories") +
    theme.qel()
  
  p7 <- 
    ggplot(details_df, aes(x = lv, y = speed, color = generation)) +
    ggtitle("G") +
    geom_point(size = 0.2) +
    labs(x = expression(bold(l[v])), y = "Speed (m/s)") +
    theme.qel()
  
  p8 <- 
    ggplot(details_df, aes(x = lv, y = offspring, color = generation)) +
    ggtitle("F") +
    geom_point(size = 0.2) +
    labs(x = expression(bold(l[v])), y = "Offspring") +
    theme.qel()
  
  final <- grid.arrange(p1, p2, p3, p4,
                        p5, p6, p7, p8,
                        ncol = 4)
  return(final)
}

