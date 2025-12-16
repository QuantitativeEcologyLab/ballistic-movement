# this script is designed to create the theme for ggplots in the prey search behaviour simulation project

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