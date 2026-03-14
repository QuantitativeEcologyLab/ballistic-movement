# script for making plot of ballistic length scale vs number of patches in the habitat

library(tidyverse)

cmbnd <- list.files(path = "simulations/prey_results/num_patches", 
                          pattern = "prey_details\\.Rds$", 
                          full.names = TRUE) %>% 
  map(~ {
    id_val <- as.numeric(str_extract(basename(.x), "\\d+"))
    
    data_list <- readRDS(.x) %>% map_dfr(., ~ mutate(.x, file_id = id_val))
    
    data_list <- data_list %>%
      group_by(file_id) %>%
      filter(generation <= max(generation - 500)) %>%
      summarise(mean_lv = mean(lv))
    
    return(data_list)
  }) %>%
  list_rbind()

cmbnd <- cmbnd %>% 
  filter(!file_id %in% c(1000, 30800, 40400))

mod <- glm(mean_lv ~ file_id, 
           data = cmbnd, 
           family = Gamma(link = "log"))

New_Data <- data.frame(file_id = seq(0, 55000, length.out = 100))

#generate predictions from GLM
preds <- predict(mod, newdata = New_Data, type = "link", se = TRUE)

New_Data$fit <- exp(preds$fit)
New_Data$lowerci <- exp(preds$fit - preds$se.fit * 1.96)
New_Data$upperci <- exp(preds$fit + preds$se.fit * 1.96)

ggplot() +
  geom_point(data = cmbnd, aes(x = file_id, y = mean_lv)) +
  geom_ribbon(data = New_Data, aes(ymin = lowerci, ymax = upperci, x = file_id), alpha = 0.3) +
  geom_line(data = New_Data, aes(x = file_id, y = fit), col = "red") +
  labs(x = "Number of Patches in Landscape", y = "Mean ballistic length scale") +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=12, family = "sans", face = "bold"),
        axis.title.x = element_text(size=12, family = "sans", face = "bold"),
        axis.text.y = element_text(size=10, family = "sans"),
        axis.text.x  = element_text(size=10, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 16, family = "sans", face = "bold"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
