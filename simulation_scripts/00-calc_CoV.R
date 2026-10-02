## script to calculate CoV of landscape and save to Rds data file
## needed for the first set of caloric-variance simulations, which did not have the CoV saved
## should not be needed for simulations run with `03-variance-simulation.R`

library(spatstat)
library(tidyverse)

setwd("~/hdrive/GitHub/ballistic-movement/simulations/prey_results/calorie-variance/duplicates")

data_dir <- "prey_details"
habitat_dir <- "habitats"

habitat_files <- list.files(
  habitat_dir,
  pattern = "^[0-9]+\\.[0-9]+cv\\.Rds$",
  full.names = TRUE
)

data_files <- list.files(data_dir, pattern = "^[0-9]+\\.[0-9]+cv.*\\.Rds$")
data_nums <- as.numeric(sub("^([0-9]+\\.[0-9]+)cv.*$", "\\1", data_files))

for(f in habitat_files){
  habitat_path <- file.path(f)
  
  habitat_num <- as.numeric(sub("^.*?([0-9]+\\.[0-9]+)cv.*$", "\\1", f))
  match_idx <- which(abs(data_nums - habitat_num) < 1e-8)
  
  data_path <- file.path(data_dir, data_files[match_idx])
  
  pp <- readRDS(habitat_path)
  m <- marks(pp)
  
  cov_marks <- sd(m) / mean(m)
  
  data <-  readRDS(data_path)
  
  data <- data %>% 
    bind_rows() %>% 
    mutate(cov = cov_marks)
  
  data <- data %>% 
    group_split(generation)
  
  saveRDS(data, data_path)
  
  message(f, ": CoV =", round(cov_marks, 4), "-- updated")
}
