library(tidyverse)

# list files that are stored in "/simulations"
files <- list.files("simulations", pattern = "\\.Rds$",
                    recursive = TRUE, full.names = TRUE)

# collect relevant information on those files
manifest <- tibble(path = files) %>% 
  mutate(
    source = case_when(
      str_detect(path, "/archive/")      ~ "archive",
      str_detect(path, "/testing/")      ~ "testing",
      str_detect(path, "/sensitivity/")  ~ "sensitivity",
      str_detect(path, "/prey_results/") ~ "prey_results"),
    sweep     = str_extract(path, "calorie-variance|clustering|patches|npatches|variance"),
    replicate = if_else(str_detect(path, "duplicates"), "dup", "orig"),
    type = case_when(
      str_detect(basename(path), "_prey_details") ~ "prey_details",
      str_detect(basename(path), "_prey_res")     ~ "prey_res",
      TRUE                                        ~ "habitat"),
    param   = parse_number(basename(path)),
    size_mb = file.size(path) / 1e6,
    created = file.info(path)$mtime   
  )

runs <- manifest %>% 
  filter(type == "prey_details", !is.na(sweep), !is.na(param)) %>% 
  mutate(param = round(param, 3)) %>% 
  group_by(source, sweep, replicate, param) %>% 
  summarise(created = min(created), .groups = "drop") %>% 
  mutate(done = TRUE,
         seed_issue = if_else(created < as.POSIXct("2026-08-05"), "yes", "no"))

planned <- bind_rows(
  tibble(sweep = "calorie-variance", param = seq(0, 1, by = 0.025)),
  tibble(sweep = "clustering",       param = seq(1, 30, by = 0.5)),
  tibble(sweep = "patches",          param = seq(250, 5000, by = 50))) %>% 
  mutate(param = round(param, 3), in_plan = TRUE) %>% 
  crossing(replicate = c("orig", "dup", "trip"))

status <- planned %>% 
  full_join(filter(runs, source == "prey_results") %>% select(-source),
            by = c("sweep", "replicate", "param")) %>% 
  mutate(in_plan = coalesce(in_plan, FALSE),
         status = case_when(
           !in_plan   ~ "extra (not in plan)",
           is.na(done) ~ "todo",
           TRUE        ~ "done"))

status %>% count(sweep, replicate, status)

todo <- status %>% 
  filter(status == "todo")

# write to csv
write.csv(status, "simulations/manifest.csv")

write.csv(todo, "simulations/todo.csv")
