# Folder: figures

## Purpose:

This folder contains visual analyses of simulation results.

## Subfolders:

`maintext/` contains figures used in the manuscript.

`sensivity/` contains figures used in the sensitivity analysis.

`supplementary/` contains figures used for the supplementary documents of the study.

## How these were generated:

These plots were generated using the following scripts:
```r
source("scripts/plotting.R")
source("scripts/dianostics.R")
source("scripts/cost_analsis.R")
source("simulations/sensitivity/ballistic_motion_sensitivity_script.R")