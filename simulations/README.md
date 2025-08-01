# Folder: simulations

## Purpose:

This directory contains all simulation outputs for the ballistic movement modeling
project. Simulations test how predator and prey adjust their movement strategies in response 
to environmental conditions, such as food distribution, predator presence, and movement cost.

## Structure:
- `prey_results/`: Outputs for simulations of prey across the mass spectrum, in the absense of predators
- `sensitivity/`: Simulations exploring the effects of varying parameters, like patches per home range area and kcalorie density.
- `supplementary/`: Simulation outputs for comparing the effect of metabolic costs on the relationship between speed and ballistic lengthscale, and its effect on calorie gain. 

## How to use:

Each subfolder contains its own `README.md` and associated outputs. 
To regenerate specific simulations or outputs, refer to the subfolders as appropriate.

## Notes:
- Large results files are stored as `.Rda` and named with parameter values.