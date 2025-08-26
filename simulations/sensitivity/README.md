# simulations/sensitivity

## Purpose:

This folder contains the outputs from senstivity analyses that vary ecological parameters to assess
their impact on prey fitness, movement parameters, and model behaviour.

## Parameters tested:

- Patches present in the 95% Home Range Area
- kcal per patch density
- the effect of the sampling interval on the number of patch encounters detected


## Structure:

- `ballistic_motion_sensitivity_script.R`: Main script to generate all sensitivity simulations.
- `ballistic_motion_senstivity_report.R`: rMarkdown file to generate figures and create `.pdf` report for manuscript..
- `data/`: contains `.Rda` files storing simulation outputs by mass and test conducted.
- `kcal_per_patch`: contains figures from kcalorie per patch density analysis for masses to test at full scale.

## How to reproduce:

To rerun all sensitivity analyses:
```r
source("simulations/sensitivity/ballistic_motion_sensitivity_script.R") 
```
