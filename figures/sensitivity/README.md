# figures/sensitivity

## Purpose:

This folder contains figures for the sensitivity analysis of the ballistic lengthscale 
simulation models.

## How these were generated
These plots were generated using the following scripts:
```r
source("simulations/sensitivity/ballistic_motion_sensitivity_script.R")
```

`<mass>_movetracks.png` are movement tracks from sensitivity simulations testing the number of patches
per HR area needed to detect the maximum number of encounters.

`kcal_per_patch` contains figures from testing the bare minimum number of kcalories per 
patch needed to prevent extinction in each prey mass value to be tested. The value from
these analyses were used as a baseline, to detemine the lowest kcalorie value possible to prevent
extinction in a full simulation setup. 