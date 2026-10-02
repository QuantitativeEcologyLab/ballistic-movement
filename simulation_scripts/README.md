This folder contains all scripts used for simulations and analyses. Details of each script can be found below.
- `00-calc-CoV.R`: uses `habitat` file from caloric variance simulations to calculate the CoV of patches and adds to simulation save file for any simulations that did not include this in the original file.
- `00-prey-functions.R`: contains the functions necessary for running the simulations. Functions calculate movement parameters, habitat, sampling interval, foraging, and fitness based on body size (in grams)
- `01-patches-simulation.R`: contains code to run resource abundance simulations.
- `02-clustering-simulation. R`: contains code to run resource distribution simulations. 
- `03-variance-simulation.R`: contains code to run resource unpredictability simulations. 
- `04-resume-preysim.R`: contains code to resume a simulation in the case of failure or the desire to extend simulation length without restarting. 
- `05-samp-int-sensitivity.R`: contains code for testing the required sampling interval to capture all patches encountered along a movement track. Setting the appropriate sampling interval prevents the animal entity from missing the consumption of a patch that was indeed encountered.