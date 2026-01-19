# A Computational Approach to Identifying Evolutionarily Stable Strategies in Mammalian Search Behaviour

**Authors:** L. A. Terpsma, R. Martinez-Garcia, C. H. Fleming, W. F. Fagan, J. M. Calabrese, M. J. Noonan

This repository contains the code supporting the results of Lynndsay Terpsma's 2025 URA Project and Honours thesis, which investigates 
the evolutionary trends in ballistic length scale under spatial and richness variation in landscape resources. 
This project focuses on prey terrestrial mammals, evaluating behaviour across the body mass spectrum. 

# Repository Contents 

Each folder contains another README file which describes each component in greater detail. 

* `simulation_scripts/` contains the R script files for the necessary functions and simulation code.
* `figure_scripts/` contains functions and scripts necessary to create all figures presented in the figures folder.
* `figures/` contains all figures of simulation results.
* `presentations/` contains files used for presentation of the project, such as posters and slideshows.

# Start Here

The below provides details on the workflow needed to reproduce simulations. The R files listed below are found in the `simulation_scripts` folder.

* `01_prey_functions.R`: functions for generating prey only simulations
* `02_prey_simulation.R`: workflow for simulating evolution of prey search behaviour over evolutionary timescales.

# R Enviroment and Packages

Simulations, models, and generation of all figures were conducted in the R statistical package (v.4.5.2 R Core Team 2025) using the RANN (v. 2.6.2), spatstat.random (v.3.4-2), 
spatstat.geom (v. 3.6-0), ctmm (v. 1.2.0), extraDistr (v. 1.10-0), mgcv (v. 1.9-3), tictoc (v. 1.2.1), tidyverse (v. 2.0.0) packages. 

Baddeley A, Rubak E, Turner R (2015). _Spatial Point Patterns: Methodology and Applications with R_. Chapman and Hall/CRC Press, London. ISBN 9781482210200,
<https://www.routledge.com/Spatial-Point-Patterns-Methodology-and-Applications-with-R/Baddeley-Rubak-Turner/p/book/9781482210200/>.

Fleming CH, Calabrese JM (2023). _ctmm: Continuous-Time Movement Modeling_. doi:10.32614/CRAN.package.ctmm <https://doi.org/10.32614/CRAN.package.ctmm>, R package
version 1.2.0, <https://CRAN.R-project.org/package=ctmm>.

Izrailev S (2024). _tictoc: Functions for Timing R Scripts, as Well as Implementations of "Stack" and "StackList" Structures_. doi:10.32614/CRAN.package.tictoc
<https://doi.org/10.32614/CRAN.package.tictoc>, R package version 1.2.1, <https://CRAN.R-project.org/package=tictoc>.

Jefferis G, Kemp SE, Arya S, Mount D (2024). _RANN: Fast Nearest Neighbour Search (Wraps ANN Library) Using L2 Metric_. doi:10.32614/CRAN.package.RANN
<https://doi.org/10.32614/CRAN.package.RANN>, R package version 2.6.2, <https://CRAN.R-project.org/package=RANN>.

Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, Grolemund G, Hayes A, Henry L, Hester J, Kuhn M, Pedersen TL, Miller E, Bache SM, Müller K, Ooms J,
Robinson D, Seidel DP, Spinu V, Takahashi K, Vaughan D, Wilke C, Woo K, Yutani H (2019). “Welcome to the tidyverse.” _Journal of Open Source Software_, *4*(43), 1686.
doi:10.21105/joss.01686 <https://doi.org/10.21105/joss.01686>.

Wolodzko T (2023). _extraDistr: Additional Univariate and Multivariate Distributions_. doi:10.32614/CRAN.package.extraDistr
<https://doi.org/10.32614/CRAN.package.extraDistr>, R package version 1.10.0, <https://CRAN.R-project.org/package=extraDistr>.

Wood, S.N. (2011) Fast stable restricted maximum likelihood and marginal likelihood estimation of semiparametric generalized linear models. Journal of the Royal
Statistical Society (B) 73(1):3-36