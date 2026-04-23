## A Computational Approach to Identifying Evolutionarily Stable Strategies in Mammalian Search Behaviour

**Authors:** L. A. Terpsma, R. Martinez-Garcia, C. H. Fleming, W. F. Fagan, J. M. Calabrese, M. J. Noonan

Understanding how organisms navigate their environment under landscape changes is 
increasingly important under rapid human-induced changes, yet predictive models 
often fail to incorporate ecological constraints which govern animal movement. This 
study aims to investigate how mammalian movement strategies adapt to variation in resource 
distribution using a computational approach. We developed a stochastic system, grounded 
in allometric scaling relationships, where individual movement behaviour emerges from interaction 
with variable resource landscapes and energetic costs, evolving over evolutionary time. 
Our results indicate that search strategies in a population stabilise over time, with 
emergent optima dependent on landscape composition. A fundamental trade-off between search 
efficiency and movement speed, emerges as a result of energetic constraints. Higher resource 
availability reduced ballistic length scale ($l_v$) (Gamma GLM: $\beta = -1.35 \times 10^{-4}$, 
$p = 0.00448$) and increased movement speed non-linearly (Nonlinear least squares (Michaelis-Menten): 
$a = 5.22$, $p < 0.001$, $b = 1079$, $p < 0.001$), suggesting individuals rely less on slow 
long-distance movement when resources are energetically rich. In contrast, when resources are 
highly clustered, we found $l_v$ decreased (Gamma GLM: $\beta = -0.00794$, $p = 0.183$), 
along with speed (Gamma GLM: $\beta = -0.0136$, $p < 0.001$), suggesting search strategy 
relies on maintaining position within a resource cluster.  Resource unpredictability has 
found to have no significant effect on movement, but further investigation is warranted in 
this domain. Overall, our findings highlight the importance of incorporating ecological 
constraints into predictive movement models to better capture biologically plausible dynamics. 
This framework provides a generalisable approach for predicting how environmental change may 
influence animal movement, with potential applications in conservation planning and broader 
relevance to understanding movement processes across various ecological systems.

### Repository Contents 

Each folder contains another README file which describes each component in greater detail. 

* `simulation_scripts/` contains the R script files for the necessary functions and simulation code.
* `figure_scripts/` contains functions and scripts necessary to create all figures presented in the figures folder.
* `figures/` contains all figures of simulation results.
* `presentations/` contains files used for presentation of the project, such as posters and slideshows.

### Start Here

The below provides details on the workflow needed to reproduce simulations. The R files listed below are found in the `simulation_scripts` folder.

* `01-prey-functions.R`: functions for generating prey only simulations
* `02-prey-simulation.R`: workflow for simulating evolution of prey search behaviour over evolutionary timescales.

### R Environment and Packages

Simulations, models, and generation of all figures were conducted in the R statistical package (v.4.5.2 R Core Team 2025) using the RANN (v. 2.6.2), spatstat.random (v.3.4-2), 
spatstat.geom (v. 3.6-0), ctmm (v. 1.2.0), extraDistr (v. 1.10-0), mgcv (v. 1.9-3), tictoc (v. 1.2.1), tidyverse (v. 2.0.0), gridExtra (v. 2.3), viridis (v. 0.6.5), 
propagate (v. 1.1-0), patchwork (v. 1.3.2), and scico (v. 1.5.0.9000) packages. 

Auguie B (2017). _gridExtra: Miscellaneous Functions for "Grid" Graphics_. doi:10.32614/CRAN.package.gridExtra
<https://doi.org/10.32614/CRAN.package.gridExtra>, R package version 2.3, <https://CRAN.R-project.org/package=gridExtra>.

Baddeley A, Rubak E, Turner R (2015). _Spatial Point Patterns: Methodology and Applications with R_. Chapman and Hall/CRC Press, London. ISBN 9781482210200,
<https://www.routledge.com/Spatial-Point-Patterns-Methodology-and-Applications-with-R/Baddeley-Rubak-Turner/p/book/9781482210200/>.

Fleming CH, Calabrese JM (2023). _ctmm: Continuous-Time Movement Modeling_. doi:10.32614/CRAN.package.ctmm <https://doi.org/10.32614/CRAN.package.ctmm>, R package
version 1.2.0, <https://CRAN.R-project.org/package=ctmm>.

Izrailev S (2024). _tictoc: Functions for Timing R Scripts, as Well as Implementations of "Stack" and "StackList" Structures_. doi:10.32614/CRAN.package.tictoc
<https://doi.org/10.32614/CRAN.package.tictoc>, R package version 1.2.1, <https://CRAN.R-project.org/package=tictoc>.

Jefferis G, Kemp SE, Arya S, Mount D (2024). _RANN: Fast Nearest Neighbour Search (Wraps ANN Library) Using L2 Metric_. doi:10.32614/CRAN.package.RANN
<https://doi.org/10.32614/CRAN.package.RANN>, R package version 2.6.2, <https://CRAN.R-project.org/package=RANN>.

Pedersen T, Crameri F (2025). _scico: Colour Palettes Based on the Scientific Colour-Maps_. R package version 1.5.0.9000, commit
e94d08c334c8de7ba5dd0c405baeb578a5d2651c, <https://github.com/thomasp85/scico>.

Pedersen T (2025). _patchwork: The Composer of Plots_. doi:10.32614/CRAN.package.patchwork
<https://doi.org/10.32614/CRAN.package.patchwork>, R package version 1.3.2, <https://CRAN.R-project.org/package=patchwork>.

Simon Garnier, Noam Ross, Robert Rudis, Antônio P. Camargo, Marco Sciaini, and Cédric Scherer (2024). viridis(Lite) - Colorblind-Friendly Color Maps for R. viridis package version 0.6.5.

Spiess A (2026). _propagate: Propagation of Uncertainty_. doi:10.32614/CRAN.package.propagate
  <https://doi.org/10.32614/CRAN.package.propagate>, R package version 1.1-0, <https://CRAN.R-project.org/package=propagate>.

Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, Grolemund G, Hayes A, Henry L, Hester J, Kuhn M, Pedersen TL, Miller E, Bache SM, Müller K, Ooms J,
Robinson D, Seidel DP, Spinu V, Takahashi K, Vaughan D, Wilke C, Woo K, Yutani H (2019). “Welcome to the tidyverse.” _Journal of Open Source Software_, *4*(43), 1686.
doi:10.21105/joss.01686 <https://doi.org/10.21105/joss.01686>.

Wolodzko T (2023). _extraDistr: Additional Univariate and Multivariate Distributions_. doi:10.32614/CRAN.package.extraDistr
<https://doi.org/10.32614/CRAN.package.extraDistr>, R package version 1.10.0, <https://CRAN.R-project.org/package=extraDistr>.

Wood, S.N. (2011) Fast stable restricted maximum likelihood and marginal likelihood estimation of semiparametric generalized linear models. Journal of the Royal
Statistical Society (B) 73(1):3-36