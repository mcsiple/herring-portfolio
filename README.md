# herring-portfolio
Multivariate autoregressive state-space (MARSS) models for population structure, correlations between subpopulations, and portfolio effects in Puget Sound Pacific herring.

This folder contains two files with code:
1) The model structures that we fit for the first part of the paper. 
The data used to fit the model are from WDFW surveys of egg biomass. They are available from the WDFW from their most recent stock status report: http://wdfw.wa.gov/publications/01628/

2) The model structures we used to estimate biomass using both egg and acoustic surveys. We used model fitting to determine the error structure that best fit the data. The states estimated by that best fit model are estimates of biomass based on egg and acoustic surveys.


The derivation of the EM algorithm for  with linear constraints is here: http://cran.r-project.org/web/packages/MARSS/vignettes/EMDerivation.pdf

The User Guide/vignette for the MARSS package is available here: http://cran.r-project.org/web/packages/MARSS/vignettes/UserGuide.pdf
