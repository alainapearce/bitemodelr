# bitemodelr

In its current form, the bitemodelr package is written to support a
simulation study that aimed determine if observational coded bite timing
and average bite size could be used to model cumulative intake curves.
Additionally, the simulation study compared the Logistic Ordinary
Differential Equation (LODE) model (Thomas et al., 2017) to Kissielf’s
Quadratic model. For more details about the simulation study see:
<https://github.com/alainapearce/LODEModel_SimStudy>

Currently, this package contains functions to simulate cumulative intake
datasets with process and/or measurement noise added. It is also able to
fit parameter values for both the LODE and Quadratic models and provide
confidence intervals for each parameter. Additionally, it includes
functions to evaluate the performance of model fit including
distinctness and error (goodness-of-fit, RMSE, psuedo-*R*<sup>2</sup>).
This package is still being updated and refined. Further documentation
is also planned for the future.

## Prerequisites

Included dependencies: stats, truncnorm, MASS, faux, methods, utils,
truncdist, HelpersMG

## Installation

library(devtools) devtools::install\_github(“alainapearce/bitemodelr”)

## Citation

[![DOI](https://zenodo.org/badge/265935094.svg)](https://zenodo.org/badge/latestdoi/265935094)

## License

GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007

## Contact

Alaina Pearce (<azp271@psu.edu>)
