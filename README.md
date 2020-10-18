# PlasmoMAPI
### Version 0.1.0
[![Travis build status](https://travis-ci.org/mrc-ide/PlasmoMAPI.svg?branch=master)](https://travis-ci.org/mrc-ide/PlasmoMAPI)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/mrc-ide/PlasmoMAPI?branch=master&svg=true)](https://ci.appveyor.com/project/mrc-ide/PlasmoMAPI)

[![R build status](https://github.com/mrc-ide/PlasmoMAPI/workflows/R-CMD-check/badge.svg)](https://github.com/mrc-ide/PlasmoMAPI/actions)

--------------------------------------------------------------------------------------------------------------------------------

Genetic data from Plasmodium parasites can tell us something about the spatial connectivity of parasite populations. Infections that trace back to a common ancestor in the recent past tend to be highly related, meaning we can use estimates of pairwise relatedness between sampling locations to infer levels of gene flow in the population, and to identify barriers and corridors of migration. PlasmoMAPI does this in a statistically robust way, accounting for the natural fall-off of relatedness with distance to identify significant discontinuities in the overall spatial pattern.

A complete description of the method, including worked example, can be found on the [PlasmoMAPI webpage](https://mrc-ide.github.io/PlasmoMAPI/).

