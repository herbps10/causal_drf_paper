# Causal-DRF: Conditional Kernel Treatment Effect Estimation using Distributional Random Forest

This repository contains the code to replicate the results presented in the paper:  
**[Causal-DRF: Conditional Kernel Treatment Effect Estimation using Distributional Random Forest](https://arxiv.org/abs/2411.08778)**.

## Requirements

To run this code, you will need to install a new version of the `drf` R package from the `causal-clean` branch on GitHub.

## Installation

Follow these steps to install the required version of the `drf` package:

The development version can be installed from github

```R
devtools::install_github("herbps10/drf/causal-clean",subdir = "r-package/drf")
```

For this one may clone the repo, and then within the r-package folder run

```R
Rscript build_package.R
```
