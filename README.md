# NiNBayes

## About

The purpose of this package is to implemetn *N*on-*i*gnorable *N*onparametric *Bayes* estimation procedures. Emphasis is given to conducting methods which allow the user to conduct a sensitivity analysis by using different assumptions. We allow for both fully-Bayes estimation of aspects of the full-data distribution and (more usefully?) multiple imputation of data under different assumptions about the missing data. 

## Installation

```
library(devtools)
install_github("theodds/NiNBayes")
```

## Using the package

This package is currently in a very preliminary status. Right now, the package allows for analyzing simple multivariate binary data, though there are plans in the future to (1) allow for the inclusion of covariates and (2) allow for continuous responses. 
