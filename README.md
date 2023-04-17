# MediationFDR

Mediator Selection Controlling For FDR

## About 

This package contains the function to select the mediators in the high-dimensional mediators setting with FDR control, by utilizing the recent development of the knockoff method for FDR-controlled variable selection. 


## Installation

```{r}
## install package
# install.packages("devtools")
devtools::install_github("ruiyangli1/MediationFDR")

## load package
library(MediationFDR)
```


## Usage

```{r}
MediationFDR(X, Y, M)
```
For more example, please see [here](https://ruiyangli1.github.io/MediationFDR/articles/Example.html).


