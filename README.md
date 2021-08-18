# Hi-LMM -Early Access

## 1. Getting started
####	Downloading Hi-LMM
HiRRM can be downloaded https://github.com/RunKingProgram/Hi-LMM. It can be installed as a regular R package.
####	Installing Hi-LMM
Hi-GLMM links to R packages Rcpp, RcppArmadillo, RcppEigen, snow, parallel,data.table, nlme and BEDMatrix. These dependencies should be installed before installing Hi-GLMM. In addition, **Hi-LMM requires GEMMA software (chmod 777) under your working directory**. Here is an example for installing Hi-GLMM and all its dependencies in an R session(assuming none of the R packages other than the default has been installed):
```
install.packages( c( "Rcpp", "RcppArmadillo", "RcppEigen", "snow", "parallel", "data.table", "BEDMatrix" ))
system( “R CMD install HiGLMM_v0.9_MacOS.tgz” )
```
## 2. Main functions
The current version of HiRRM includes two main functions:
```
gdata = Data_HiLMM(filename, h2 ,ebv = NULL) 
HiLMM(gdata,Test=c("Jiont","Separate") ,thread=NULL,QQ=F,Manh=F)
```
#### Arguments
#### filename
An object class of character，which consists of three PLINK BED files with the same name. For example, Genotype.bed, Genotype.bim and Genotype.fam.

#### h2
An object class of numeric. Heritability will be used for estimating breeding values. Must provide, if ebv = NULL.
#### ebv
Estimated breeding values. NULL by default.
#### gdata
An object class of list generated from the Data_HiGLMM function.
####Test

An optional for association test, a test at once or joint analysis. You can choose “Separate” for a test at once and “Joint” for further joint analysis based on the result of “Separate”.

#### QQ

logicals. If TURE, Q-Q plot would be drawn.

#### Manh

logicals. If TURE, Manh plot would be drawn.


## 3.Example
```
library(HiLMM)
library(parallel)
library(BEDMatrix)
library(data.table)
library(snow)
library(RcppArmadillo)


setwd("./example")
gdata = Data_HiLMM("geno", 0.656) 
HiLMM(gdata,QQ=T,Manh=T)
```


