# ICScure #
ICScure stand for <ins>**I**</ins>nformative <ins>**C**</ins>uster <ins>**S**</ins>ize <ins>**cure**</ins> model. ICScure is a package to perform estimation and inference for clustered current status data with informative cluster size using the method proposed by Lam et al. (2021) <DOI: [10.1002/sim.891](https://doi.org/10.1002/sim.8910)>. It is implemented in the R-package **ICScure**.

**ICScure** relies on the R-package `stats`, `stabledist`, `numDeriv`, `MASS`, which is hosted on CRAN.

# Installation #
**IBoost** can be installed from github directly:
```
install.packages("devtools")
library(devtools)
install_github("alexwky/ICScure")
```

# Usage #
The package contains 3 functions:
Functions  | Description
------------- | -------------
BernsteinPolynomial  | Calculate the values of the Bernstein polynomial at time `t`.
ICDASim  | Generate a data set according to the simulation studies in Lam et al. (2021) <DOI: [10.1002/sim.891](https://doi.org/10.1002/sim.8910)>
Est.ICScure  |  Perform the cluster-weighted GEE or GEE estimation of Lam et al. (2021) <DOI: [10.1002/sim.891](https://doi.org/10.1002/sim.8910)>


<ins>**BernsteinPolynomial**</ins>

```
BernsteinPolynomial(t, psi, mint, maxt)
```
This function evaluate the Bernstein polynomial at given time points `t` using user-provided coefficients `psi`, degree, and support `mint` `maxt`.

Example:
```
t <- seq(1,10,1)
psi <- seq(0.25,1,0.25)
BernsteinPolynomial(t = t, psi = psi, mint = 0, maxt = 5)
# 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0
```

<ins>**ICDASim**</ins>

```
ICDASim(seed = NA, n, beta00, beta10, beta20, rho, PoAlpha, cs, zeta)
```
This function generate a data set according to the simulation studies in Lam et al. (2021) <DOI: 10.1002/sim.8910>. From the simulator, only 2 covariates will be generated. In order to reproduce the same dataset, seed can be input optionally. 

Example:
```
data<-ICDASim(seed = 1942, n = 100, beta00 = 0.5, beta10 = -1, beta20 = 1, cs = 40, rho=0.5,PoAlpha = 0.5, zeta=0.75)
head(data)

#   family        Ci delta x1         x2
# 1      1 4.0000000     0  0  0.5042357
# 2      1 2.2748313     1  0  1.8434784
# 3      1 4.0000000     0  1  0.2093028
# 4      1 0.4409429     0  0 -0.1799730
# 5      1 4.0000000     1  0  0.5579359
# 6      1 1.4709212     0  0  0.4700130
```

<ins>**Est.ICScure**</ins>
  
```
Est.ICScure(data, rho, degree,
            weighting = TRUE, t.min = 0, t.max = NA, reltol = 10^-6, maxit = 1000, Hessian = TRUE)
```
This function perform the cluster-weighted GEE or GEE estimation of Lam et al. (2021) <DOI: 10.1002/sim.8910>
`data` is a n x (p+3) matrix, where `n` is the sample size and `p` is the numbers of covariates. Since intercept does not include in the covariates, `data` with no covariates is also allowed. The format of `data` is as follow:
**Cluster Number**  | **Observation Time**  | **Event Indicator** | **1<sup>st</sup> covariates** | **2<sup>nd</sup> covariates** | ... | **n<sup>th</sup> covariates**
------------- | ------------- | ------------- | ------------- | ------------- | ------------- | -------------
1  | 3.7322 | 1 | 1 | 0.0888 | ... | 1
1  | 4.0000 | 1 | 0 | -0.4965 | ... | 0



Example:
```
Dataset <- ICDASim(seed = 1942, n = 100, beta00 = 0.5, beta10 = -1, beta20 = 1, cs = 40, rho=0.5,PoAlpha = 0.5, zeta=0.75)
Result <- Est.ICScure(data=Dataset, rho=0.5, degree=3, weighting=FALSE)
Result

# $degree
# [1] 3
# 
# $psi
# [1] 0.9758965 0.9772273 1.0000000
# 
# $beta
# [1]  0.2244137 -1.0920922  1.0339562
# 
# $betaSE
# [1] 0.13324418 0.11512944 0.08043015
# 
# $iteration
# [1] 91
# 
# $covergence
# [1] "TRUE"
# 
# $AIC
# [1] 1793.856
# 
# $ct
# [1] 4

```



# Help #

**Est.ICScure()** is the main function that implements the I-Boost procedure. Details about the function can be found in the user manual:
```
?ICScure
```

# Contact #
Dr. Wong Kin Yau, Alex <<kin-yau.wong@polyu.edu.hk>>

# Reference #
Lam, K. F., Lee, C. Y., Wong, K. Y., & Bandyopadhyay, D. (2021). Marginal analysis of current status data with informative cluster size using a class of semiparametric transformation cure models. Statistics in Medicine, 1–13. https://doi.org/10.1002/sim.8910
