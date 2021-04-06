# ICScure #
ICScure (which stands for <ins>**I**</ins>nformative <ins>**C**</ins>uster <ins>**S**</ins>ize <ins>**cure**</ins>) is a package that performs estimation and inference for clustered current status data with informative cluster size using the method proposed by Lam et al. (2021) <DOI: [10.1002/sim.891](https://doi.org/10.1002/sim.8910)>.

**ICScure** relies on the R-packages `stats`, `stabledist`, `numDeriv`, `MASS`, which is hosted on CRAN.

# Installation #
**ICScure** can be installed from github directly:
```
install.packages("devtools")
library(devtools)
install_github("alexwky/ICScure")
```

# Usage #
The package contains 3 functions:
Functions  | Description
------------- | -------------
BernsteinPolynomial  | Calculate the values of a Bernstein polynomial at given time points.
ICDASim  | Generate a data set according to the simulation studies in Lam et al. (2021) <DOI: [10.1002/sim.891](https://doi.org/10.1002/sim.8910)>
Est.ICScure  |  Perform the cluster-weighted GEE or GEE estimation of Lam et al. (2021) <DOI: [10.1002/sim.891](https://doi.org/10.1002/sim.8910)>


<ins>**BernsteinPolynomial**</ins>

```
BernsteinPolynomial(t, psi, mint, maxt)
```
This function evaluate the Bernstein polynomial at given time points `t` using user-provided coefficients `psi` and support (`mint`,`maxt`).

Example:
```
t <- seq(1,10,1)
psi <- seq(0.25,1,0.25)
BernsteinPolynomial(t = t, psi = psi, mint = 0, maxt = 5)
# 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0
```

<ins>**ICDASim**</ins>

```
ICDASim(seed = NA, n, beta00, beta10, beta20, rho, gamma, cs)
```
This function generates a data set according to the simulation studies in Lam et al. (2021) <DOI: 10.1002/sim.8910>.

Example:
```
data <- ICDASim(seed = 1942, n = 100, beta00 = 0.5, beta10 = -1, beta20 = 1, cs = 40, rho = 0.5, gamma = 0.5)
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
            weighting = TRUE, t.min = 0, t.max = NA, reltol = 10^-6, maxit = 1000, calSE = TRUE)
```
This function performs the cluster-weighted GEE or GEE estimation of Lam et al. (2021) <DOI: 10.1002/sim.8910>

`data` is a `n x (p+3)` matrix, where `n` is the sample size and `p` is the number of covariates. The first column consists of cluster indices, the second column consists of the observation time, the third column consists of the event indicator, and the fourth to the last columns consist of the covariates (not including the intercept). The set of covariates can be empty. The format of `data` is as follow:

**Cluster Index**  | **Observation Time**  | **Event Indicator** | **1<sup>st</sup> covariate** | **2<sup>nd</sup> covariate** | ... | **p<sup>th</sup> covariate**
------------- | ------------- | ------------- | ------------- | ------------- | ------------- | -------------
1  | 3.7322 | 1 | 1 | 0.0888 | ... | 1
1  | 4.0000 | 1 | 0 | -0.4965 | ... | 0



Example:
```
Dataset <- ICDASim(seed = 1942, n = 100, beta00 = 0.5, beta10 = -1, beta20 = 1, cs = 40, rho = 0.5, gamma = 0.5)
Result <- Est.ICScure(data = Dataset, rho = 0.5, degree = 3, weighting = TRUE)
Result

# $degree
# [1] 3
#
# $psi
# [1] 0.9917128 0.9955245 1.0000000
#
# $beta
# [1]  0.7091996 -1.0196841  0.9260482
#
# $betaSE
# [1] 0.13000664 0.11749335 0.07599541
#
# $iteration
# [1] 44
#
# $covergence
# [1] "TRUE"
```



# Help #

Details about the package can be found in the user manual:
```
?ICScure
```

# Contact #
Wong Kin Yau, Alex <<kin-yau.wong@polyu.edu.hk>>

# Reference #
Lam, K. F., Lee, C. Y., Wong, K. Y., & Bandyopadhyay, D. (2021). Marginal analysis of current status data with informative cluster size using a class of semiparametric transformation cure models. Statistics in Medicine [online], DOI:10.1002/sim.8910

