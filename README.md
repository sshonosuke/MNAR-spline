# MNAR-spline (Bayesian semiparametric estimaton with nonignorable nonresponse)

This package implements Bayesian semiparametric estimaton for response models, as proposed by the following papers.

Sugasawa, S., Morikawa, K. and Takahata, K. (2019). Bayesian semiparametric estimaton with nonignorable nonresponse. arXiv

Functions are implemented in BSS-function.R available in the repository.

```{r}
source("BSS-function.R")   # require "MCMCpack" and "pgdraw" packages
```

Load demo dataset
```{r}
load("data.RData")
```

Fit the semiparametric reponse model

Input of `MNAR`

- `oby`: response vector (missing part is `NA`)
- `X`: matrix of covariates (the first column should be 1's) in the outcome model
- `Z`: matrix of covariates in the response model
- `S`: vector of missing indicator (1 for observed; 0 for missing )
- `q`: order of spline 
- `K`: number of knots
- `mc`: length of MCMC 
- `burn`: burn-in period
- `Knot`: locations of knots are adaptively esimated if `T`

Output of `MNAR`: List object of MCMC results

- `Beta`: regression coeffieicnts in the outcome model
- `Sigma`: squared value of error variance in the outcome model
- `Phi`: coefficients of the polynomial term of `Y` in the response model
- `Gamma`: coefficients of the spline term of `Y` in the response model
- `Delta`: coefficients for the covaraites in the response model
- `Lam`: precision parameter in the prior for `Gamma`
- `a`: scale parameter controlling locations of knots

```{r}
set.seed(1)
qq=2
K=10
mc=7000
bn=2000
fit=BSS.LM(Y,X,Z,S,q=qq,K=K,mc=mc,burn=bn)
```

Posterior mean of regression coeffieicnts
```{r}
apply(fit$Beta,2,mean)
```



