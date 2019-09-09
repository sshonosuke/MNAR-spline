# MNAR-spline (Bayesian semiparametric estimaton with nonignorable nonresponse)

This package implements Bayesian semiparametric estimaton for response models, as proposed by the following papers.

Sugasawa, S., Morikawa, K. and Takahata, K. (2019). Bayesian semiparametric estimaton with nonignorable nonresponse. https://arxiv.org/abs/1909.02878

Functions are implemented in BSS-function.R available in the repository.

```{r}
source("BSS-function.R")   # require "MCMCpack" and "pgdraw" packages
```


# Example: linear regression 

Load demo dataset
```{r}
load("data.RData")
```

Fit the semiparametric reponse model with linear regression for the outcome model

Input of `BSS.LM`

- `oby`: response vector (missing part is `NA`)
- `X`: matrix of covariates (the first column should be 1's) in the outcome model
- `Z`: matrix of covariates in the response model
- `S`: vector of missing indicator (1 for observed; 0 for missing )
- `q`: order of spline 
- `K`: number of knots
- `mc`: length of MCMC 
- `burn`: burn-in period
- `Knot`: locations of knots are adaptively esimated if `T`

Output of `BSS.LM`: List object of MCMC results

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

Estimated response model and 95% point-wise credible intervals
```{r}
aa=mean(fit$a)  
kappa=quantile(na.omit(Y),prob=c(0.05,0.95))
kappa1=kappa[1]-aa*diff(kappa)/2
kappa2=kappa[2]+aa*diff(kappa)/2
knots=seq(kappa1,kappa2,length=K)

nn=mc-bn      # number of posterior samples
L=100     # number of evalaution points for response values
ran=range(na.omit(Y))
yy=seq(ran[1],ran[2],length=L)
zz=1     # evaluation point of covariate
add=as.vector(fit$Delta*zz)
logit=function(x){ 1/(1+exp(-x)) }

Mis=matrix(NA,nn,L)
for(k in 1:L){
  YY=yy[k]^(0:qq)
  SP=( (yy[k]-knots)*ifelse(yy[k]-knots>0,1,0) )^qq  
  Mis[,k]=logit( as.vector(fit$Phi%*%YY)+as.vector(fit$Gamma%*%SP)+add )
}

quant=function(x){ quantile(x,prob=c(0.025,0.975)) }
CI=apply(Mis,2,quant)     # point-wise credible intervals
plot(yy,apply(Mis,2,mean),type="l",lty=1,ylab="Probability",xlab="Response value",main="Response model",ylim=range(CI))
polygon(c(yy,rev(yy)),c(CI[1,],rev(CI[2,])),col="#30303020",border=NA)
```


# Example: linear mixed model (random intercept)

Load PANSS dataset (downloaded from https://cran.r-project.org/web/packages/Surrogate/index.html)
```{r}
load("Schizo_PANSS.rda")
data=Schizo_PANSS
n=dim(data)[1]
TT=5     # number of measurements

Y=data[,4:8]
obY=as.numeric(as.vector(t(Y)))    # response variable
S=ifelse(is.na(obY),0,1)     # missing indicator
mean(S)

Tr=as.vector( t(matrix(rep(data$Treat,TT),n,TT)) )
tr=ifelse(Tr==1,1,0)     # treatment indicator (1 for treatment and 0 for placebo)
time=rep(c(1,2,4,6,8),n)     # time of measurement
SS=as.vector( rbind(1,matrix(S,TT,n)[1:(TT-1),]) )      # one time before missing indicator 
X=cbind(1,time,time^2,time^3,tr,tr*time,tr*time^2,tr*time^3)      # matrix of covariates in the outcome model
Z=cbind(time,SS,tr)     # matrix of covariates in the response model
ID=sort(rep(1:n,TT))     # patients ID
```

Fit the semiparametric reponse model with linear mixed model for the outcome model

Input of `BSS.LMM`

- `oby`: response vector (missing part is `NA`)
- `X`: matrix of covariates (the first column should be 1's) in the outcome model
- `Z`: matrix of covariates in the response model
- `S`: vector of missing indicator (1 for observed; 0 for missing )
- `ID`: vector of cluster ID
- `q`: order of spline 
- `K`: number of knots
- `mc`: length of MCMC 
- `burn`: burn-in period
- `Knot`: locations of knots are adaptively esimated if `T`

Output of `BSS.LM`: List object of MCMC results

- `Beta`: regression coeffieicnts in the outcome model
- `Sigma`: squared value of error variance in the outcome model
- `RE`: cluster-specific random effects (intercepts) in the outcome model
- `Tau`: squared value of random effects variance in the outcome model
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
fit=BSS.LMM(obY,X,Z,S,ID,mc=mc,burn=bn,q=qq,K=K)
```

Time varying PANSS scores for two groups
```{r}
tt=c(1,2,4,6,8)
XX=cbind(1,tt,tt^2,tt^3)

pos1=XX%*%t(fit$Beta[,1:4])
pm1=apply(pos1,1,mean)
CI1=apply(pos1,1,quant)
pos2=pos1+XX%*%t(fit$Beta[,5:8])
pm2=apply(pos2,1,mean)
CI2=apply(pos2,1,quant)

ran=range(CI1,CI2)
matplot(tt,cbind(pm1,pm2),type="l",col=c(1,2),lty=1,ylim=ran,main="Semiparametric selection",ylab="PANSS score",xlab="Time (week)")
points(tt,pm1,col=1,pch=8)
points(tt,pm2,col=2,pch=8)
polygon(c(tt,rev(tt)),c(CI1[1,],rev(CI1[2,])),col="#30303020",border=NA)
polygon(c(tt,rev(tt)),c(CI2[1,],rev(CI2[2,])),col="#ff000020",border=NA)
legend("topright",c("control","treatment"),col=c(1,2),lty=1)
```

Estimated response model and 95% point-wise credible intervals
```{r}
aa=mean(fit$a)
kappa=quantile(na.omit(obY),prob=c(0.05,0.95))
kappa1=kappa[1]-aa*diff(kappa)/2
kappa2=kappa[2]+aa*diff(kappa)/2
knots=seq(kappa1,kappa2,length=K)

nn=mc-bn
L=100
eval.time=4
SS=1     # set 1 or 0 
if(SS==0){ yy=seq(-30,20,length=L) }
if(SS==1){ yy=seq(-40,20,length=L) }
add.c=as.vector(fit$Delta%*%c(eval.time,SS,0))
add.t=as.vector(fit$Delta%*%c(eval.time,SS,1))

Mis.c=matrix(NA,nn,L); Mis.t=matrix(NA,nn,L)
for(k in 1:L){
  YY=yy[k]^(0:qq)
  SP=( (yy[k]-knots)*ifelse(yy[k]-knots>0,1,0) )^qq  
  Mis.c[,k]=logit( as.vector(fit$Phi%*%YY)+as.vector(fit$Gamma%*%SP)+add.c )
  Mis.t[,k]=logit( as.vector(fit$Phi%*%YY)+as.vector(fit$Gamma%*%SP)+add.t )
}

ran=range(CI1,CI2)
est.Mis=cbind(apply(Mis.c,2,mean),apply(Mis.t,2,mean))
CI1=apply(Mis.c,2,quant)
CI2=apply(Mis.t,2,quant)
matplot(yy,est.Mis,type="l",lty=1,main=paste0("Response model (S=",SS,")"),ylab="Probability",xlab="PANSS score",ylim=ran)
polygon(c(yy,rev(yy)),c(CI1[1,],rev(CI1[2,])),col="#30303020",border=NA)
polygon(c(yy,rev(yy)),c(CI2[1,],rev(CI2[2,])),col="#ff000020",border=NA)
legend("bottomleft",c("Control","Treatment"),col=1:2,lty=1) 
```



