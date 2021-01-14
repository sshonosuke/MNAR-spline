##  preparation 
rm(list=ls())
set.seed(1)


##  load function
source("MNAR-LM-function.R")
quant <- function(x){ quantile(x,prob=c(0.025,0.975)) }

## simualted data
n <- 500

Beta <- c(0.8, 0.8, -0.5)
sig <- 1
rho <- 0.2

mat <- matrix(rho, 2, 2)
diag(mat) <- 1
X <- cbind(1, mvrnorm(n, rep(0, 2), mat))   # sample X
Y <- as.vector(X%*%Beta) + sig*rnorm(n)


# response mechanism 
sc <- 3    # from 1 to 7
if(sc==1){ U=1.5-0.5*Y+0.2*X[,2] }
if(sc==2){ U=2.5-0.2*Y-0.4*Y^2+0.2*X[,2] }
if(sc==3){ U=1.5-0.2*Y-0.4*Y^2+0.2*X[,2] }
if(sc==4){ U=0.7*Y^2+0.2*X[,2] }
if(sc==5){ U=0.5*Y^2+X[,2]^2 }
if(sc==6){ U=1.5-2*sin(Y)+X[,2]^2 }
if(sc==7){ U=0.7*Y^2+0.2*X[,3] }
prob <- 1/(1+exp(-U))
if(sc==3){ prob=1-exp(-exp(U)) }

#  missing indicator
S <- rbinom(n,1,prob)
oby <- Y
oby[S==0]<- NA
Z <- as.matrix(X[,2])


# Linear selection 
fit.LS <- LR.LM(oby, X, Z, S, mc=5000, burn=2000)
apply(fit.LS$Beta, 2, mean)     # posterior mean
t(apply(fit.LS$Beta, 2, quant))     # posterior credible intervals 
fit.LS$Aux$DIC     # DIC

# Spline selection (semipara) 
fit.SR <- SR.LM(oby, X, Z, S, q=2, K=10, mc=5000, burn=2000, NP=F)
apply(fit.SR$Beta, 2, mean)   # posterior mean
t(apply(fit.SR$Beta, 2, quant))   # posterior credible intervals
fit.SR$Aux$DIC    # DIC


# Spline selection (fully nonpara) 
fit.NR <- SR.LM(oby, X, Z, S, q=2, K=10, mc=5000, burn=2000, NP=T)
apply(fit.NR$Beta, 2, mean)   # posterior mean
t(apply(fit.NR$Beta, 2, quant))    # posterior credible intervals
fit.NR$Aux$DIC    # DIC

