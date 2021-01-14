##  preparation 
rm(list=ls())
set.seed(1)

## load dataset 
load("Schizo_PANSS.rda")

## load functions
source("MNAR-LMM-function.R")
set.seed(1)


## data 
data <- Schizo_PANSS
n <- dim(data)[1]
TT <- 5

Y <- data[,4:8]
obY <- as.numeric(as.vector(t(Y)))
S <- ifelse(is.na(obY), 0, 1)
mean(S)      # overall response probability 
tr <- data$Treat    # treatment indicator
tr[tr==-1] <- 0
time <- c(1, 2, 4, 6, 8)    # observed time 


# plot
matplot(time, t(Y), type="l", col=tr+1, lty=1, ylim=c(-100,100))
legend("top", c("Control","Treatment"), col=1:2, lty=1, ncol=2)


# treatment indicator 
Tr <- matrix(rep(data$Treat, TT), n, TT)
TTr <- as.vector(t(Tr))
tr <- ifelse(TTr==1, 1, 0)
time <- rep(c(1,2,4,6,8), n)

# response indicator at previous time 
SS <- rbind(1, matrix(S, TT, n)[1:(TT-1),])
SS <- as.vector(SS)

# covariate matrix
X <- cbind(1, time, time^2, time^3, tr, tr*time, tr*time^2, tr*time^3)
Z <- cbind(time, SS, tr)
Z2 <- cbind(time, time^2, SS, tr)
ID <- sort(rep(1:n, TT))



## method 
mc <- 5000
bn <- 2000

# linear response model 
LR.fit <- LR.LMM(obY, X, Z, S, ID, mc=mc, burn=bn)
LR.fit$DIC

# semiparametric response model (it may take a long time)
SR.fit <- SR.LMM(obY, X, Z, S, ID, mc=mc, burn=bn, q=2, K=10)
SR.fit$DIC


##  result 
CI <- function(x){ quantile(x,prob=c(0.025,0.975)) }
tt <- c(1, 2, 4, 6, 8)
XX <- cbind(1, tt, tt^2, tt^3)




## plot (treatment effect)
par(mfcol=c(1,2))
ran <- c(-25, -5)
# LR
pos1 <- XX%*%t(LR.fit$Beta[,1:4])
pm1 <- apply(pos1, 1, mean)
CI1 <- apply(pos1, 1, CI)
pos2 <- pos1+XX%*%t(LR.fit$Beta[,5:8])
pm2 <- apply(pos2, 1, mean)
CI2 <- apply(pos2, 1, CI)
matplot(tt, cbind(pm1,pm2), type="l", col=c(1,2), lty=1, ylim=ran, main="Linear selection", ylab="PANSS score", xlab="Time (week)")
points(tt, pm1, col=1, pch=8)
points(tt, pm2, col=2, pch=8)
polygon(c(tt,rev(tt)), c(CI1[1,],rev(CI1[2,])), col="#30303020", border=NA)
polygon(c(tt,rev(tt)), c(CI2[1,],rev(CI2[2,])), col="#ff000020", border=NA)
legend("topright", c("control","treatment"), col=c(1,2), lty=1)
# SR
fit <- SR.fit
pos1 <- XX%*%t(SR.fit$Beta[,1:4])
pm1 <- apply(pos1, 1, mean)
CI1 <- apply(pos1, 1, CI)
pos2 <- pos1+XX%*%t(SR.fit$Beta[,5:8])
pm2 <- apply(pos2, 1, mean)
CI2 <- apply(pos2, 1, CI)
matplot(tt, cbind(pm1,pm2), type="l", col=c(1,2), lty=1, ylim=ran, main="Semiparametric selection", ylab="PANSS score", xlab="Time (week)")
points(tt, pm1, col=1, pch=8)
points(tt, pm2, col=2, pch=8)
polygon(c(tt,rev(tt)), c(CI1[1,],rev(CI1[2,])), col="#30303020", border=NA)
polygon(c(tt,rev(tt)), c(CI2[1,],rev(CI2[2,])), col="#ff000020", border=NA)
legend("topright", c("control","treatment"), col=c(1,2), lty=1)



##  estimated response model
SS <- 0     #  0 or 1
aa <- mean(SR.fit$a)
K <- 10
kappa <- quantile(na.omit(obY), prob=c(0.05,0.95))
kappa1 <- kappa[1]-aa*diff(kappa)/2
kappa2 <- kappa[2]+aa*diff(kappa)/2
knots <- seq(kappa1,kappa2,length=K)


mc <- mc - bn
L <- 100
eval.time <- 4
if(SS==0){ yy <- seq(-30, 20, length=L) }
if(SS==1){ yy <- seq(-40, 20, length=L) }

add1c <- as.vector(LR.fit$Delta%*%c(eval.time,SS,0))
add2c <- as.vector(SR.fit$Delta%*%c(eval.time,SS,0))
add1t <- as.vector(LR.fit$Delta%*%c(eval.time,SS,1))
add2t <- as.vector(SR.fit$Delta%*%c(eval.time,SS,1))
Mis1c <- matrix(NA, mc, L)
Mis1t <- matrix(NA, mc, L)
Mis2c <- matrix(NA, mc, L)
Mis2t <- matrix(NA, mc, L)
for(k in 1:L){
  Mis1c[,k] <- as.vector(LR.fit$Phi%*%c(1,yy[k]))+add1c
  Mis1t[,k] <- as.vector(LR.fit$Phi%*%c(1,yy[k]))+add1t
  SP <- ( (yy[k]-knots)*ifelse(yy[k]-knots>0,1,0) )^2  
  Mis2c[,k] <- as.vector(SR.fit$Phi%*%c(1,yy[k],yy[k]^2))+as.vector(SR.fit$Gam%*%SP)+add2c
  Mis2t[,k] <- as.vector(SR.fit$Phi%*%c(1,yy[k],yy[k]^2))+as.vector(SR.fit$Gam%*%SP)+add2t
}

logit <- function(x){ 1/(1+exp(-x)) }
ran <- cbind(c(0,0.3), c(0.5,1))



## plot (response model) 
par(mfcol=c(1,2))

# control group
est.Mis <- cbind(apply(logit(Mis1c),2,mean), apply(logit(Mis2c),2,mean))
CI1 <- apply(logit(Mis1c), 2, CI)
CI2 <- apply(logit(Mis2c), 2, CI)
matplot(yy, est.Mis, type="l", lty=1, main=paste0("Control group (S=",SS,")"), ylab="Probability", xlab="PANSS score", ylim=ran[,SS+1])
polygon(c(yy,rev(yy)), c(CI1[1,],rev(CI1[2,])), col="#30303020", border=NA)
polygon(c(yy,rev(yy)), c(CI2[1,],rev(CI2[2,])), col="#ff000020", border=NA)
if(SS==0){ legend("topright",c("Linear selection","Semiparametric selection"),col=1:2,lty=1) }
if(SS==1){ legend("bottomleft",c("Linear selection","Semiparametric selection"),col=1:2,lty=1) }

# treatment group
est.Mis <- cbind(apply(logit(Mis1t),2,mean), apply(logit(Mis2t),2,mean))
CI1 <- apply(logit(Mis1t), 2, CI)
CI2 <- apply(logit(Mis2t), 2, CI)
matplot(yy, est.Mis, type="l", lty=1, main=paste0("Treatment group (S=",SS,")"), ylab="Probability", xlab="PANSS score", ylim=ran[,SS+1])
polygon(c(yy,rev(yy)), c(CI1[1,],rev(CI1[2,])), col="#30303020", border=NA)
polygon(c(yy,rev(yy)), c(CI2[1,],rev(CI2[2,])), col="#ff000020", border=NA)
if(SS==0){ legend("topright",c("Linear selection","Semiparametric selection"),col=1:2,lty=1) }
if(SS==1){ legend("bottomleft",c("Linear selection","Semiparametric selection"),col=1:2,lty=1) }

