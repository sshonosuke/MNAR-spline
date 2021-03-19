library(MCMCpack)
library(pgdraw)
library(mvtnorm)





###  Bayesian semiparametric modeling  ###
# outcome model: linear regression 
SR.LM <- function(obY, X, Z, S, q=2, K=10, mc=5000, burn=2000, Knot=T, NP=T){
  n <- length(obY)
  p <- dim(X)[2]
  rr <- dim(Z)[2]
  kappa <- quantile(na.omit(obY), prob=c(0.05,0.95))
  IX <- solve(t(X)%*%X)
  
  Set.Knots <- function(a){
    kappa1 <- kappa[1]-a*diff(kappa)/2
    kappa2 <- kappa[2]+a*diff(kappa)/2
    seq(kappa1, kappa2, length=K)
  }
  
  logistic <- function(x){ 1/(1+exp(-x)) }
  
  # Knots (initial)
  knots <- Set.Knots(0.2) 
  
  ccn <- 0.0001   # normal prior
  ccg <- 1   # gamma prior
  
  Y.pos <- matrix(NA, mc, n)
  Beta.pos <- matrix(NA, mc, p)
  Sig.pos <- c()
  Phi.pos <- matrix(NA, mc, q+1)
  Gam.pos <- matrix(NA, mc, K)
  Delta.pos <- matrix(NA, mc, rr)
  Lam.pos <- c()
  Om.pos <- matrix(NA, mc, n)
  a.pos <- c()
  Xi.pos <- matrix(NA, mc, K)
  Pi.pos <- matrix(NA, mc, n)
  
  # Initial values
  YY <- c()
  YY[S==0] <- 0
  YY[S==1] <- obY[S==1]
  Beta <- rep(0, p)
  Sig <- 1
  Phi <- rep(0.1, q+1)
  Gam <- rep(0.01, K)
  Delta <- rep(0, rr)
  Lam <- 1
  Om <- rep(1,n)
  aa <- 0.2
  
  # nonparametric term
  if(NP){
    Xi <- rep(0, K)
    Lam2 <- 1
    fit <- kmeans(Z, K)
    KAP <- fit$centers
    ad <- 1.5    # adjustment constant for scale
    Scale <- ad*sqrt( fit$withinss/fit$size )
    ZZ <- matrix(NA, n, K)
    for(k in 1:K){
      ss <- apply((Z-KAP[k,])^2, 1, sum)
      ZZ[,k] <- dmvnorm(Z, KAP[k,], Scale[k]^2*diag(rr))
    }
  }
  RB <- 0
  
  # functions for splines
  SPterm <- function(x, knots){
    nn <- length(x)
    mat1 <- matrix(NA, nn, q+1)
    for(j in 0:q){ mat1[,j+1] <- x^j }
    mat2 <- matrix(NA,nn,K)
    for(l in 1:K){
      val <- x-knots[l]
      mat2[,l] <- (val*ifelse(val>0, 1, 0))^q
    }
    return(cbind(mat1, mat2))
  }
  
  dSPterm <- function(x, knots){
    nn <- length(x)
    mat1 <- matrix(NA, nn, q+1)
    mat1[,1] <- 0
    for(j in 1:q){ mat1[,j+1] <- j*x^(j-1) }
    mat2 <- matrix(NA,nn,K)
    for(l in 1:K){
      val <- x-knots[l]
      mat2[,l] <- q*(val*ifelse(val>0,1,0))^(q-1)
    }
    return(cbind(mat1, mat2))
  }
  
  ## MCMC
  for(k in 1:mc){
    if(NP){ RB <- as.vector(ZZ%*%Xi) }
    Zreg <- as.vector(Z%*%Delta)
    # Omega
    U <- as.vector( SPterm(YY,knots)%*%c(Phi,Gam) + Zreg + RB )
    Om <- pgdraw(rep(1,n), U)
    Om.pos[k,] <- Om
    # Y  (Langevin MC)
    h <- 0.2   # step size 
    mu <- as.vector(X%*%Beta)
    V <- 0.5*(YY-mu)^2/Sig^2 + 0.5*U + 0.5*U^2*Om
    dV <- (YY-mu)/Sig^2 + (0.5+U*Om)*as.vector( dSPterm(YY,knots)%*%c(Phi,Gam) )
    YY.prop <- rnorm(n, YY-h*dV, sqrt(2*h))
    U.prop <- as.vector( SPterm(YY.prop,knots)%*%c(Phi,Gam) + Zreg + RB )
    V.prop <- 0.5*(YY.prop-mu)^2/Sig^2 + 0.5*U.prop + 0.5*U.prop^2*Om
    dV.prop <- (YY.prop-mu)/Sig^2 + (0.5+U.prop*Om)*as.vector( dSPterm(YY.prop,knots)%*%c(Phi,Gam) )
    Prob <- exp(-V.prop+V)*dnorm(YY, YY.prop-h*dV.prop, sqrt(2*h)) / dnorm(YY.prop, YY-h*dV, sqrt(2*h))
    Prob[Prob>1] <- 1
    ch <- rbinom(n, 1, Prob)
    YY[S==0 & ch==1] <- YY.prop[S==0 & ch==1]
    Y.pos[k,] <- YY
    # Beta
    m1 <- IX%*%t(X)%*%YY
    m2 <- Sig^2*IX
    Beta <- mvrnorm(1,m1,m2)
    Beta.pos[k,] <- Beta
    # Sig
    Mu <- as.vector(X%*%Beta)
    resid <- YY-Mu
    Sig <- sqrt(rinvgamma(1, ccg+n/2, ccg+sum(resid^2)/2))
    Sig.pos[k] <- Sig
    # Phi & Gam
    mat <- diag(c(rep(ccn, q+1), rep(Lam, K)))
    PP <- SPterm(YY, knots)
    A <- solve(t(PP)%*%diag(Om)%*%PP+mat)
    mm <- t(PP)%*%( (S-0.5)-diag(Om)%*%(Zreg+RB) )
    rn <- mvrnorm(1, as.vector(A%*%mm), A)
    Phi <- rn[1:(q+1)]
    Gam <- rn[-(1:(q+1))]
    Phi.pos[k,] <- Phi
    Gam.pos[k,] <- Gam
    # Delta and Xi
    if(NP){
      mat <- diag(c(rep(ccn, rr), rep(Lam2, K)))
      EZ <- cbind(Z, ZZ)
      A <- solve(t(EZ)%*%diag(Om)%*%EZ + mat)
      mm <- t(EZ)%*%( (S-0.5) - diag(Om)%*%PP%*%c(Phi, Gam) )
      rn <- mvrnorm(1, as.vector(A%*%mm), A)
      Delta <- rn[1:rr]
      Xi <- rn[-(1:rr)]
      Delta.pos[k,] <- Delta
      Xi.pos[k,] <- Xi
    }else{
      A <- solve( t(Z)%*%diag(Om)%*%Z + ccn*diag(rr) )
      mm <- t(Z)%*%( (S-0.5) - diag(Om)%*%PP%*%c(Phi,Gam) )
      Delta <- mvrnorm(1, as.vector(A%*%mm), A)
      Delta.pos[k,] <- Delta
    }
    
    if(NP){ RB <- as.vector(ZZ%*%Xi) }
    Zreg <- as.vector(Z%*%Delta)
    
    # response probability
    Pi.pos[k,] <- logistic( as.vector( SPterm(YY, knots)%*%c(Phi,Gam) + Zreg + RB ) )
    
    # Lam (tuning parameter in spline)
    Lam <- rgamma(1, ccg+0.5*K, ccg+0.5*sum(Gam^2))
    Lam.pos[k] <- Lam
    
    # Lam2 (tuning parameters in radial basis part)
    if(NP){
      Lam2 <- rgamma(1, ccg+0.5*K, ccg+0.5*sum(Xi^2))
    }
    
    # Knots
    if(Knot){
      cc <- 0.1
      new.aa <- aa + cc*rnorm(1)
      if(new.aa<0){ new.aa <- 0 }
      if(new.aa>1){ new.aa <- 1 }
      new.knots <- Set.Knots(new.aa)
      U1 <- as.vector( SPterm(YY,knots)%*%c(Phi,Gam) + Zreg + RB )
      U2 <- as.vector( SPterm(YY,new.knots)%*%c(Phi,Gam) + Zreg + RB )
      L1 <- sum((S-0.5)*U1-0.5*Om*U1^2)
      L2 <- sum((S-0.5)*U2-0.5*Om*U2^2)
      prob <- min(1,exp(L2-L1))
      ch <- rbinom(1,1,prob)
      aa <- ch*new.aa+(1-ch)*aa
      a.pos[k] <- aa
      knots <- Set.Knots(aa)
    }
  }
  
  ## Summary
  om <- 1:burn
  Beta.pos <- Beta.pos[-om,]
  Sig.pos <- Sig.pos[-om]
  Y.pos <- Y.pos[-om,]
  Phi.pos <- Phi.pos[-om,]
  Delta.pos <- as.matrix(Delta.pos[-om,])
  Gam.pos <- Gam.pos[-om,]
  a.pos <- a.pos[-om]
  Pi.pos <- Pi.pos[-om,]
  if(NP){ Xi.pos <- Xi.pos[-om,] }
  
  # complete DIC
  hBeta <- apply(Beta.pos, 2, mean)
  hsig <- mean(Sig.pos)
  hY <- apply(Y.pos, 2, mean)
  hPhi <- apply(Phi.pos, 2, mean)
  hDelta <- apply(Delta.pos, 2, mean)
  hGam <- apply(Gam.pos, 2, mean)
  ha <- mean(a.pos)
  RB <- 0
  if(NP){ 
    hXi <- apply(Xi.pos, 2, mean) 
    RB <- as.vector(ZZ%*%hXi)
  }
  knots <- Set.Knots(ha)
  Zreg <- as.vector(Z%*%hDelta)
  hPi <- logistic( as.vector( SPterm(hY, knots)%*%c(hPhi,hGam) + Zreg + RB ) )
  ep <- 10^(-10)
  hPi[hPi<ep] <- ep
  hPi[hPi>1-ep] <- 1-ep
  
  hY <- apply(Y.pos, 2, mean)
  mu.pos <- Beta.pos%*%t(X)
  hmu <- apply(mu.pos, 2, mean)
  hsig <- mean(Sig.pos)
  
  D1 <- -2*( dnorm(Y.pos, mu.pos, Sig.pos, log=T) + t(t(log(Pi.pos))*S) + t(t(log(1-Pi.pos))*(1-S)) )
  D1 <- apply(na.omit(D1), 2, mean)
  mc2 <- mc-burn
  Pi2.pos <- matrix(NA, n, mc2)
  for(k in 1:mc2){
    Pi2.pos[,k] <- logistic( as.vector( SPterm(Y.pos[k,], knots)%*%c(hPhi,hGam) + Zreg + RB ) )
  }
  Pi2.pos[Pi2.pos<ep] <- ep
  Pi2.pos[Pi2.pos>1-ep] <- 1-ep
  D2 <- -2*( dnorm(t(Y.pos), hmu, hsig, log=T) + log(Pi2.pos)*S + log(1-Pi2.pos)*(1-S) )
  D2 <- apply(D2, 1, mean)
  DIC <- sum(2*D1 - D2)
  
  ## summary 
  Aux <- list(Phi=Phi.pos, Delta=Delta.pos, Gam=Gam.pos, Xi=Xi.pos, aa=a.pos, prob=Pi.pos, DIC=DIC)
  Result <- list(Beta=Beta.pos, Sig=Sig.pos, YY=Y.pos, Aux=Aux)
  return(Result)
}









###  Linear selection model  ###
LR.LM <- function(obY, X, Z, S, mc=5000, burn=2000){
  n <- length(obY)
  p <- dim(X)[2]
  rr <- dim(Z)[2]
  
  logistic <- function(x){ 1/(1+exp(-x)) }
  
  ccn <- 0.0001   # hyperparameter (normal prior)
  ccg <- 1   # hyperparameter (gamma prior)
  IX <- solve(t(X)%*%X)
  
  Y.pos <- matrix(NA,mc,n)
  Beta.pos <- matrix(NA,mc,p)
  Sig.pos <- c()
  Phi.pos <- matrix(NA,mc,2)
  Delta.pos <- matrix(NA,mc,rr)
  Om.pos <- matrix(NA,mc,n)
  
  # Initial values
  YY <- c()
  YY[S==0] <- 0
  YY[S==1] <- obY[S==1]
  Beta <- rep(0,p)
  Sig <- 1
  Phi <- rep(0.1, 2)
  Delta <- rep(0, rr)
  Om <- rep(1,n)
  
  # MCMC
  for(k in 1:mc){
    # Omega
    U <- as.vector( Phi[1]+Phi[2]*YY+Z%*%Delta )
    Om <- pgdraw(rep(1,n), U)
    Om.pos[k,] <- Om
    # Y
    mu <- as.vector(X%*%Beta)
    VV <- 1/(Om*Phi[2]^2+1/Sig^2)
    mm <- VV*(-Om*Phi[1]*Phi[2]-0.5*Phi[2]+mu/Sig^2)
    prop <- rnorm(n,mm,sqrt(VV))
    YY[S==0] <- prop[S==0]
    Y.pos[k,] <- YY
    # Beta
    m1 <- IX%*%t(X)%*%YY
    m2 <- Sig^2*IX
    Beta <- mvrnorm(1,m1,m2)
    Beta.pos[k,] <- Beta
    # Sig
    Mu <- as.vector(X%*%Beta)
    resid <- YY-Mu
    Sig <- sqrt(rinvgamma(1,ccg+n/2,ccg+sum(resid^2)/2))
    Sig.pos[k] <- Sig
    # Phi 
    PP <- cbind(1,YY)
    A <- solve(t(PP)%*%diag(Om)%*%PP+ccn*diag(2))
    mm <- t(PP)%*%((S-0.5)-diag(Om)%*%Z%*%Delta)
    Phi <- mvrnorm(1,as.vector(A%*%mm),A)
    Phi.pos[k,] <- Phi
    # Delta
    A <- solve(t(Z)%*%diag(Om)%*%Z+ccn*diag(rr))
    mm <- t(Z)%*%((S-0.5)-diag(Om)%*%PP%*%Phi)
    Delta <- mvrnorm(1,as.vector(A%*%mm),A)
    Delta.pos[k,] <- Delta
    reg <- as.vector(Z%*%Delta)
  }
  
  ## Summary
  om <- 1:burn
  Beta.pos <- Beta.pos[-om,]
  Sig.pos <- Sig.pos[-om]
  Y.pos <- Y.pos[-om,]
  Phi.pos <- Phi.pos[-om,]
  Delta.pos <- as.matrix(Delta.pos[-om,])
  
  # complete DIC
  hPhi <- apply(Phi.pos, 2, mean)
  hDelta <- apply(Delta.pos, 2, mean)
  mu.pos <- Beta.pos%*%t(X)
  hmu <- apply(mu.pos, 2, mean)
  hsig <- mean(Sig.pos)
  hY <- apply(Y.pos, 2, mean)
  PP.pos <- logistic( Phi.pos[,1] + Phi.pos[,2]*Y.pos + as.matrix(Delta.pos)%*%t(Z) )
  hPP <- as.vector(logistic( hPhi[1] + hPhi[2]*hY + Z%*%hDelta ))
  
  D1 <- -2*( dnorm(Y.pos, mu.pos, Sig.pos, log=T) + t(t(log(PP.pos))*S) + t(t(log(1-PP.pos))*(1-S)) )
  D1 <- apply(D1, 2, mean)
  D2 <- -2*( dnorm(hY, hmu, hsig, log=T) + log(hPP)*S + log(1-hPP)*(1-S) )
  DIC <- sum(2*D1 - D2)
  
  # Summary
  Aux <- list(Phi=Phi.pos, Delta=Delta.pos, DIC=DIC)
  Result <- list(Beta=Beta.pos, Sig=Sig.pos, YY=Y.pos, Aux=Aux)
  return(Result)
}







