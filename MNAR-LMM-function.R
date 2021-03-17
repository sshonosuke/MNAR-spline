library(MCMCpack)
library(pgdraw)





###  Bayesian semiparametric modeling  ###
# gam.hp: hyperparameter (shape and scale) of gamma prior
# beta.hp: hyperparameter (precision) of normal prior
SR.LMM <- function(oby, X, Z, S, ID, q=2, K=10, mc=5000, burn=2000, Knot=T, gam.hp=1, beta.hp=10^(-4)){
  n <- length(oby)
  p <- dim(X)[2]
  rr <- dim(Z)[2]
  m <- max(ID)
  qq <- K+q+1
  logistic <- function(x){ 1/(1+exp(-x)) }
  
  kappa <- quantile(na.omit(oby),prob=c(0.05,0.95))
  Set.Knots <- function(a){
    kappa1 <- kappa[1]-a*diff(kappa)/2
    kappa2 <- kappa[2]+a*diff(kappa)/2
    seq(kappa1,kappa2,length=K)
  }
  
  # Knots (initial)
  knots <- Set.Knots(0.2) 
  
  ccn <- beta.hp   # precision of normal prior
  ccg <- gam.hp     # shape and scale of gamma prior
  
  ## MCMC box
  Y.pos <- matrix(NA, mc, n)
  RE.pos <- matrix(NA, mc, m)
  Beta.pos <- matrix(NA, mc, p)
  Sig.pos <- c()
  Tau.pos <- c()
  Phi.pos <- matrix(NA, mc, q+1)
  Gam.pos <- matrix(NA, mc, K)
  Delta.pos <- matrix(NA, mc, rr)
  Lam.pos <- c()
  a.pos <- c()
  Pi.pos <- matrix(NA, mc, n)
  
  ## Initial values
  YY <- c()
  YY[S==0] <- 0
  YY[S==1] <- oby[S==1]
  Beta <- rep(0, p)
  Sig <- 1
  RE <- rep(0,m)
  Tau <- 1
  Phi <- rep(0.1, q+1)
  Gam <- rep(0.01, K)
  Delta <- rep(0, rr)
  Lam <- 1
  Om <- rep(1, n)
  aa <- 0.2
  
  ## functions for MCMC
  SPterm <- function(x, knots){
    nn <- length(x)
    mat1 <- matrix(NA, nn, q+1)
    for(j in 0:q){ mat1[,j+1] <- x^j }
    mat2 <- matrix(NA, nn, K)
    for(l in 1:K){
      val <- (x-knots[l])
      mat2[,l] <- (val*ifelse(val>0, 1, 0))^q
    }
    return(cbind(mat1, mat2))
  }
  
  dSPterm <- function(x, knots){
    nn <- length(x)
    mat1 <- matrix(NA, nn, q+1)
    mat1[,1] <- 0
    for(j in 1:q){ mat1[,j+1] <- j*x^(j-1) }
    mat2 <- matrix(NA, nn, K)
    for(l in 1:K){
      val <- (x-knots[l])
      mat2[,l] <- q*(val*ifelse(val>0, 1, 0))^(q-1)
    }
    return(cbind(mat1, mat2))
  }
  
  
  ## MCMC replications
  for(k in 1:mc){
    # Omega
    U <- as.vector( SPterm(YY,knots)%*%c(Phi,Gam)+Z%*%Delta )
    Om <- pgdraw(rep(1,n), U)
    
    # Y  (Langevin MC)
    h <- 0.2
    mu <- as.vector(X%*%Beta+RE[ID])
    V <- 0.5*(YY-mu)^2/Sig^2+0.5*U+0.5*U^2*Om
    dV <- (YY-mu)/Sig^2+(0.5+U*Om)*as.vector( dSPterm(YY,knots)%*%c(Phi,Gam) )
    YY.prop <- rnorm(n,YY-h*dV,sqrt(2*h))
    U.prop <- as.vector( SPterm(YY.prop,knots)%*%c(Phi,Gam)+Z%*%Delta )
    V.prop <- 0.5*(YY.prop-mu)^2/Sig^2+0.5*U.prop+0.5*U.prop^2*Om
    dV.prop <- (YY.prop-mu)/Sig^2+(0.5+U.prop*Om)*as.vector( dSPterm(YY.prop,knots)%*%c(Phi,Gam) )
    d1 <- -0.5*(YY.prop-h*dV.prop-YY)^2/(2*h)
    d2 <- -0.5*(YY-h*dV-YY.prop)^2/(2*h)
    LP <- exp(-V.prop+V+d1-d2)
    LP[LP>0] <- 0
    Prob <- exp(LP)
    ch <- rbinom(n,1,Prob)
    YY[S==0 & ch==1] <- YY.prop[S==0 & ch==1]
    Y.pos[k,] <- YY
    
    # Beta
    m1 <- solve(t(X)%*%X)%*%t(X)%*%(YY-RE[ID])
    m2 <- Sig^2*solve(t(X)%*%X)
    Beta <- mvrnorm(1, m1, m2)
    Beta.pos[k,] <- Beta
    
    # Sig
    Mu <- as.vector(X%*%Beta)
    resid <- YY-Mu-RE[ID]
    Sig <- sqrt(rinvgamma(1,ccg+n/2,ccg+sum(resid^2)/2))
    Sig.pos[k] <- Sig
    
    # RE
    resid <- YY-Mu
    m1 <- TT*Tau^2/(Sig^2+TT*Tau^2)*apply(matrix(resid, TT, m), 2, mean)
    m2 <- Sig^2*Tau^2/(Sig^2+TT*Tau^2)
    RE <- rnorm(m, m1, sqrt(m2))
    RE.pos[k,] <- RE
    
    # Tau
    Tau <- sqrt( rinvgamma(1, ccg+m/2, ccg+sum(RE^2)/2) )
    Tau.pos[k] <- Tau
    
    # Phi & Gam
    mat <- diag( c(rep(ccn,q+1), rep(Lam,K)) )
    PP <- SPterm(YY, knots)
    A <- matrix(NA, qq, qq)
    for(i in 1:qq){
      for(j in 1:qq){ A[i,j] <- sum(Om*PP[,i]*PP[,j]) }
    }
    A <- ginv(A+mat)
    mm <- t(PP)%*%((S-0.5)-Om*as.vector(Z%*%Delta))
    rn <- mvrnorm(1, as.vector(A%*%mm), A)
    Phi <- rn[1:(q+1)]
    Gam <- rn[-(1:(q+1))]
    Phi.pos[k,] <- Phi
    Gam.pos[k,] <- Gam
    
    # Delta
    A <- matrix(NA,rr,rr)
    for(i in 1:rr){
      for(j in 1:rr){ A[i,j] <- sum(Om*Z[,i]*Z[,j]) }
    }
    A <- solve(A+ccn*diag(rr))
    mm <- t(Z)%*%((S-0.5)-Om*as.vector(PP%*%c(Phi,Gam)))
    Delta <- mvrnorm(1, as.vector(A%*%mm), A)
    Delta.pos[k,] <- Delta
    Zreg <- as.vector(Z%*%Delta)
    
    # response probability
    Pi.pos[k,] <- logistic( as.vector( SPterm(YY, knots)%*%c(Phi,Gam) + Zreg ) )
    
    # Lam
    Lam <- rgamma(1, ccg+0.5*K, ccg+0.5*sum(Gam^2))
    Lam.pos[k] <- Lam
    
    # Knots
    if(Knot){
      cc <- 0.05
      new.aa <- aa+cc*rnorm(1)
      if(new.aa<0){ new.aa <- 0 }
      if(new.aa>1){ new.aa <- 1 }
      new.knots <- Set.Knots(new.aa)
      U1 <- as.vector( SPterm(YY,knots)%*%c(Phi,Gam)+Z%*%Delta )
      U2 <- as.vector( SPterm(YY,new.knots)%*%c(Phi,Gam)+Z%*%Delta )
      L1 <- sum((S-0.5)*U1-0.5*Om*U1^2)
      L2 <- sum((S-0.5)*U2-0.5*Om*U2^2)
      prob <- min(1, exp(L2-L1))
      ch <- rbinom(1, 1, prob)
      aa <- ch*new.aa+(1-ch)*aa
      a.pos[k] <- aa
      knots <- Set.Knots(aa) 
    }
    
    if(round(k/1000)==(k/1000)){ print(k) }
  }
  
  ## Summary
  omit <- 1:burn
  Beta.pos <- Beta.pos[-omit,]
  Sig.pos <- Sig.pos[-omit]
  RE.pos <- RE.pos[-omit,]
  Tau.pos <- Tau.pos[-omit]
  Y.pos <- Y.pos[-omit,]
  Phi.pos <- Phi.pos[-omit,]
  Delta.pos <- Delta.pos[-omit,]
  Gam.pos <- Gam.pos[-omit,]
  a.pos <- a.pos[-omit]
  Pi.pos <- Pi.pos[-omit,]
  
  ## DIC
  hBeta <- apply(Beta.pos, 2, mean)
  hsig <- mean(Sig.pos)
  hY <- apply(Y.pos, 2, mean)
  hPhi <- apply(Phi.pos, 2, mean)
  hDelta <- apply(Delta.pos, 2, mean)
  hGam <- apply(Gam.pos, 2, mean)
  ha <- mean(a.pos)
  knots <- Set.Knots(ha)
  Zreg <- as.vector(Z%*%hDelta)
  hPi <- logistic( as.vector( SPterm(hY, knots)%*%c(hPhi,hGam) + Zreg ) )
  ep <- 10^(-5)
  hPi[hPi<ep] <- ep
  hPi[hPi>1-ep] <- 1-ep
  Pi.pos[Pi.pos<ep] <- ep
  Pi.pos[Pi.pos>1-ep] <- 1-ep
  
  hY <- apply(Y.pos, 2, mean)
  mu.pos <- Beta.pos%*%t(X) + RE.pos[,ID]
  hmu <- apply(mu.pos, 2, mean)
  hsig <- mean(Sig.pos)
  
  D1 <- -2*( dnorm(Y.pos, mu.pos, Sig.pos, log=T) + t(t(log(Pi.pos))*S) + t(t(log(1-Pi.pos))*(1-S)) )
  D1 <- apply(na.omit(D1), 2, mean)
  D2 <- -2*( dnorm(hY, hmu, hsig, log=T) + log(hPi)*S + log(1-hPi)*(1-S) )
  DIC <- sum(2*D1 - D2)
  
  Result <- list(Beta=Beta.pos, Sig=Sig.pos, RE=RE.pos, Tau=Tau.pos, Y=Y.pos, 
                 Phi=Phi.pos, Delta=Delta.pos, Gam=Gam.pos, a=a.pos, DIC=DIC) 
  return(Result)
}







###  Linear selection model  ###
LR.LMM <- function(oby, X, Z, S, ID, mc=5000, burn=2000){
  n <- length(oby)
  p <- dim(X)[2]
  rr <- dim(Z)[2]
  m <- max(ID)
  
  ccn <- 0.0001   # normal prior
  ccg <- 1   # gamma prior
  
  Y.pos <- matrix(NA,mc,n)
  Beta.pos <- matrix(NA,mc,p)
  Sig.pos <- c()
  RE.pos <- matrix(NA,mc,m)
  Tau.pos <- c()
  Phi.pos <- matrix(NA,mc,2)
  Delta.pos <- matrix(NA,mc,rr)
  
  # Initial values
  YY <- rep(0,n)
  YY[S==0] <- 0
  YY[S==1] <- oby[S==1]
  Beta <- rep(0,p)
  Sig <- 1
  RE <- rep(0,m)
  Tau <- 1
  Phi <- rep(0.1,2)
  Delta <- rep(0,rr)
  Om <- rep(1,n)
  
  # MCMC
  for(k in 1:mc){
    # Omega
    U <- as.vector( Phi[1]+Phi[2]*YY+Z%*%Delta )
    Om <- pgdraw(rep(1,n),U)
    # Y
    mu <- as.vector(X%*%Beta+RE[ID])
    VV <- 1/(Om*Phi[2]^2+1/Sig^2)
    mm <- VV*(-Om*Phi[1]*Phi[2]-0.5*Phi[2]+mu/Sig^2)
    prop <- rnorm(n,mm,sqrt(VV))
    YY[S==0] <- prop[S==0]
    Y.pos[k,] <- YY
    # Beta
    m1 <- solve(t(X)%*%X)%*%t(X)%*%(YY-RE[ID])
    m2 <- Sig^2*solve(t(X)%*%X)
    Beta <- mvrnorm(1,m1,m2)
    Beta.pos[k,] <- Beta
    # Sig
    Mu <- as.vector(X%*%Beta)
    resid <- YY-Mu-RE[ID]
    Sig <- sqrt(rinvgamma(1,ccg+n/2,ccg+sum(resid^2)/2))
    Sig.pos[k] <- Sig
    # RE
    resid <- YY-Mu
    m1 <- TT*Tau^2/(Sig^2+TT*Tau^2)*apply(matrix(resid,TT,m),2,mean)
    m2 <- Sig^2*Tau^2/(Sig^2+TT*Tau^2)
    RE <- rnorm(m1,sqrt(m2))
    RE.pos[k,] <- RE
    # Tau
    Tau <- sqrt(rinvgamma(1,ccg+m/2,ccg+sum(RE^2)/2))
    Tau.pos[k] <- Tau
    # Phi 
    PP <- cbind(1,YY)
    mat <- matrix(NA,2,2)
    mat[1,1] <- sum(PP[,1]^2*Om)
    mat[2,2] <- sum(PP[,2]^2*Om)
    mat[1,2] <- mat[2,1] <- sum(PP[,1]*PP[,2]*Om)
    A <- ginv(mat+ccn*diag(2))
    mm <- apply(PP*((S-0.5)-Om*as.vector(Z%*%Delta)),2,sum)
    Phi <- mvrnorm(1,as.vector(A%*%mm),A)
    Phi.pos[k,] <- Phi
    # Delta
    A <- matrix(NA,rr,rr)
    for(i in 1:rr){
      for(j in 1:rr){ A[i,j] <- sum(Om*Z[,i]*Z[,j]) }
    }
    A <- solve(A+ccn*diag(rr))
    mm <- t(Z)%*%((S-0.5)-Om*as.vector(PP%*%Phi))
    Delta <- mvrnorm(1,as.vector(A%*%mm),A)
    Delta.pos[k,] <- Delta
    if(round(k/1000)==(k/1000)){ print(k) }
  }
  
  ## Summary
  omit <- 1:burn
  Beta.pos <- Beta.pos[-omit,]
  Sig.pos <- Sig.pos[-omit]
  RE.pos <- RE.pos[-omit,]
  Tau.pos <- Tau.pos[-omit]
  Y.pos <- Y.pos[-omit,]
  Phi.pos <- Phi.pos[-omit,]
  Delta.pos <- Delta.pos[-omit,]
  
  # DIC
  logistic <- function(x){ 1/(1+exp(-x)) }
  hPhi <- apply(Phi.pos, 2, mean)
  hDelta <- apply(Delta.pos, 2, mean)
  mu.pos <- Beta.pos%*%t(X) + RE.pos[,ID]
  hmu <- apply(mu.pos, 2, mean)
  hsig <- mean(Sig.pos)
  hY <- apply(Y.pos, 2, mean)
  PP.pos <- logistic( Phi.pos[,1] + Phi.pos[,2]*Y.pos + as.matrix(Delta.pos)%*%t(Z) )
  hPP <- as.vector(logistic( hPhi[1] + hPhi[2]*hY + Z%*%hDelta ))
  
  D1 <- -2*( dnorm(Y.pos, mu.pos, Sig.pos, log=T) + t(t(log(PP.pos))*S) + t(t(log(1-PP.pos))*(1-S)) )
  D1 <- apply(D1, 2, mean)
  D2 <- -2*( dnorm(hY, hmu, hsig, log=T) + log(hPP)*S + log(1-hPP)*(1-S) )
  DIC <- sum(2*D1 - D2)
  
  Result <- list(Beta=Beta.pos, Sig=Sig.pos, RE=RE.pos, Tau=Tau.pos, Y=Y.pos, Phi=Phi.pos, Delta=Delta.pos, DIC=DIC) 
  return(Result)
}









