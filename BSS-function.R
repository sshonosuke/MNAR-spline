library(MCMCpack)
library(pgdraw)


###  Linear regression models  ###
## Inut
# oby: response vector (missing part is `NA`)
# X: matrix of covariates (the first column should be 1's) in the outcome model
# Z: matrix of covariates in the response model
# S: vector of missing indicator (1 for observed; 0 for missing )
# q: order of spline 
# K: number of knots
# mc: length of MCMC 
# burn: burn-in period
# Knot: locations of knots are adaptively esimated if `T`

BSS.LM=function(oby,X,Z,S,q=2,K=10,mc=5000,burn=2000,Knot=T){
  n=length(oby)
  p=dim(X)[2]
  rr=dim(Z)[2]
  kappa=quantile(na.omit(oby),prob=c(0.05,0.95))
  qq=K+q+1
  
  Set.Knots=function(a){
    kappa1=kappa[1]-a*diff(kappa)/2
    kappa2=kappa[2]+a*diff(kappa)/2
    seq(kappa1,kappa2,length=K)
  }
  
  # Knots (initial)
  knots=Set.Knots(0.2) 
  
  ccn=0.0001   # normal prior
  ccg=1   # gamma prior
  
  Beta.pos=matrix(NA,mc,p)
  Sig.pos=c()
  Phi.pos=matrix(NA,mc,q+1)
  Gam.pos=matrix(NA,mc,K)
  Delta.pos=matrix(NA,mc,rr)
  Lam.pos=c()
  a.pos=c()
  
  # Initial values
  YY=c()
  YY[S==0]=0; YY[S==1]=oby[S==1]
  Beta=rep(0,p); Sig=1
  Phi=rep(0.1,q+1)
  Gam=rep(0.01,K)
  Delta=rep(0,rr)
  Lam=1
  Om=rep(1,n)
  aa=0.2
  
  SPterm=function(x,knots){
    nn=length(x)
    mat1=matrix(NA,nn,q+1)
    for(j in 0:q){ mat1[,j+1]=x^j }
    mat2=matrix(NA,nn,K)
    for(l in 1:K){
      val=(x-knots[l])
      mat2[,l]=(val*ifelse(val>0,1,0))^q
    }
    return(cbind(mat1,mat2))
  }
  
  dSPterm=function(x,knots){
    nn=length(x)
    mat1=matrix(NA,nn,q+1)
    mat1[,1]=0
    for(j in 1:q){ mat1[,j+1]=j*x^(j-1) }
    mat2=matrix(NA,nn,K)
    for(l in 1:K){
      val=(x-knots[l])
      mat2[,l]=q*(val*ifelse(val>0,1,0))^(q-1)
    }
    return(cbind(mat1,mat2))
  }
  
  ## MCMC
  for(k in 1:mc){
    # Omega
    U=as.vector( SPterm(YY,knots)%*%c(Phi,Gam)+Z%*%Delta )
    Om=pgdraw(rep(1,n),U)
    # Y  (Langevin MC)
    h=0.2   # step size
    mu=as.vector(X%*%Beta)
    V=0.5*(YY-mu)^2/Sig^2+0.5*U+0.5*U^2*Om
    dV=(YY-mu)/Sig^2+(0.5+U*Om)*as.vector( dSPterm(YY,knots)%*%c(Phi,Gam) )
    YY.prop=rnorm(n,YY-h*dV,sqrt(2*h))
    U.prop=as.vector( SPterm(YY.prop,knots)%*%c(Phi,Gam)+Z%*%Delta )
    V.prop=0.5*(YY.prop-mu)^2/Sig^2+0.5*U.prop+0.5*U.prop^2*Om
    dV.prop=(YY.prop-mu)/Sig^2+(0.5+U.prop*Om)*as.vector( dSPterm(YY.prop,knots)%*%c(Phi,Gam) )
    Prob=exp(-V.prop+V)*dnorm(YY,YY.prop-h*dV.prop,sqrt(2*h))/dnorm(YY.prop,YY-h*dV,sqrt(2*h))
    Prob[Prob>1]=1; ch=rbinom(n,1,Prob)
    YY[S==0 & ch==1]=YY.prop[S==0 & ch==1]
    # Beta
    m1=solve(t(X)%*%X)%*%t(X)%*%YY
    m2=Sig^2*solve(t(X)%*%X)
    Beta=mvrnorm(1,m1,m2)
    Beta.pos[k,]=Beta
    # Sig
    Mu=as.vector(X%*%Beta)
    resid=YY-Mu
    Sig=sqrt(rinvgamma(1,ccg+n/2,ccg+sum(resid^2)/2))
    Sig.pos[k]=Sig
    # Phi & Gam
    mat=diag(c(rep(ccn,q+1),rep(Lam,K)))
    PP=SPterm(YY,knots)
    A=matrix(NA,qq,qq)
    for(i in 1:qq){
      for(j in 1:qq){ A[i,j]=sum(Om*PP[,i]*PP[,j]) }
    }
    A=solve(A+mat)
    mm=t(PP)%*%((S-0.5)-Om*as.vector(Z%*%Delta))
    rn=mvrnorm(1,as.vector(A%*%mm),A)
    Phi=rn[1:(q+1)]; Gam=rn[-(1:(q+1))]
    Phi.pos[k,]=Phi; Gam.pos[k,]=Gam
    # Delta
    A=matrix(NA,rr,rr)
    for(i in 1:rr){
      for(j in 1:rr){ A[i,j]=sum(Om*Z[,i]*Z[,j]) }
    }
    A=solve(A+ccn*diag(rr))
    mm=t(Z)%*%((S-0.5)-Om*as.vector(PP%*%c(Phi,Gam)))
    Delta=mvrnorm(1,as.vector(A%*%mm),A)
    Delta.pos[k,]=Delta
    # Lam
    Lam=rgamma(1,ccg+0.5*K,ccg+0.5*sum(Gam^2))
    Lam.pos[k]=Lam
    # Knots
    if(Knot){
      cc=0.05; new.aa=aa+cc*rnorm(1)
      if(new.aa<0){ new.aa=0 }
      if(new.aa>1){ new.aa=1 }
      new.knots=Set.Knots(new.aa)
      U1=as.vector( SPterm(YY,knots)%*%c(Phi,Gam)+Z%*%Delta )
      U2=as.vector( SPterm(YY,new.knots)%*%c(Phi,Gam)+Z%*%Delta )
      L1=sum((S-0.5)*U1-0.5*Om*U1^2)
      L2=sum((S-0.5)*U2-0.5*Om*U2^2)
      prob=min(1,exp(L2-L1))
      ch=rbinom(1,1,prob)
      aa=ch*new.aa+(1-ch)*aa
      a.pos[k]=aa
    }
  }
  
  om=1:burn
  Beta.pos=Beta.pos[-om,]
  Sig.pos=Sig.pos[-om]
  Phi.pos=Phi.pos[-om,]
  Gam.pos=Gam.pos[-om,]
  Delta.pos=Delta.pos[-om,]
  Lam.pos=Lam.pos[-om]
  a.pos=a.pos[-om]
  
  Result=list(Beta.pos,Sig.pos,Phi.pos,Gam.pos,Delta.pos,Lam.pos,a.pos)
  names(Result)=c("Beta","Sigma","Phi","Gamma","Delta","Lam","a")
  return(Result)
}











###  Linear Mixed Models  ###
## Inut
# oby: response vector (missing part is `NA`)
# X: matrix of covariates (the first column should be 1's) in the outcome model
# Z: matrix of covariates in the response model
# S: vector of missing indicator (1 for observed; 0 for missing )
# ID: vector of cluster ID
# q: order of spline 
# K: number of knots
# mc: length of MCMC 
# burn: burn-in period
# Knot: locations of knots are adaptively esimated if `T`

BSS.LMM=function(oby,X,Z,S,ID,q=2,K=10,mc=5000,burn=2000,Knot=T){
  n=length(oby)
  p=dim(X)[2]
  rr=dim(Z)[2]
  m=max(ID)
  kappa=quantile(na.omit(oby),prob=c(0.05,0.95))
  qq=K+q+1
  
  Set.Knots=function(a){
    kappa1=kappa[1]-a*diff(kappa)/2
    kappa2=kappa[2]+a*diff(kappa)/2
    seq(kappa1,kappa2,length=K)
  }
  
  # Knots (initial)
  knots=Set.Knots(0.2) 
  
  ccn=0.0001   # normal prior
  ccg=1   # gamma prior
  
  RE.pos=matrix(NA,mc,m)
  Beta.pos=matrix(NA,mc,p)
  Sig.pos=c()
  Tau.pos=c()
  Phi.pos=matrix(NA,mc,q+1)
  Gam.pos=matrix(NA,mc,K)
  Delta.pos=matrix(NA,mc,rr)
  Lam.pos=c()
  a.pos=c()
  
  # Initial values
  YY=c()
  YY[S==0]=0; YY[S==1]=oby[S==1]
  Beta=rep(0,p); Sig=1
  RE=rep(0,m); Tau=1
  Phi=rep(0.1,q+1)
  Gam=rep(0.01,K)
  Delta=rep(0,rr)
  Lam=1
  Om=rep(1,n)
  aa=0.2
  
  SPterm=function(x,knots){
    nn=length(x)
    mat1=matrix(NA,nn,q+1)
    for(j in 0:q){ mat1[,j+1]=x^j }
    mat2=matrix(NA,nn,K)
    for(l in 1:K){
      val=(x-knots[l])
      mat2[,l]=(val*ifelse(val>0,1,0))^q
    }
    return(cbind(mat1,mat2))
  }
  
  dSPterm=function(x,knots){
    nn=length(x)
    mat1=matrix(NA,nn,q+1)
    mat1[,1]=0
    for(j in 1:q){ mat1[,j+1]=j*x^(j-1) }
    mat2=matrix(NA,nn,K)
    for(l in 1:K){
      val=(x-knots[l])
      mat2[,l]=q*(val*ifelse(val>0,1,0))^(q-1)
    }
    return(cbind(mat1,mat2))
  }
  
  ## MCMC
  for(k in 1:mc){
    # Omega
    U=as.vector( SPterm(YY,knots)%*%c(Phi,Gam)+Z%*%Delta )
    Om=pgdraw(rep(1,n),U)
    # Y  (Langevin MC)
    h=0.2
    mu=as.vector(X%*%Beta+RE[ID])
    V=0.5*(YY-mu)^2/Sig^2+0.5*U+0.5*U^2*Om
    dV=(YY-mu)/Sig^2+(0.5+U*Om)*as.vector( dSPterm(YY,knots)%*%c(Phi,Gam) )
    YY.prop=rnorm(n,YY-h*dV,sqrt(2*h))
    U.prop=as.vector( SPterm(YY.prop,knots)%*%c(Phi,Gam)+Z%*%Delta )
    V.prop=0.5*(YY.prop-mu)^2/Sig^2+0.5*U.prop+0.5*U.prop^2*Om
    dV.prop=(YY.prop-mu)/Sig^2+(0.5+U.prop*Om)*as.vector( dSPterm(YY.prop,knots)%*%c(Phi,Gam) )
    d1=-0.5*(YY.prop-h*dV.prop-YY)^2/(2*h)
    d2=-0.5*(YY-h*dV-YY.prop)^2/(2*h)
    LP=exp(-V.prop+V+d1-d2)
    LP[LP>0]=0; Prob=exp(LP)
    ch=rbinom(n,1,Prob)
    YY[S==0 & ch==1]=YY.prop[S==0 & ch==1]
    # Beta
    m1=solve(t(X)%*%X)%*%t(X)%*%(YY-RE[ID])
    m2=Sig^2*solve(t(X)%*%X)
    Beta=mvrnorm(1,m1,m2)
    Beta.pos[k,]=Beta
    # Sig
    Mu=as.vector(X%*%Beta+RE[ID])
    resid=YY-Mu
    Sig=sqrt(rinvgamma(1,ccg+n/2,ccg+sum(resid^2)/2))
    Sig.pos[k]=Sig
    # RE
    resid=YY-Mu
    m1=TT*Tau^2/(Sig^2+TT*Tau^2)*apply(matrix(resid,TT,m),2,mean)
    m2=Sig^2*Tau^2/(Sig^2+TT*Tau^2)
    RE=rnorm(m1,sqrt(m2))
    RE.pos[k,]=RE
    # Tau
    Tau=sqrt(rinvgamma(1,ccg+m/2,ccg+sum(RE^2)/2))
    Tau.pos[k]=Tau
    # Phi & Gam
    mat=diag(c(rep(ccn,q+1),rep(Lam,K)))
    PP=SPterm(YY,knots)
    A=matrix(NA,qq,qq)
    for(i in 1:qq){
      for(j in 1:qq){ A[i,j]=sum(Om*PP[,i]*PP[,j]) }
    }
    A=ginv(A+mat)
    mm=t(PP)%*%((S-0.5)-Om*as.vector(Z%*%Delta))
    rn=mvrnorm(1,as.vector(A%*%mm),A)
    Phi=rn[1:(q+1)]; Gam=rn[-(1:(q+1))]
    Phi.pos[k,]=Phi; Gam.pos[k,]=Gam
    # Delta
    A=matrix(NA,rr,rr)
    for(i in 1:rr){
      for(j in 1:rr){ A[i,j]=sum(Om*Z[,i]*Z[,j]) }
    }
    A=solve(A+ccn*diag(rr))
    mm=t(Z)%*%((S-0.5)-Om*as.vector(PP%*%c(Phi,Gam)))
    Delta=mvrnorm(1,as.vector(A%*%mm),A)
    Delta.pos[k,]=Delta
    # Lam
    Lam=rgamma(1,ccg+0.5*K,ccg+0.5*sum(Gam^2))
    Lam.pos[k]=Lam
    # Knots
    if(Knot){
      cc=0.05; new.aa=aa+cc*rnorm(1)
      if(new.aa<0){ new.aa=0 }
      if(new.aa>1){ new.aa=1 }
      new.knots=Set.Knots(new.aa)
      U1=as.vector( SPterm(YY,knots)%*%c(Phi,Gam)+Z%*%Delta )
      U2=as.vector( SPterm(YY,new.knots)%*%c(Phi,Gam)+Z%*%Delta )
      L1=sum((S-0.5)*U1-0.5*Om*U1^2)
      L2=sum((S-0.5)*U2-0.5*Om*U2^2)
      prob=min(1,exp(L2-L1))
      ch=rbinom(1,1,prob)
      aa=ch*new.aa+(1-ch)*aa
      a.pos[k]=aa
    }
    if(round(k/100)==(k/100)){ print(k) }
  }
  
  om=1:burn
  Beta.pos=Beta.pos[-om,]
  Sig.pos=Sig.pos[-om]
  RE.pos=RE.pos[-om,]
  Tau.pos=Tau.pos[-om]
  Phi.pos=Phi.pos[-om,]
  Gam.pos=Gam.pos[-om,]
  Delta.pos=Delta.pos[-om,]
  Lam.pos=Lam.pos[-om]
  a.pos=a.pos[-om]
  
  Result=list(Beta.pos,Sig.pos,RE.pos,Tau.pos,Phi.pos,Gam.pos,Delta.pos,Lam.pos,a.pos)
  names(Result)=c("Beta","Sigma","RE","Tau","Phi","Gamma","Delta","Lam","a")
  return(Result)
}


