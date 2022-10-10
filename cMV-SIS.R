# PAPER : cMV Survival for Survial Data
# Codes in the matrix version
# For both continous and discrete case

###################################################################
# packages required
###################################################################
rm(list=ls())
gc()
library(MASS)#mvrnorm
library(survival)
library(Matrix)#needed for ahaz
library(ahaz)#compute "fast"
library(survPresmooth)#for IPOD by Hong et al
library(mvtnorm)#for Lq-norm screening
library(prodlim)#for Lq-norm screening


#################################################################################
# functions to compute the criteria in the simulation.
#################################################################################
########################################################################################################################



# 1. M: to compute the minimum model size to ensure the inclusion of all active predictors. 
# 2. mqtl: to compute the 5%, 25%, 50%, 75% and 95% quantiles of the minimum model size out of 1,000 replications.
# 3. Sel.rate: to compute the proportion that every single active predictor is selected 
#    for a given model size, which is defauted c[n/log(n)], in the 1,000 replications.
# 4. Merge.func : to merge two columns by times
# 5. Impute.func : to impute a data in each column with its former value if it's NA




mqtl<-function(M) {
  # Input    
  # M        :  a vector of the minimum model sizes to ensure the inclusion of all active predictors 
  # Output
  # 5%,25%,50%,75%,95% quantiles of minimum model sizes out of 1000 replications
  quantile(M, probs =c(0.05,0.25,0.5,0.75,0.95))
}


M<-function(true.v,rank.mtx) {
  # Input
  # true.v   :  the true variables index
  # rank.mtx :  the ranked index matrix by screening method for the 1000 replications
  #             each column corresponds the ranked index in one replication.
  # Output
  # M        :  a vector of the minimum model sizes to ensure the inclusion of all active predictors 
  r<-min(dim(rank.mtx)[2],length(rank.mtx))##column number
  M<-c()
  for (j in 1:r) {M[j]<-max(match(true.v,rank.mtx[,j]))}
  return(M)
}



Sel.rate<-function(n,c,true.v,rank.mtx) {
  # Input
  # n        :  the sample size
  # c        :  coeficient of cutoffs, for example c=2, cutoff=2[n/log(n)]
  # true.v   :  the true variables index
  # rank.mtx :  the ranked index matrix by screening method for the 1000 replications
  #             each column corresponds the ranked index in one replication.
  # Outputi'y'g't'g'f'f'f'f'f'f'f'f'f'f'f'f'f'f
  # rate     :  the proportions that every single active predictor is selected 
  #             for a given model size, which is defauted c[n/log(n)], in the 1,000 replications.
  d<-c*floor(n/log(n))
  rank.mtx.sel<-rank.mtx[1:d,]
  r<-min(dim(rank.mtx)[2],length(rank.mtx))
  p0<-length(true.v) 
  R<-matrix(0,p0,r)
  rate<-c()
  for (i in 1:p0) {
    for (j in 1:r) {R[i,j]<-(min(abs(rank.mtx.sel[,j]-true.v[i]))==0) }
    rate[i]<-mean(R[i,])
  }
  return(rate)
}


KM= function(Y,de)
{
  n = length(Y)
  
  index<-order(Y,-de)
  
  Y.sort =Y[index]
  de.sort =de[index]
  KM.est<-rep(0,n)
  
  KM.estsort = rep(0,n)
  
  p.vector<-sapply(1:n,function(x) 1-de.sort[x]/(n-x+1))
  
  KM.estsort[1]=ifelse(sum(Y==Y.sort[1])==1,p.vector[1],prod(p.vector[1:sum(Y==Y.sort[1])]))
  
  i=2
  while(i<=n)
    
  {
    
    if(sum(Y==Y.sort[i])==1){
      
      
      KM.estsort[i] = KM.estsort[i-1]*p.vector[i]
      i<- i+1
      
    }else{
      
      l<-sum(Y==Y.sort[i])
      
      pi<-prod(p.vector[i:(i+l-1)])
      
      if(Y.sort[i]==Y.sort[1]) {
        l <- sum(Y==Y.sort[1])-1
        pi <- 1 }
      
      KM.estsort[i:(i+l-1)]<- KM.estsort[i-1]*pi 
      
      i<-i+l
    }
    
  }
  
  KM.estsort[which(KM.estsort==0)]<- 1e-05
  
  KM.est[index] = KM.estsort  # Estimator for survival function of survival time
  
  return(list(KM.est=KM.est,KM.estsort=KM.estsort))
}


#################################################################################
# functions for cMV.SIS
#################################################################################

disc.CMV<- function(j,X,Y,de)
{
  N=nrow(X)
  R=sort(unique(X[,j]))
  ##Compute the fraction
  R.frac = as.numeric(table(X[,j]))/N
  ##define the time interval
  time.pool=numeric(length(R))
  for(r in 1:length(R)){##r
    time.pool[r]=max(Y[X[,j]==R[r]])
  }
  times=seq(min(Y),min(time.pool),.1)
  T0=length(times)
  ##obtain dS=S(t)-S(t-1)
  
  km.all=survfit(Surv(Y,de)~1)
  
  km.all.t=summary(km.all, times=times)
  km.t<-km.all.t$surv
  dS=km.t[-T0]-km.t[-1]
  #get conditional S(t|xj)
  temp.km=array(0,dim=c(length(times),length(R)))
  for(r in 1:length(R)){##r
    sub.Y=Y[X[,j]==R[r]]; sub.de=de[X[,j]==R[r]]
    km=survfit(Surv(sub.Y,sub.de)~ 1)
    km.at.time.t=summary(km, times=times)
    temp.km[,r]<-km.at.time.t$surv
  }
  KM.diff = temp.km[-T0,]-matrix(rep(km.t[-T0],length(R)),nrow=T0-1)
  Intg.KM = colSums( (KM.diff)^2*dS )
  cMVj = sum(Intg.KM*R.frac)
  return(cMVj)
  
}



cont.CMV<-function(j,X,Y,de){
  N=nrow(X)
  Lambda=(3:round(log(N),0))
  out=rep(0,length(Lambda))
  for( i in 1: length(Lambda)){##i
    R=Lambda[i]
    qt=c(quantile(X[,j],probs=c((1:(R-1))/R)))
    #Slice Xj into 1,2,..,R
    index=rep(1,nrow(X))
    for (r in 1:R){##r
      ind=which(X[,j]>=qt[r])
      index[ind]=r+1
    }##r
    #obtain dS
    time.pool=numeric(R)
    for(r in 1:R){##r
      time.pool[r]=max(Y[index==r])
    }
    times=seq(min(Y),min(time.pool),.1)
    T0=length(times)
    km.all=survfit(Surv(Y, de)~ 1)
    km.all.t=summary(km.all, times=times)
    km.t<-km.all.t$surv
    dS=km.t[-T0]-km.t[-1]
    
    #obtain S(t|R=r)
    temp.km=array(0,dim=c(length(times),R))
    for(r in 1:R){##r
      sub.Y=Y[index==r]; sub.de=de[index==r]
      km=survfit(Surv(sub.Y, sub.de)~ 1)
      km.at.time.t=summary(km, times=times)
      temp.km[,r]<-km.at.time.t$surv
    }##r
    
    KM.diff = temp.km[-T0,]-matrix(rep(km.t[-T0],R),nrow=T0-1)
    Intg.KM = colSums( (KM.diff)^2*dS )
    out[i] = sum(Intg.KM*1/R)
    
  }##i
  cMVj=sum(out)
  return(cMVj)
  
}

###Wrapper function
cMVSIS=function(X,Y,de){
  X=scale(X)
  p = ncol(X)
  cep=numeric(p)
  cMV.func=function(j){
    if( length(unique(X[,j]))<10) {cMV.j=disc.CMV(j,X,Y,de)
    } else{cMV.j=cont.CMV(j,X,Y,de)
    }
    return(cMV.j)
  }
  cmv.sis = sapply(1:p,cMV.func)
  return( cmv.sis)
}

#################################################################################
# functions for PSIS: ZhaoLi2012
#################################################################################

PSIS = function(X,Y,de)
{
  n = nrow(X)
  p = ncol(X)
  psis = rep(0,p)
  
  for(j in 1:p)
  {
    xj = as.vector(X[,j])
    cox.fit = coxph(Surv(Y,de)~xj)
    psis[j] = abs(cox.fit$coef)*cox.fit$var^(1/2)
  }
  return(psis)
}





#################################################################################
# functions for CRIS: SongLuMaJeng2014
#################################################################################

CRIS = function(X,Y,de,KM.cenest)
{
  n = nrow(X)
  p = ncol(X)
  rs = rep(0,p)
  
  temp.weight = as.vector(de/KM.cenest^2)
  temp.Y = (matrix(rep(Y,n),nrow=n,byrow=F)<matrix(rep(Y,n),nrow=n,byrow=T))
  
  temp.Y = diag(temp.weight)%*%matrix(as.numeric(temp.Y),nrow=n)
  
  for(j in 1:p)
  {
    x = as.vector(X[,j])
    temp.x = (matrix(rep(x,n),nrow=n,byrow=F)<matrix(rep(x,n),nrow=n,byrow=T))
    temp.x = matrix(as.numeric(temp.x),nrow=n)
    temp.mat = temp.x*temp.Y
    temp.mat[upper.tri(temp.mat)] = 0
    rs[j] = abs(2*sum(temp.mat)/(n*(n-1))-0.25)
  }
  
  
  return(rs)
}



#################################################################################
# functions for CSIRS: ZhouZhu2017
#################################################################################

CSIRS= function(X,Y,de,KM.cenest)
{
  p = dim(X)[2]
  n = length(Y)
  w = rep(0,p)
  X = apply(X,2,function(x) (x-mean(x))/sd(x))
  
  temp.Y = (matrix(rep(Y,n),nrow=n,byrow=F)<matrix(rep(Y,n),nrow=n,byrow=T))
  
  KM.G <- matrix(0,n,p)
  
  for(j in 1:p){
    H<-3
    r<-length(unique(X[,j]))
    
    if(r<5){
      ##Categorical or discrete
      H<-r
      
      xj<-sort(unique(X[,j]))

      for(h in 1:H){
        sub.Y = Y[X[,j]==xj[h]]
        
        sub.de = de[X[,j]==xj[h]]
        
        sub.de.cen=1-sub.de
        
        
        KM.G[X[,j]==xj[h],j]= KM(sub.Y,sub.de.cen)$KM.est
      }
      
    }
    
    else{
      q=c(quantile(X[,j],probs=c((1:(H-1))/3),na.rm = TRUE))
      
      ##create index for subgroup if xj isn't  categorical or discrete
      index<-rep(1,n)
      
      for (r in 1:length(q)){
        ind=which(X[,j]>=q[r])
        index[ind]=r+1
      }
      
      for(h in 1:H){
        
        sub.Y = Y[index==h]
        sub.de = de[index==h]
        
        sub.de.cen=1-sub.de
        
        
        KM.G[index==h,j]= KM(sub.Y,sub.de.cen)$KM.est
      }
      
    }    
    
    w[j]<-mean( apply(t(X[,j]*(de/KM.G[,j]))%*%temp.Y,2,mean)^2)
  }
  
  
  return(w)
}


#################################################################################
# functions for CRSIS: ZhangLiuWu2017
#################################################################################

CRSIS = function(X,KM.est)
{
  n = nrow(X)
  p = ncol(X)
  cr = rep(0,p)
  
  KM.Fest = 1-KM.est
  X = apply(X,2,function(x) (x-mean(x))/sd(x))
  
  for(j in 1:p)
  {
    cr[j] = mean(X[,j]*KM.Fest)^2
  }
  
  return(cr)
}




#################################################################################
# functions for IPOD: written by Hongetal 2017
#################################################################################
## Inputs:
## time: time to death or time to censored
## delta: cenosring indicator
## x: covariates
## gamma: the power of interest 
## note: the bandwidth was set as ((90th percentile of time)- min(time)) /10
## Outputs: screening statisitcs
#---------------------------------------------------------#

#define pairwise.difference ft
pairwise.difference <- function(m){
  npairs <- choose( ncol(m), 2 )
  results <- matrix( NA, nc=npairs, nr=nrow(m) )
  cnames <- rep(NA, npairs)
  if(is.null(colnames(m))) colnames(m) <- paste("col", 1:ncol(m), sep="")
  
  k <- 1
  for(i in 1:ncol(m)){
    for(j in 1:ncol(m)){
      if(j <= i) next;
      results[ ,k] <- m[ ,i] - m[ ,j]
      cnames[k] <- paste(colnames(m)[ c(i, j) ], collapse=".vs.")
      k <- k + 1
    }
  }
  
  colnames(results) <- cnames
  rownames(results) <- rownames(m)
  return(results)
}

##beginning of the program

IPOD.cont= function(j,x,delta,time,gamma){
  #tau=quantile(time, prob=.9)
  
  N=nrow(x); Lambda=(3:round(log(N),0))
  out=array(0,dim=c(length(gamma), length(Lambda)))
  
  for( i in 1: length(Lambda)){##i
    R=Lambda[i] #of slicings, R=3,4,5,6,...
    q=c(quantile(x[,j],probs=c((1:(R-1))/R),na.rm = TRUE))
    
    ##create index for subgroup
    index=rep(1,nrow(x))
    for (r in 1:length(q)){
      ind=which(x[,j]>=q[r])
      index[ind]=r+1
    }
    
    time.pool=numeric(R)
    for(r in 1:R){##r
      time.pool[r]=max(time[index==r])
    }
    t=seq(min(time),min(time.pool),.1)
    h=c(0,diff(range(t))/10)
    
    temp.a=array(0,dim=c(length(t)-1,length(gamma),R))
    R.result=numeric(length(gamma))
    for(r in 1:R){##r
      sub.time=time[index==r]; sub.delta=delta[index==r]
      f=presmooth(sub.time, sub.delta, x.est=t, estimand = "f", bw.selec="fixed",fixed.bw =h)
      
      for(a in 1:length(gamma)){ ##a
        g=f$estimate^gamma[a]
        rec=(t[-1]-t[-length(t)])*g[-length(t)]
        temp.a[,a,r]<- cumsum(rec[1:length(rec)])
      } ##a
    }##r
    
    for( a in 1:length(gamma)){#a
      R.result[a]=max(apply(abs(pairwise.difference (temp.a[,a,])),1,max))
    }##a
    out[,i]=R.result
  }##i
  out=apply(out,1,sum)
  return(out)
}#function

###discrete covariate
IPOD.disc= function(j,x,delta,time,gamma){
  R=sort(unique(x[,j]))
  time.pool=length(R)
  for(r in 1:length(R)){##r
    time.pool[r]=max(time[x[,j]==R[r]])
  }
  t=seq(min(time),min(time.pool),.1)
  h=c(0,diff(range(t))/10)
  out=numeric(length(gamma))
  temp.a=array(0,dim=c(length(t)-1,length(gamma),length(R)))
  for(r in 1:length(R)){##r
    sub.time=time[x[,j]==R[r]]; sub.delta=delta[x[,j]==R[r]]
    f=presmooth(sub.time, sub.delta, x.est=t, estimand = "f", bw.selec="fixed",fixed.bw =h)
    
    for(a in 1:length(gamma)){ ##a
      g=f$estimate^gamma[a]
      rec=(t[-1]-t[-length(t)])*g[-length(t)]
      temp.a[,a,r]<- cumsum(rec[1:length(rec)])
    }##a
  }##r
  for( a in 1:length(gamma)){#a
    out[a]=max(apply(abs(pairwise.difference (temp.a[,a,])),1,max))
  }##a
  return(out)
} 

### wrapper function for continuous and discrete covariates
IPOD=function(x,delta,time,gamma){
  p = ncol(x)
  cep=numeric(p)
  one_model=function(j){
    if( length(unique(x[,j]))<5) {res=IPOD.disc(j,x,delta,time,gamma)
    }else  {res=IPOD.cont(j,x,delta,time,gamma)
    }
  }
  cep = sapply(1:p,one_model)
  return(cep)
}



##==================================================
##Main function
##=================================================


Main1<-function(r,n,p) {# r = simulation times
  MCW<-NULL;   MCWord<-NULL    # MVC-SIS
  PW<-NULL;    PWord<-NULL    # PSIS
  FWord<-NULL    # FAST
  CRW<-NULL;   CRWord<-NULL   #CRIS
  CSW<-NULL;  CSWord<-NULL   #CSIRS1
  CRSW<-NULL;  CRSWord<-NULL#CRSIS
  IP1W<-NULL;  IP1Word<-NULL #IPOD1,gamma=0.8
  IP2W<-NULL;  IP2Word<-NULL#IPOD2,gamma=1
  IP3W<-NULL;  IP3Word<-NULL#IPOD3,gamma=1.7
  
  set.seed(100)
  for (re in 1:r) { 
    ##Create data
    print(re)
    data = get.data(n,p)
    Y = data$Y
    X = data$X
    de = data$de
    Yobj<-Surv(Y,de)
    
    #Break ties
    
    Yobj.breakties<-Surv(Y+runif(n)*1e-5,de)
    
    ##Define the order of Y 
    index = order(Y,-de)
    de.sort = de[index]
    Y.sort = Y[index]
    X.sort = X[index,]
    
    KM.est = KM(Y,de)$KM.est
    
    de.cen = 1-de
    
    KM.cenest = KM(Y,de.cen)$KM.est
    
    ### 2. Independence Screening ### 
    # CMV.SIS
    
    MCW.v<-cMVSIS(X,Y,de)
    MCW<-cbind(MCW,MCW.v)
    MCWord<-cbind( MCWord, order(abs(MCW.v),decreasing = T))
    
    # 
    # PSIS
    PW.v<-PSIS(X,Y,de)
    PW<-cbind(PW,PW.v)
    PWord<-cbind(PWord, order(abs(PW.v),decreasing = T) )
    
    
    ##FAST
    k<-floor(n/log(n))
    FWord.v = ahazisis(Yobj.breakties,X,nsis=k,do.isis=FALSE,rank="FAST")
    FWord<-cbind(FWord ,FWord.v$initRANKorder)
    
    
    #CRIS
    CRW.v<-CRIS(X,Y,de,KM.cenest)
    CRW<-cbind(CRW,CRW.v)
    CRWord<-cbind(CRWord, order(abs(CRW.v),decreasing = T) )
    
    
    #CSIR
    CSW.v<- CSIRS(X,Y,de,KM.cenest)
    CSW<-cbind(CSW,CSW.v)
    CSWord<-cbind(CSWord, order(abs(CSW.v),decreasing = T) )
  
    #CRSIS
    CRSW.v<- CRSIS(X,KM.est)
    
    CRSW<-cbind(CRSW,CRSW.v)
    
    CRSWord<-cbind(CRSWord, order(abs(CRSW.v),decreasing = T) )
    
    
    #IPOD1
    # 
    gamma<-c(0.8,1,1.2)
    
    IPODW.v<-IPOD(X,de,Y,gamma)
    
    IP1W<-cbind(IP1W, IPODW.v[1,])
    
    IP1Word<-cbind(IP1Word, order(abs(IPODW.v[1,]),decreasing = T) )
    
    
    #IPOD2
    IP2W<-cbind(IP2W, IPODW.v[2,])
    IP2Word<-cbind(IP2Word, order(abs(IPODW.v[2,]),decreasing = T) )
    
    #IPOD3
    IP3W<-cbind(IP3W, IPODW.v[3,])
    IP3Word<-cbind(IP3Word, order(abs(IPODW.v[3,]),decreasing = T) )
   }
  print(re)
  
  result=list(MCW=MCW,MCWord=MCWord,PW=PW,PWord=PWord,FWord=FWord,CRW=CRW,CRWord=CRWord,
              CRW=CSW,CSWord=CSWord,CRSW=CRSW,CRSWord=CRSWord,
              IP1W=IP1W,IP1Word=IP1Word,IP2W=IP2W,IP2Word=IP2Word,IP3W=IP3W,IP3Word=IP3Word)
  
  return(result)
}




Main2<-function(n,sim,indx){
  MCM<-M(indx,sim$MCWord)
  PM<-M(indx,sim$PWord)
  FM<-M(indx,sim$FWord)
  CRM<-M(indx,sim$CRWord)
  CSM<-M(indx,sim$CSWord)
  CRSM<-M(indx,sim$CRSWord)
  IP1M<-M(indx,sim$IP1Word)
  IP2M<-M(indx,sim$IP2Word)
  IP3M<-M(indx,sim$IP3Word)

  MC.SR <- Sel.rate(n,1,indx,sim$MCWord)
  PSIS.SR<- Sel.rate(n,1,indx,sim$PWord)
  FAST.SR<-Sel.rate(n,1,indx,sim$FWord)
  CRIS.SR<- Sel.rate(n,1,indx,sim$CRWord)
 
  CSIRS.SR <- Sel.rate(n,1,indx,sim$CSWord)
  
  CRSIS.SR <- Sel.rate(n,1,indx,sim$CRSWord)
  IPod1.SR <- Sel.rate(n,1,indx,sim$IP1Word)
  IPod2.SR <- Sel.rate(n,1,indx,sim$IP2Word)
  IPod3.SR <- Sel.rate(n,1,indx,sim$IP3Word)
 
  
  MC.TR  <- mean(MCM<1*n/log(n))

  P.TR<- mean(PM<1*n/log(n))
  F.TR<-mean(FM<1*n/log(n))
  CRM.TR <- mean(CRM<1*n/log(n))
  CSM.TR <- mean(CSM<1*n/log(n))

  CRSM.TR <- mean(CRSM<1*n/log(n))
  IP1.TR  <- mean(IP1M<1*n/log(n))
  IP2.TR <- mean(IP2M<1*n/log(n))
  IP3.TR <- mean(IP3M<1*n/log(n))

  
  result1<-rbind( cMVSIS =mqtl(MCM),PSIS=mqtl(PM),FAST=mqtl(FM),CRIS=mqtl(CRM),CSIRS=mqtl(CSM),CRSIS=mqtl(CRSM),
                  IPod1=mqtl(IP1M),IPod2=mqtl(IP2M),IPod3=mqtl(IP3M))
  result2<-rbind( cMVSIS =c(MC.SR,MC.TR),PSIS=c(PSIS.SR,P.TR),FAST=c(FAST.SR, F.TR),CRIS=c(CRIS.SR,CRM.TR),CSIRS=c(CSIRS.SR,CSM.TR),
                  CRSIS=c( CRSIS.SR,CRSM.TR),IPod1=c(IPod1.SR,IP1.TR),IPod2=c(IPod2.SR,IP2.TR),IPod3=c(IPod3.SR,IP3.TR))
  
  return(result=list(result1=result1,result2=result2))
}











