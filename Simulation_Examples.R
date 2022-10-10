

##################################################################
##### Example 1 ######
##################################################################
get.data <- function(n, p) {
  X <- matrix(0, n, p)
  for (j in 1:p) {
    X[, j] <- sample(0:3, n,replace = TRUE, prob = c(0.15, 0.2, 0.5, 0.15))
  }
  D <- matrix(0, n, 9)
  
  for (l in 1:3) {
    D[, 3 * l - 2][X[, l] == 3] <- 1
    D[, 3 * l - 1][X[, l] == 2] <- 1
    D[, 3 * l][X[, l] == 1] <- 1
    
  }
  
  beta0 <- c(-1, 1, 0.7,  1, -2, -0.8, 1, -0.9, 0.8)
  Time <- -10 * log(runif(n)) * exp(-D %*% beta0)
  C <- runif(n, 0, 187)
  de = (Time <= C)
  
  print(mean(1 - de))
  Y <- pmin(Time, C)
  return(list(X = X, de = de, Y = Y))
}

R.times = 500
n = 300
p = 1000

indx <- 1:3

Sim5 <- Main1(R.times, n, p)

rzlt5 <- Main2(n, Sim5, indx)

write.csv(rzlt5$result1,file="Z:/User/Documents/SimulationCox12.10/rzlt5_1.csv",row.names=T)

write.table(rzlt5$result2,file="Z:/User/Documents/SimulationCox12.10/rzlt5_2.csv",row.names = T,sep=",",col.names=c("P1","P2","P3","Pa"))




get.data <- function(n, p) {
  X <- matrix(0, n, p)
  for (j in 1:p) {
    X[, j] <- sample(0:3, n,replace = TRUE, prob = c(0.15, 0.2, 0.5,0.15))
  }
  D <- matrix(0, n, 9)
  
  for (l in 1:3) {
    D[, 3 * l - 2][X[, l] == 3] <- 1
    D[, 3 * l - 1][X[, l] == 2] <- 1
    D[, 3 * l][X[, l] == 1] <- 1
    
  }
  
  beta0 <- c(-1, 1, 0.7,  1, -2, -0.8, 1, -0.9, 0.8)
  Time <- -10 * log(runif(n)) * exp(-D %*% beta0)
  C <- runif(n, 0, 53.3)
  de = (Time <= C)
  
  print(mean(1 - de))
  Y <- pmin(Time, C)
  return(list(X = X, de = de, Y = Y))
}

R.times = 500
n = 300
p = 1000

indx <- 1:3

Sim5 <- Main1(R.times, n, p)

per40.rzlt5 <- Main2(n, Sim5, indx)

write.csv(per40.rzlt5$result1,file="Z:/User/Documents/SimulationCox12.10/per40.rzlt5_1.csv",row.names=T)

write.table(per40.rzlt5$result2,file="Z:/User/Documents/SimulationCox12.10/per40.rzlt5_2.csv",row.names = T,sep=",",col.names=c("P1","P2","P3","Pa"))



##################################################################
##### Example 2 ######
##################################################################

######Four categories with different probabilies in a different Cox model

get.data<-function(n,p){
  X <- matrix(0, n, p)
  for (j in 1:p) {
    X[, j] <- sample(0:3, n,replace = TRUE, prob = c(0.2, 0.2, 0.1, 0.5))
  }
  D <- matrix(0, n, 9)
  
  for (l in 1:3) {
    D[, 3 * l - 2][X[, l] == 3] <- 1
    D[, 3 * l - 1][X[, l] == 2] <- 1
    D[, 3 * l][X[, l] == 1] <- 1
    
  }
  
  beta0 <- c(-0.8, 0.9, 1,  -0.6, -1, 0.8, 0.7, -0.9, 0.7)
  
  B<--30*log(runif(n))-8
  Time <- (sign(B)*abs(B)^(1/3)+2) * exp(-D %*% beta0)
  
  
  C <- runif(n, 0, 28.4)
  de = (Time <= C)
  
  print(mean(1 - de))
  Y <- pmin(Time, C)
  return(list(X = X, de = de, Y = Y))
}

n = 300
p = 1000
R.times = 500
indx <- 1:3
Sim6 <- Main1(R.times, n, p)

rzlt6 <- Main2(n, Sim6, indx)

write.csv(rzlt6$result1,file="Z:/User/Documents/SimulationCox12.10/rzlt6_1.csv",row.names=T)

write.table(rzlt6$result2,file="Z:/User/Documents/SimulationCox12.10/rzlt6_2.csv",row.names = T,sep=",",col.names=c("P1","P2","P3","Pa"))




######Four categories with different probabilies in a different Cox model

get.data<-function(n,p){
  X <- matrix(0, n, p)
  for (j in 1:p) {
    X[, j] <- sample(0:3, n,replace = TRUE, prob = c(0.2, 0.2, 0.1, 0.5))
  }
  D <- matrix(0, n, 9)
  
  for (l in 1:3) {
    D[, 3 * l - 2][X[, l] == 3] <- 1
    D[, 3 * l - 1][X[, l] == 2] <- 1
    D[, 3 * l][X[, l] == 1] <- 1
    
  }
  
  beta0 <- c(-0.8, 0.9, 1,  -0.6, -1, 0.8, 0.7, -0.9, 0.7)
  
  B<--30*log(runif(n))-8
  Time <- (sign(B)*abs(B)^(1/3)+2) * exp(-D %*% beta0)
  
  
  C <- runif(n, 0, 11.25)
  de = (Time <= C)
  
  print(mean(1 - de))
  Y <- pmin(Time, C)
  return(list(X = X, de = de, Y = Y))
}

n = 300
p = 1000
R.times = 500
indx <- 1:3
Sim6 <- Main1(R.times, n, p)

rzlt6 <- Main2(n, Sim6, indx)

write.csv(rzlt6$result1,file="Z:/User/Documents/SimulationCox12.10/per40.rzlt6_1.csv",row.names=T)

write.table(rzlt6$result2,file="Z:/User/Documents/SimulationCox12.10/per40.rzlt6_2.csv",row.names = T,sep=",",col.names=c("P1","P2","P3","Pa"))




##################################################################
##### Example 3 ######
##################################################################

creat.sigma<-function(rho,p){
  Sigma<-matrix(0,p,p)
  for (i in 1:p){
    for (j in 1:p){
      Sigma[i,j]<-rho^(abs(i-j))
    }
  }	
  return(Sigma)
}


get.data<-function(n,p){
  rho=0.5
  sigma = creat.sigma(rho,p)
  X = rmvnorm(n,mean=rep(0,p), sigma=sigma, method="chol")
  for (i in 1:n) {
    for (j in 1:p) {
      tmp <- X[i,j]
      X[i,j] <-0 
      X[i,j] <- -1*(tmp<qnorm(.25))+(tmp>qnorm(.75))
    }
    
  }
  
  
  beta0<-c(1,2,0.9,-1.2,1)
  X0 <- cbind(X[,1],X[,2],2*X[,10],2*X[,20],-2*abs(X[,100]))
  Time<-   exp(-X0 %*% beta0+rnorm(n))
  C <- runif(n, 0, 122)
  de = (Time <= C)
  
  print(mean(1 - de))
  Y<-pmin(Time,C)
  return(list(X=X,de=de,Y=Y))
}
n = 300
p = 1000
R.times = 500

indx <- c(1,2,10,20,100)

Sim7 <- Main1(R.times, n, p)

rzlt7 <- Main2(n, Sim7, indx)

write.csv(rzlt7$result1,file="Z:/User/Documents/Simulation12.10/rzlt7_1.csv",row.names=T)

write.table(rzlt7$result2,file="Z:/User/Documents/Simulation12.10/rzlt7_2.csv",row.names = T,sep=",",col.names=c("P1","P2","P3","P4","P5","Pa"))



creat.sigma<-function(rho,p){
  Sigma<-matrix(0,p,p)
  for (i in 1:p){
    for (j in 1:p){
      Sigma[i,j]<-rho^(abs(i-j))
    }
  }	
  return(Sigma)
}


get.data<-function(n,p){
  rho=0.5
  sigma = creat.sigma(rho,p)
  X = rmvnorm(n,mean=rep(0,p), sigma=sigma, method="chol")
  for (i in 1:n) {
    for (j in 1:p) {
      tmp <- X[i,j]
      X[i,j] <-0 
      X[i,j] <- -1*(tmp<qnorm(.25))+(tmp>qnorm(.75))
    }
    
  }
  
  
  beta0<-c(1,2,0.9,-1.2,1)
  X0 <- cbind(X[,1],X[,2],2*X[,10],2*X[,20],-2*abs(X[,100]))
  Time<-   exp(-X0 %*% beta0+rnorm(n))
  C <- runif(n, 0, 17.5)
  de = (Time <= C)
  
  print(mean(1 - de))
  Y<-pmin(Time,C)
  return(list(X=X,de=de,Y=Y))
}
n = 300
p = 1000
R.times = 500

indx <- c(1,2,10,20,100)

Sim7 <- Main1(R.times, n, p)

per40.rzlt7 <- Main2(n, Sim7, indx)

write.csv(per40.rzlt7$result1,file="Z:/User/Documents/Simulation12.10/per40.rzlt7_1.csv",row.names=T)

write.table(per40.rzlt7$result2,file="Z:/User/Documents/Simulation12.10/per40.rzlt7_2.csv",row.names = T,sep=",",col.names=c("P1","P2","P10","P20","P100","Pa"))




creat.sigma<-function(rho,p){
  Sigma<-matrix(0,p,p)
  for (i in 1:p){
    for (j in 1:p){
      Sigma[i,j]<-rho^(abs(i-j))
    }
  }	
  return(Sigma)
}


get.data<-function(n,p){
  rho=0.5
  sigma = creat.sigma(rho,p)
  X = rmvnorm(n,mean=rep(0,p), sigma=sigma, method="chol")
  for (i in 1:n) {
    for (j in 1:p) {
      tmp <- X[i,j]
      X[i,j] <-0 
      X[i,j] <- -1*(tmp<qnorm(.25))+(tmp>qnorm(.75))
    }
  }
  
  #beta0<-(-1)^(sample(c(0,1),1,replace=T,prob=c( a0.6,0.4)))*(a+abs(rnorm(5)))
  beta0<-c(1,2,0.9,-1.2,1)
  X0 <- cbind(X[,1],X[,2],2*X[,10],2*X[,20],-2*abs(X[,100]))
  Time<-   exp(-X0 %*% beta0+rt(n,2))
  C <- runif(n, 0, 156)
  de = (Time <= C)
  
  print(mean(1 - de))
  Y<-pmin(Time,C)
  return(list(X=X,de=de,Y=Y))
}
n = 300
p = 1000
R.times = 500

indx <- c(1,2,10,20,100)


Sim7 <- Main1(R.times, n, p)

rzlt7 <- Main2(n, Sim7, indx)

write.csv(rzlt7$result1,file="Z:/User/Documents/Simulation12.10.t/rzlt7_1.csv",row.names=T)

write.table(rzlt7$result2,file="Z:/User/Documents/Simulation12.10.t/rzlt7_2.csv",row.names = T,sep=",",col.names=c("P1","P2","P3","P4","P5","Pa"))


creat.sigma<-function(rho,p){
  Sigma<-matrix(0,p,p)
  for (i in 1:p){
    for (j in 1:p){
      Sigma[i,j]<-rho^(abs(i-j))
    }
  }	
  return(Sigma)
}


get.data<-function(n,p){
  rho=0.5
  sigma = creat.sigma(rho,p)
  X = rmvnorm(n,mean=rep(0,p), sigma=sigma, method="chol")
  for (i in 1:n) {
    for (j in 1:p) {
      tmp <- X[i,j]
      X[i,j] <-0 
      X[i,j] <- -1*(tmp<qnorm(.25))+(tmp>qnorm(.75))
    }
  }
  
  #beta0<-(-1)^(sample(c(0,1),1,replace=T,prob=c( a0.6,0.4)))*(a+abs(rnorm(5)))
  beta0<-c(1,2,0.9,-1.2,1)
  X0 <- cbind(X[,1],X[,2],2*X[,10],2*X[,20],-2*abs(X[,100]))
  Time<-   exp(-X0 %*% beta0+rt(n,2))
  C <- runif(n, 0, 18)
  de = (Time <= C)
  
  print(mean(1 - de))
  Y<-pmin(Time,C)
  return(list(X=X,de=de,Y=Y))
}
n = 300
p = 1000
R.times = 500

indx <- c(1,2,10,20,100)

Sim7 <- Main1(R.times, n, p)

per40.rzlt7 <- Main2(n, Sim7, indx)

write.csv(per40.rzlt7$result1,file="Z:/User/Documents/Simulation12.10.t/per40.rzlt7_1.csv",row.names=T)

write.table(per40.rzlt7$result2,file="Z:/User/Documents/Simulation12.10.t/per40.rzlt7_2.csv",row.names = T,sep=",",col.names=c("P1","P2","P3","P4","P5","Pa"))



##################################################################
##### Example 4 ######
##################################################################

creat.sigma<-function(rho,p){
  Sigma<-matrix(0,p,p)
  for (i in 1:p){
    for (j in 1:p){
      Sigma[i,j]<-rho^(abs(i-j))
    }
  }	
  return(Sigma)
}
g1<-function(x){
  4*cos(2*pi*x)
}

g2<-function(x){
  -0.3*(x-1)^3
} 

g3<-function(x){
  5*exp(1.2*x-1)-3
}
g4<-function(x) {
  2*atan(3*x-2)
}

g5<-function(x){
  10*(exp(-(3*x-1)^2)+exp(-4*(x-3)^2))-1.5
}


get.data<-function(n,p){
  
  rho=0.5
  sigma = creat.sigma(rho,p)
  X = rmvnorm(n,mean=rep(0,p), sigma=sigma, method="chol")
  Time<- exp( g1(X[,1])+g2(X[,2])+g3(X[,3])+g4(X[,4])+g5(X[,5])+rt(n,2) )
  
  C <- runif(n, 0, 2500)##
  de = (Time <= C)
  
  print(mean(1 - de))          
  Y<-pmin(Time,C)
  return(list(X=X,de=de,Y=Y))
}

n = 300
p = 1000
R.times = 500

indx <- 1:5

Sim10 <- Main1(R.times, n, p)

rzlt10<- Main2(n, Sim10, indx)

write.csv(rzlt10$result1,file="Z:/User/Documents/Simulation12.10.t/rzlt10_1.csv",row.names=T)

write.table(rzlt10$result2,file="Z:/User/Documents/Simulation12.10.t/rzlt10_2.csv",row.names = T,sep=",",col.names=c("P1","P2","P3","P4","P5","Pa"))




##################################################################
##Continous+Another nonlinear model  ######
##################################################################

creat.sigma<-function(rho,p){
  Sigma<-matrix(0,p,p)
  for (i in 1:p){
    for (j in 1:p){
      Sigma[i,j]<-rho^(abs(i-j))
    }
  }	
  return(Sigma)
}
g1<-function(x){
  4*cos(2*pi*x)
}

g2<-function(x){
  -0.3*(x-1)^3
} 

g3<-function(x){
  5*exp(1.2*x-1)-3
}
g4<-function(x) {
  2*atan(3*x-2)
}

g5<-function(x){
  10*(exp(-(3*x-1)^2)+exp(-4*(x-3)^2))-1.5
}


get.data<-function(n,p){
  
  rho=0.5
  sigma = creat.sigma(rho,p)
  X = rmvnorm(n,mean=rep(0,p), sigma=sigma, method="chol")
  Time<- exp( g1(X[,1])+g2(X[,2])+g3(X[,3])+g4(X[,4])+g5(X[,5])+rt(n,2) )
  
  C <- runif(n, 0, 20)##
  de = (Time <= C)
  
  print(mean(1 - de))          
  Y<-pmin(Time,C)
  return(list(X=X,de=de,Y=Y))
}

n = 300
p = 1000
R.times = 500

indx <- 1:5

Sim10 <- Main1(R.times, n, p)

per40.rzlt10<- Main2(n, Sim10, indx)

write.csv(per40.rzlt10$result1,file="Z:/User/Documents/Simulation12.10.t/per40.rzlt10_1.csv",row.names=T)

write.table(per40.rzlt10$result2,file="Z:/User/Documents/Simulation12.10.t/per40.rzlt10_2.csv",row.names = T,sep=",",col.names=c("P1","P2","P3","P4","P5","Pa"))


creat.sigma<-function(rho,p){
  Sigma<-matrix(0,p,p)
  for (i in 1:p){
    for (j in 1:p){
      Sigma[i,j]<-rho^(abs(i-j))
    }
  }	
  return(Sigma)
}
g1<-function(x){
  4*cos(2*pi*x)
}

g2<-function(x){
  -0.3*(x-1)^3
} 

g3<-function(x){
  5*exp(1.2*x-1)-3
}
g4<-function(x) {
  2*atan(3*x-2)
}

g5<-function(x){
  10*(exp(-(3*x-1)^2)+exp(-4*(x-3)^2))-1.5
}


get.data<-function(n,p){
  
  rho=0.5
  sigma = creat.sigma(rho,p)
  X = rmvnorm(n,mean=rep(0,p), sigma=sigma, method="chol")
  Time<- exp( g1(X[,1])+g2(X[,2])+g3(X[,3])+g4(X[,4])+g5(X[,5])+rt(n,2) )
  
  C <- runif(n, 0, 2500)##
  de = (Time <= C)
  
  print(mean(1 - de))          
  Y<-pmin(Time,C)
  return(list(X=X,de=de,Y=Y))
}

n = 300
p = 1000


R.times = 500

indx <- 1:5

Sim10 <- Main1(R.times, n, p)

rzlt10<- Main2(n, Sim10, indx)

write.csv(rzlt10$result1,file="Z:/User/Documents/Simulation12.10.t/rzlt10_1.csv",row.names=T)

write.table(rzlt10$result2,file="Z:/User/Documents/Simulation12.10.t/rzlt10_2.csv",row.names = T,sep=",",col.names=c("P1","P2","P3","P4","P5","Pa"))




creat.sigma<-function(rho,p){
  Sigma<-matrix(0,p,p)
  for (i in 1:p){
    for (j in 1:p){
      Sigma[i,j]<-rho^(abs(i-j))
    }
  }	
  return(Sigma)
}
g1<-function(x){
  4*cos(2*pi*x)
}

g2<-function(x){
  -0.3*(x-1)^3
} 

g3<-function(x){
  5*exp(1.2*x-1)-3
}
g4<-function(x) {
  2*atan(3*x-2)
}

g5<-function(x){
  10*(exp(-(3*x-1)^2)+exp(-4*(x-3)^2))-1.5
}


get.data<-function(n,p){
  
  rho=0.5
  sigma = creat.sigma(rho,p)
  X = rmvnorm(n,mean=rep(0,p), sigma=sigma, method="chol")
  Time<- exp( g1(X[,1])+g2(X[,2])+g3(X[,3])+g4(X[,4])+g5(X[,5])+rt(n,2) )
  
  C <- runif(n, 0, 20)##
  de = (Time <= C)
  
  print(mean(1 - de))          
  Y<-pmin(Time,C)
  return(list(X=X,de=de,Y=Y))
}

n = 300
p = 1000
R.times = 500

indx <- 1:5

Sim10 <- Main1(R.times, n, p)


per40.rzlt10<- Main2(n, Sim10, indx)

write.csv(per40.rzlt10$result1,file="Z:/User/Documents/Simulation12.10.t/per40.rzlt10_1.csv",row.names=T)

write.table(per40.rzlt10$result2,file="Z:/User/Documents/Simulation12.10.t/per40.rzlt10_2.csv",row.names = T,sep=",",col.names=c("P1","P2","P3","P4","P5","Pa"))







