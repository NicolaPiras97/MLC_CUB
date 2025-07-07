#Code for simulation scenario in the hierarchical case (Simulation 1 and 2)
library(Rcpp)
library(RcppEigen)
library(RcppArmadillo)
Rcpp::sourceCpp("MLC_CUB.cpp") 

L <- 3
H <- 2
K <- 50
C=7
ph <- c(0.4,0.6)
pxw<-matrix(c(0.7,0.1,0.2,0.1,0.3,0.6),nrow=H,ncol=L,byrow=T)
pp1<-c(0.7,0.6,0.8)
pp2<-c(0.6,0.8,0.7)
pp3<-c(0.9,0.6,0.8)
pp4<-c(0.5,0.7,0.6)
pp5<-c(0.8,0.6,0.9)
x1<-c(0.7,0.6,0.2)
x2<-c(0.6,0.8,0.4)
x3<-c(0.1,0.3,0.8)
x4<-c(0.3,0.4,0.6)
x5<-c(0.2,0.1,0.9)

S <- 100
for(s in 1:S){
  
  components<-rep(0,K)
  while(length(which(components==1))!=round(K*ph[1])){ #remove for simulation with random sampling of memberships of level-2 units
    components <- sample(1:H,prob=ph,size=K,replace=TRUE)      
  }
  
  nk<-50 # Simulation 1
  #nk<-10 # Simulation 2
  
  n=K*nk
  
  data <- matrix(nrow=n,ncol=2) 
  data[,1] <- seq(1:n)
  d1<-NULL
  for(k in 1:K){
    d1<-c(d1,rep(k,nk))
  }
  
  data[,2]<-d1
  
  colh<-NULL
  for(i in 1:(K)){
    colh<-c(colh,rep(components[i],nk))
  }
  datac<-cbind(data,colh)
  
  datacc<-datac[order(datac[,3]),]
  data<-datacc[,1:2]
  
  count<-c(length(which(components==1))*nk,length(which(components==2))*nk)
  
  
  data2<-NULL
  samples2 <- NULL
  for(j in (1:H)){
    samples <- sample(1:L,prob=pxw[j,],size=count[j],replace=TRUE) 
    for(i in (1:count[j])){
      mis=c(rbinom(1,1,pp1[samples[i]]),rbinom(1,1,pp2[samples[i]]),rbinom(1,1,pp3[samples[i]]),rbinom(1,1,pp4[samples[i]]),rbinom(1,1,pp5[samples[i]]))
      data2p = cbind(mis[1]*(rbinom(1,size=(C-1),prob=(1-x1[samples[i]]))+1)+(1-mis[1])*sample(1:C,prob=rep(1/C,C),size=1),mis[2]*(rbinom(1,size=(C-1),prob=(1-x2[samples[i]]))+1)+(1-mis[2])*sample(1:C,prob=rep(1/C,C),size=1),mis[3]*(rbinom(1,size=(C-1),prob=(1-x3[samples[i]]))+1)+(1-mis[3])*sample(1:C,prob=rep(1/C,C),size=1),mis[4]*(rbinom(1,size=(C-1),prob=(1-x4[samples[i]]))+1)+(1-mis[4])*sample(1:C,prob=rep(1/C,C),size=1),mis[5]*(rbinom(1,size=(C-1),prob=(1-x5[samples[i]]))+1)+(1-mis[5])*sample(1:C,prob=rep(1/C,C),size=1))
      data2=rbind(data2,data2p)
    } 
    samples2<-c(samples2,samples)  
  }
  
  data<-cbind(data,data2)  
  datac<-cbind(data,samples2)
 
  
  a<-order(data[,2])
  s<-NULL
  for(i in 1:(dim(data)[1])){
    s<-rbind(s,data[a[i],])
  }
  s<-as.matrix(s)
  data<-s
  
  
  y<-list(s[1,])
  for(i in (2:dim(s)[1])){
    y[[i]]<-s[i,]
  }
  
  main2(y)
  
  
}
