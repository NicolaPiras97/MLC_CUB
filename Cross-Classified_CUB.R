#Code for simulation scenario in the Cross-Classified case
library(Rcpp)
library(RcppEigen)
library(RcppArmadillo)
Rcpp::sourceCpp("CrossClassified_CUB.cpp") 

L <- 4
R <- 2
H <- 3
K <- 50
Q <- 15
pr <- c(0.4,0.6)
ph <- c(0.2,0.3,0.5)

pxwz<-matrix(c(0.4,0.3,0.2,0.1,0.15,0.1,0.35,0.4,0.15,0.35,0.2,0.3,0.3,0.1,0.4,0.2,0.3,0.35,0.15,0.2,0.05,0.25,0.4,0.3),nrow=H*R,ncol=L,byrow=T)

c1=7
c2=5

pp1<-c(0.7,0.6,0.8,0.9)
pp2<-c(0.6,0.8,0.7,0.9)
pp3<-c(0.9,0.6,0.8,0.7)
pp4<-c(0.5,0.7,0.6,0.8)
pp5<-c(0.8,0.6,0.9,0.7)
pp6<-c(0.7,0.6,0.8,0.9)
x1<-c(0.7,0.6,0.2,0.3)
x2<-c(0.6,0.8,0.4,0.2)
x3<-c(0.1,0.3,0.8,0.7)
x4<-c(0.3,0.4,0.6,0.8)
x5<-c(0.2,0.9,0.3,0.7)
x6<-c(0.3,0.8,0.4,0.6)

S <- 100
for(s in 1:S){
  
  components<-rep(0,K)
  while(length(which(components==1))!=round(K*ph[1])){ #remove for simulation scheme with random sampling of memberships of level-2 units
    components <- sample(1:H,prob=ph,size=K,replace=TRUE)      
  }
  
  components2<-rep(0,Q)
  while(length(which(components2==1))!=round(Q*pr[1])){ #remove for simulation scheme with random sampling of memberships of level-2 units
    components2 <- sample(1:R,prob=pr,size=Q,replace=TRUE)      
  }
  
  
  nk<-8
  nq<-8
  
  n=K*nk*Q
  
  data <- matrix(nrow=n,ncol=3) 
  data[,1] <- seq(1:n)
  d1<-NULL
  for(k in 1:K){
    d1<-c(d1,rep(k,nk*Q))
  }
  
  data[,2]<-d1
  d2<-NULL
  for(k in 1:K){
    d2<-c(d2,rep(seq(1:Q),nq))
  }
  
  data[,3]<-d2
  
  colr<-NULL
  colr<-c(colr,rep(components2,nq*K))
  colh<-NULL
  for(i in 1:(K)){
    colh<-c(colh,rep(components[i],nk*Q))
  }
  datac<-cbind(data,colh,colr)
  
  
  datacc<-datac[order(datac[,4],datac[,5]),]
  data<-datacc[,1:3]
  
  count<-c(length(which(components==1))*length(which(components2==1))*nk,length(which(components==1))*length(which(components2==2))*nk,length(which(components==2))*length(which(components2==1))*nk,length(which(components==2))*length(which(components2==2))*nk,length(which(components==3))*length(which(components2==1))*nk,length(which(components==3))*length(which(components2==2))*nk)
  
  
  data2<-NULL
  samples2 <- NULL
  w <- 1
  for(j in (1:H)){
    for(m in (1:R)){
      samples <- sample(1:L,prob=pxwz[w,],size=count[w],replace=TRUE) 
      for(i in (1:count[w])){
        mis=c(rbinom(1,1,pp1[samples[i]]),rbinom(1,1,pp2[samples[i]]),rbinom(1,1,pp3[samples[i]]),rbinom(1,1,pp4[samples[i]]),rbinom(1,1,pp5[samples[i]]),rbinom(1,1,pp6[samples[i]]))
        data2p = cbind(mis[1]*(rbinom(1,size=(c1-1),prob=(1-x1[samples[i]]))+1)+(1-mis[1])*sample(1:c1,prob=rep(1/c1,c1),size=1),mis[2]*(rbinom(1,size=(c1-1),prob=(1-x2[samples[i]]))+1)+(1-mis[2])*sample(1:c1,prob=rep(1/c1,c1),size=1),mis[3]*(rbinom(1,size=(c1-1),prob=(1-x3[samples[i]]))+1)+(1-mis[3])*sample(1:c1,prob=rep(1/c1,c1),size=1),mis[4]*(rbinom(1,size=(c2-1),prob=(1-x4[samples[i]]))+1)+(1-mis[4])*sample(1:c2,prob=rep(1/c2,c2),size=1),mis[5]*(rbinom(1,size=(c2-1),prob=(1-x5[samples[i]]))+1)+(1-mis[5])*sample(1:c2,prob=rep(1/c2,c2),size=1),mis[6]*(rbinom(1,size=(c2-1),prob=(1-x6[samples[i]]))+1)+(1-mis[6])*sample(1:c2,prob=rep(1/c2,c2),size=1))
        data2=rbind(data2,data2p)
      } 
      samples2<-c(samples2,samples)  
      w=w+1    
    }
  }
  data<-cbind(data,data2)  
  datac<-cbind(data,samples2)
  
  setwd("C://Users//nicol//Desktop//statistica//Progetto")
  write.table(data,file="daticubcc.txt",sep=",",na="NA",row.names=F,col.names=F)
  
  
  a<-order(data[,2],data[,3])
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
