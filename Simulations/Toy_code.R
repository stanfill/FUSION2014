library(rotations)
library(reshape2)
library(plyr)
source("Source_Code/robustFunctions.R")
sourceCpp("Source_Code/robustCpp.cpp")
  
S2 <- genR(pi/2)
Qs<-ruarsCont(n=25,rangle=rcayley,kappa1=50,p=.1,Scont=S2,space='Q4')  

plot(Qs,center=id.SO3)

SL2<-mean(Qs)
SL1<-median(Qs)
TSL2<-trimMean(Qs,.1)
Hn<-HnFun(Qs)
WSL2<-weighted.mean(Qs,w=1/Hn)

rot.dist(SL2)
rot.dist(SL1)
rot.dist(TSL2$Shat)
rot.dist(WSL2)

#Longer simulation
S2 <- genR(pi/2)
B<-100
res<-data.frame(ps=rep(c(0,.1,.2),each=B),Mean=0,Median=0,Trim=0,Wei=0)

for(i in 1:nrow(res)){
  
  Qs<-ruarsCont(n=25,rangle=rcayley,kappa1=50,p=res$ps[i],Scont=S2,S=id.SO3,kappa2=50)  

  res$Mean[i]<-rot.dist(mean(Qs))
  res$Median[i]<-rot.dist(median(Qs))
  res$Trim[i]<-rot.dist(trimMean(Qs,res$ps[i]/2,smart=T)$Shat)
  Hns<-HnFun(Qs)
  res$Wei[i]<-rot.dist(weighted.mean(Qs,w=(1/Hns)))
}

resSum<-ddply(res,.(ps),summarize,L2=mean(Mean),L1=mean(Median),TL2=mean(Trim),WL2=mean(Wei))
resSum
MresSum<-melt(resSum,id.vars="ps")

qplot(ps,value,colour=variable,group=variable,geom='line',data=MresSum,xlab=expression(p),ylab='Bias')
