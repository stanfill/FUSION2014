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
resBias<-data.frame(ps=rep(c(0,.1,.2),each=B),Mean=0,Median=0,Trim=0,Wei=0)
resSE<-data.frame(ps=rep(c(0,.1,.2),each=B),Mean=0,Median=0,Trim=0,Wei=0)

for(i in 1:nrow(res)){
  
  Qs<-ruarsCont(n=25,rangle=rcayley,kappa1=50,p=resBias$ps[i],Scont=S2,S=id.SO3,kappa2=50)  
  
  L2<-mean(Qs)
  resBias$Mean[i]<-rot.dist(L2)
  resSE$Mean[i]<-bootSE(Qs,mean,100,L2)
  
  L1<-median(Qs)
  resBias$Median[i]<-rot.dist(L1)
  resSE$Median[i]<-bootSE(Qs,median,100,L1)
  
  tL2<-trimMean(Qs,a=resBias$ps[i])
  resBias$Trim[i]<-rot.dist(tL2)
  resSE$Trim[i]<-bootSE(Qs,trimMean,100,tL2,a=resBias$ps[i])
  
  Hns<-HnFun(Qs)
  wL2<-weighted.mean(Qs,w=(1/Hns))
  resBias$Wei[i]<-rot.dist(wL2)
  resSE$Wei[i]<-bootSE(Qs,weighted.mean,100,wL2,weight=T)
}

#Compare estimator bias
resSum<-ddply(resBias,.(ps),summarize,L2=mean(Mean),L1=mean(Median),TL2=mean(Trim),WL2=mean(Wei))
resSum
MresSum<-melt(resSum,id.vars="ps")

qplot(ps,value,colour=variable,group=variable,geom='line',data=MresSum,xlab=expression(p),ylab='Bias')

#Compare estimator standard error
resSESum<-ddply(resSE,.(ps),summarize,L2=mean(Mean),L1=mean(Median),TL2=mean(Trim),WL2=mean(Wei))
resSESum
MresSumSE<-melt(resSESum,id.vars="ps")

qplot(ps,value,colour=variable,group=variable,geom='line',data=MresSumSE,xlab=expression(p),ylab='SE')

#Compare estimator efficiency
efDf<-resSESum
efDf[,2:5]<-efDf[,2:5]^(-1)*efDf[,2]
MefDf<-melt(efDf,id.vars="ps")

qplot(ps,value,colour=variable,group=variable,geom='line',data=MefDf,xlab=expression(p),ylab='Efficiency')
