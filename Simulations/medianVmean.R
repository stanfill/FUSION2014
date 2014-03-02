library(rotations)
library(plyr)
source("Source_Code/robustFunctions.R")
sourceCpp('Source_Code/MeanMedianAsy.cpp')

Qs<-ruars(20,rcayley,space='Q4',nu=0.25)
L2<-mean(Qs)
L2cd<-cdfunsC(Qs,L2)
SEL2<-L2cd[1]/(2*L2cd[2]^2)

L1<-median(Qs)
L1cd<-cdfunsCMedian(Qs,L1)
SEL1<-L1cd[1]/(2*L2cd[2]^2)

SEL2/SEL1

#####
#Estimate efficiency
B<-500
eff<-rep(0,B)

for(i in 1:B){
  Qs<-ruars(50,rcayley,space='Q4',kappa=1)
  L2<-mean(Qs)
  L2cd<-cdfunsC(Qs,L2)
  SEL2<-L2cd[1]/(2*L2cd[2]^2)
  
  L1<-median(Qs)
  L1cd<-cdfunsCMedian(Qs,L1)
  SEL1<-L1cd[1]/(2*L1cd[2]^2)
  
  eff[i]<-SEL2/SEL1
}

hist(eff,breaks=25)
mean(eff)

####
#Efficiency under contamination
S2 <- genR(pi/2)
B<-100
n<-c(10,25,50)
ps<-c(0,.1,.2)
nps<-expand.grid(n=n,ps=ps)

res<-data.frame(n=rep(nps$n,each=B),ps=rep(nps$ps,each=B),L2Bias=0,L2AV=0,L1Bias=0,L1AV=0)

for(i in 1:nrow(res)){
  
  Rs<-ruarsCont(n=res$n[i],rangle=rcayley,kappa1=50,p=res$ps[i],Scont=S2,S=id.SO3,kappa2=50)  
  Qs<-as.Q4(Rs)
  
  L2<-mean(Qs)
  res$L2Bias[i]<-rot.dist(L2)
  L2cd<-cdfunsC(Qs,L2)
  res$L2AV[i]<-L2cd[1]/(2*L2cd[2]^2)
  
  L1<-median(Qs)
  res$L1Bias[i]<-rot.dist(L1)
  L1cd<-cdfunsCMedian(Qs,L1)
  res$L1AV[i]<-L1cd[1]/(2*L1cd[2]^2)
  
}

#Bias comparison
resSum<-ddply(res,.(n,ps),summarize,L2Bias=mean(L2Bias),L1Bias=mean(L1Bias))
resSum
MresSum<-melt(resSum,id.vars=c("n","ps"))
MresSum$variable<-factor(MresSum$variable,labels=c("Mean","Median"))

qplot(ps,value,data=MresSum,group=variable,colour=variable,geom='line',
  ylab='Estimated Bias',xlab="p",size=I(2))+facet_grid(.~n,labeller=label_bquote(n == .(x)))+
  scale_x_continuous(breaks=c(0,.1,.2))+coord_fixed()+
  guides(colour=guide_legend(title='Estimator'))
  

#Efficiency comparison
resSum2<-ddply(res,.(n,ps),summarize,L2Var=mean(L2AV),L1Var=mean(L1AV),Eff=mean(L2AV/L1AV))
resSum2
MresSum2<-melt(resSum2,id.vars=c("n","ps"))
MresSum2NoEff<-MresSum2[MresSum2$variable!='Eff',]

MresSum2NoEff$variable<-factor(MresSum2NoEff$variable,labels=c("Mean","Median"))

qplot(ps,value,data=MresSum2NoEff,group=variable,colour=variable,geom='line',
  ylab='Estimated Variance',xlab="p",size=I(2))+facet_grid(.~n,labeller=label_bquote(n == .(x)))+
  scale_x_continuous(breaks=c(0,.1,.2))+coord_fixed(2)+
  guides(colour=guide_legend(title='Estimator'))+theme_bw()

qplot(ps,value,data=MresSum2NoEff,group=variable,colour=variable,geom='line',facets=.~n, 
      xlab='Variance',size=I(2))
