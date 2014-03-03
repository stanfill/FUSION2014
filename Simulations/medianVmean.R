library(rotations)
library(plyr)
library(gridExtra)
library(reshape2)
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
S21 <- as.SO3(c(1,0,0),pi/2)
S22 <- as.SO3(c(1,0,0),pi/4)
B<-1000
n<-c(10,25,50)
ps<-c(0,.1,.2)
nu<-c(.25,.5)
nps<-expand.grid(n=n,ps=ps)

res<-data.frame(n=rep(nps$n,each=B),ps=rep(nps$ps,each=B),L2Bias=0,L2AV=0,L1Bias=0,L1AV=0,WSL2=0,GL1Bias=0,Sstar=1)
res2<-res
res2$Sstar<-2

for(i in 1:nrow(res)){
  
  #Use pi/2 contamination
  Rs<-ruarsCont(n=res$n[i],rangle=rcayley,kappa1=50,p=res$ps[i],Scont=S21,S=id.SO3,kappa2=50)  
  Qs<-as.Q4(Rs)
  
  L2<-mean(Qs)
  res$L2Bias[i]<-rot.dist(L2)
  L2cd<-cdfunsC(Qs,L2)
  res$L2AV[i]<-L2cd[1]/(2*L2cd[2]^2)
  
  L1<-median(Qs)
  res$L1Bias[i]<-rot.dist(L1)
  L1cd<-cdfunsCMedian(Qs,L1)
  res$L1AV[i]<-L1cd[1]/(2*L1cd[2]^2)
  
  Hn<-sqrt(HnFun(Qs))
  res$WSL2[i]<-rot.dist(weighted.mean(Qs,w=1/Hn))
  
  GL1<-median(Qs,type='geometric')
  res$GL1Bias[i]<-rot.dist(GL1)
  
  
  #Use pi/4 contamination
  Rs<-ruarsCont(n=res2$n[i],rangle=rcayley,kappa1=50,p=res2$ps[i],Scont=S22,S=id.SO3,kappa2=50)  
  Qs<-as.Q4(Rs)
  
  L2<-mean(Qs)
  res2$L2Bias[i]<-rot.dist(L2)
  L2cd<-cdfunsC(Qs,L2)
  res2$L2AV[i]<-L2cd[1]/(2*L2cd[2]^2)
  
  L1<-median(Qs)
  res2$L1Bias[i]<-rot.dist(L1)
  L1cd<-cdfunsCMedian(Qs,L1)
  res2$L1AV[i]<-L1cd[1]/(2*L1cd[2]^2)
  
  Hn<-sqrt(HnFun(Qs))
  res2$WSL2[i]<-rot.dist(weighted.mean(Qs,w=1/Hn))
  
  GL1<-median(Qs,type='geometric')
  res2$GL1Bias[i]<-rot.dist(GL1)
}

res<-rbind(res,res2)
#Bias comparison
resSum<-ddply(res,.(n,ps,Sstar),summarize,L2Bias=mean(L2Bias),L1Bias=mean(L1Bias),WL2Bias=mean(WSL2),GL1Bias=mean(GL1Bias))
resSum
MresSum<-melt(resSum,id.vars=c("n","ps","Sstar"))
MresSum$variable<-factor(MresSum$variable,labels=c("Proj. Mean","Weighted Mean","Proj. Median","Geom. Median"),
                         levels=c("L2Bias","WL2Bias","L1Bias","GL1Bias"))
MresSum$Sstar<-factor(MresSum$Sstar,labels=c("r==frac(pi,4)","r==frac(pi,2)"),levels=c(2,1))
MresSum$n<-factor(MresSum$n,labels=c("n==10","n==25","n==50"))


qplot(ps,value,data=MresSum,colour=variable,group=variable,geom='line',alpha=I(.75),
  ylab='Estimated Bias',xlab=expression(p),size=I(1.5))+facet_grid(Sstar~n,labeller=label_parsed)+
  scale_x_continuous(breaks=c(0,.1,.2))+coord_fixed()+theme_bw()+
  guides(colour=guide_legend(title='Estimator'))+theme(legend.position='top')+coord_fixed(1/2)+
  theme(panel.margin = unit(.75, "lines"))

#ggsave("C:/Users/Sta36z/Dropbox/Conferences/FUSION/Submission/Figures/BiasComp.pdf",width=8,height=4.8)

#library(xtable)
xtable(resSum[,c(1,2,4,6,5,7)],digits=3)

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
