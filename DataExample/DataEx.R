source("../intervalsSO3/intervals/MedianPaper/PrepareData.R")
sourceCpp('Source_Code/MeanMedianAsy.cpp')
source("Source_Code/robustFunctions.R")
source("Source_Code/modPlot.R")


#####################
#### Find all locations where mean/median differ substantially
#####################

possibles<-which(loc.stats$dE>.1 & loc.stats$dE<.15)
possibles<-sort(c(possibles,which.max(loc.stats$dE)))
m<-length(possibles)
DataExDF<-data.frame(Location=rep(0,m),MeanE=0,MedianE=0,GMedianE=0,WMeanE=0,MeanR=0,MedianR=0,GMedianR=0,WMeanR=0)

for(i in 1:m){
  ex<-loc.stats[possibles[i],]$location 
  DataExDF$Location[i]<-ex
  exRots<-as.SO3(data.matrix(dat.out[dat.out$location==ex,3:11]))
  Qs<-as.Q4(exRots)
  
  L2<-mean(Qs)
  DataExDF$MeanE[i]<-sum(rot.dist(Qs,L2))
  DataExDF$MeanR[i]<-sum(rot.dist(Qs,L2,method='intrinsic'))
  
  L1<-median(Qs)
  DataExDF$MedianE[i]<-sum(rot.dist(Qs,L1))
  DataExDF$MedianR[i]<-sum(rot.dist(Qs,L1,method='intrinsic'))
  
  G1<-median(Qs,type='geometric')
  DataExDF$GMedianE[i]<-sum(rot.dist(Qs,G1))
  DataExDF$GMedianR[i]<-sum(rot.dist(Qs,G1,method='intrinsic'))
  
  Hn<-as.vector(HnCpp(Qs))
  WL1<-weighted.mean(Qs,w=1/Hn)
  DataExDF$WMeanE[i]<-sum(rot.dist(Qs,WL1))
  DataExDF$WMeanR[i]<-sum(rot.dist(Qs,WL1,method='intrinsic'))
}

#Visulize location 698 with 
i<-3
ex<-loc.stats[possibles[i],]$location 
exRots<-as.SO3(data.matrix(dat.out[dat.out$location==ex,3:11]))
plot(exRots,center=median(exRots),col=c(1,2,3))

modPlot(exRots,center=median(exRots),col=c(1),show_estimates='all',size=I(4))+theme(legend.position='none')
#ggsave("C:/Users/Sta36z/Dropbox/Conferences/FUSION/Submission/Figures/Eye_x.png",height=3,width=3)

modPlot(exRots,center=median(exRots),col=c(2),show_estimates='all',size=I(4))+theme(legend.position='none')
#ggsave("C:/Users/Sta36z/Dropbox/Conferences/FUSION/Submission/Figures/Eye_y.png",height=3,width=3)

modPlot(exRots,center=median(exRots),col=c(3),show_estimates='all',size=I(4))+theme(legend.position='none')
#ggsave("C:/Users/Sta36z/Dropbox/Conferences/FUSION/Submission/Figures/Eye_z.png",height=3,width=3)
# DataExDF<-data.frame(Location=rep(0,m),MeanNTH=rep(0,m),MedianNTH=rep(0,m),MeanBoot=rep(0,m),MedianBoot=rep(0,m))
# for(i in 1:m){
#   ex<-loc.stats[possibles[i],]$location 
#   DataExDF$Location[i]<-ex
#   exRots<-as.SO3(data.matrix(dat.out[dat.out$location==ex,3:11]))
#   DataExDF$MeanNTH[i]<-region(exRots,method='moment',type='theory',estimator='mean',alp=.1)*180/pi
#   DataExDF$MedianNTH[i]<-region(exRots,method='moment',type='theory',estimator='median',alp=.1)*180/pi
#   DataExDF$MeanBoot[i]<-region(exRots,method='moment',type='bootstrap',estimator='mean',alp=.1,m=500)*180/pi
#   DataExDF$MedianBoot[i]<-region(exRots,method='moment',type='bootstrap',estimator='median',alp=.1,m=500)*180/pi
# }
# 
# DataExDF$Location<-as.factor(DataExDF$Location)
# xtable(DataExDF,digits=3)



###################
#Use a different method
library(rotations)
library(xtable)
sourceCpp('Source_Code/MeanMedianAsy.cpp')
source("Source_Code/robustFunctions.R")
source("Source_Code/modPlot.R")
data(nickel)
loc111<-nickel[nickel$location==111,]
#loc2975<-nickel[nickel$location==2975,]
exRots<-as.SO3(loc111[,5:13])
exRots

modPlot(exRots,center=median(exRots),col=c(1),show_estimates='all',size=I(4))+theme(legend.position='none')
#ggsave("C:/Users/Sta36z/Dropbox/Conferences/FUSION/Submission/Figures/Eye_x.png",height=3,width=3)

modPlot(exRots,center=median(exRots),col=c(2),show_estimates='all',size=I(4))+theme(legend.position='none')
#ggsave("C:/Users/Sta36z/Dropbox/Conferences/FUSION/Submission/Figures/Eye_y.png",height=3,width=3)

modPlot(exRots,center=median(exRots),col=c(3),show_estimates='all',size=I(4))+theme(legend.position='none')
#ggsave("C:/Users/Sta36z/Dropbox/Conferences/FUSION/Submission/Figures/Eye_z.png",height=3,width=3)

Qs<-as.Q4(exRots)
Hn<-sqrt(HnFun(Qs))
Qstab<-cbind(Qs[order(Hn),],sort(Hn))
xtable(Qstab,digits=3,caption="Example data points expressed as quaternions with corresponding $H_n$ values.")

