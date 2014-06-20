source("../intervalsSO3/intervals/MedianPaper/PrepareData.R")
sourceCpp('Source_Code/MeanMedianAsy.cpp')
source("Source_Code/robustFunctions.R")
source("Source_Code/modPlot.R")

#####################
#### Draw the grain map
#####################

#Grain map based on median off all 14 scans
possibles<-which(loc.stats$dE>.075)
d <- ggplot(loc.stats, aes(xpos, ypos, color=dE1))
d2 <- d + geom_point(size=4) + scale_colour_gradient(expression(d[R](tilde(S)[E], I["3x3"])), low="grey99", high="grey10", limits=c(0, pi), breaks=c( pi/4, pi/2, 3*pi/4), labels=expression( pi/4, pi/2, 3*pi/4)) + 
  theme_bw() + xlab("") + ylab("") + coord_equal() + scale_x_continuous(limits=c(0, 12.5), breaks=seq(0, 12.5, by=2.5), labels=expression(0*mu*m, 2.5*mu*m, 5*mu*m, 7.5*mu*m, 10*mu*m, 12.5*mu*m)) + 
  scale_y_continuous(limits=c(0, 10), breaks=seq(0, 10, by=2.5), labels=expression(0*mu*m, 2.5*mu*m, 5*mu*m, 7.5*mu*m, 10*mu*m)) + 
  geom_point(shape="o", colour="yellow", size=5, data=loc.stats[possibles,])  + 
  theme(plot.margin=unit(rep(0,4), "lines"))
d2
ggsave("C:/Users/Sta36z/Dropbox/Conferences/FUSION/Submission/Figures/Grain_map_with_circles.png",height=8,width=9)


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

modPlot(exRots,center=median(exRots),col=c(1),size=I(2))+theme(legend.position='none')
#ggsave("C:/Users/Sta36z/Dropbox/Conferences/FUSION/Presentation/figures/Eye_x_no_Ests.png",height=4,width=4)

modPlot(exRots,center=median(exRots),col=c(2),show_estimates='all',size=I(4))+theme(legend.position='none')
#ggsave("C:/Users/Sta36z/Dropbox/Conferences/FUSION/Submission/Figures/Eye_y.png",height=3,width=3)

modPlot(exRots,center=median(exRots),col=c(2),size=I(4))+theme(legend.position='none')
#ggsave("C:/Users/Sta36z/Dropbox/Conferences/FUSION/Presentation/figures/Eye_y_no_Ests.png",height=3,width=3)

modPlot(exRots,center=median(exRots),col=c(3),show_estimates='all',size=I(4))+theme(legend.position='none')
#ggsave("C:/Users/Sta36z/Dropbox/Conferences/FUSION/Submission/Figures/Eye_z.png",height=3,width=3)

#modPlot(exRots,center=median(exRots),col=c(3),size=I(4))+theme(legend.position='none')
#ggsave("C:/Users/Sta36z/Dropbox/Conferences/FUSION/Presentation/figures/Eye_z_no_Ests.png",height=3,width=3)

Qs<-as.Q4(exRots)
Hn<-(HnFun(Qs))
Qstab<-cbind(Qs[order(Hn),],sort(Hn))
xtable(Qstab,digits=3,caption="Example data points expressed as quaternions with corresponding $H_n$ values.")


adjLoc1<-nickel[nickel$location==112,]
mean(as.SO3(adjLoc1[,5:13]))

adjLoc2<-nickel[nickel$location==110,]
mean(as.SO3(adjLoc2[,5:13]))

adjLoc3 <- nickel[nickel$location==232,]
mean(as.SO3(adjLoc3[,5:13]))

median(exRots)


####################
#For presentation
library(rotations)
library(xtable)
library(gridExtra)
sourceCpp('Source_Code/MeanMedianAsy.cpp')
source("Source_Code/robustFunctions.R")
source("Source_Code/modPlot.R")
data(nickel)
loc111<-nickel[nickel$location==111,]
#loc2975<-nickel[nickel$location==2975,]
exRots<-as.SO3(loc111[,5:13])
Qs<-as.Q4(exRots)
Hn<-(HnFun(Qs))

modPlot(exRots,center=median(exRots),col=c(1),size=I(2))+theme(legend.position='none')
#ggsave("C:/Users/Sta36z/Dropbox/Conferences/FUSION/Presentation/figures/Eye_x_no_Ests.png",height=3,width=3)

modPlot(exRots,center=median(exRots),col=c(1),size=I(4))+theme(legend.position='none')+xlim(-.25,.4)+ylim(-.4,.4)
#ggsave("C:/Users/Sta36z/Dropbox/Conferences/FUSION/Presentation/figures/Eye_x_no_Ests_zoom.png",height=4,width=3.25)

plot(exRots,center=median(exRots),col=c(1),size=I(4))+
  geom_text(aes(x=0.17,y=.19,label="H[i]%~~%0.5"),parse=T)+
  geom_text(aes(x=0.17,y=-.19,label="H[i]%~~%1.9"),parse=T)+
  theme(legend.position='none')+xlim(-.25,.4)+ylim(-.4,.4)
ggsave("C:/Users/Sta36z/Dropbox/Conferences/FUSION/Presentation/figures/Eye_x_no_Ests_labels.png",height=4,width=3.25)


modPlot(exRots,center=median(exRots),col=c(1),show_estimates='all',size=I(2))+theme(legend.position='none')
#ggsave("C:/Users/Sta36z/Dropbox/Conferences/FUSION/Presentation/figures/Eye_x.png",height=3,width=3)

modPlot(exRots,center=median(exRots),col=c(1),show_estimates='all',size=I(4))+theme(legend.position='none')+xlim(-.25,.4)+ylim(-.4,.4)
#ggsave("C:/Users/Sta36z/Dropbox/Conferences/FUSION/Presentation/figures/Eye_x_zoom.png",height=4,width=3.25)
