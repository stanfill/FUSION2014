#Empirical geometric median influence function
library(rotations)
Rs<-ruars(20,rcayley,kappa=1)

ri<-seq(0,pi,length=100)
med<-median(Rs,type='geometric')
medE<-median(Rs)

diffG<-rep(0,length(ri))
diffE<-rep(0,length(ri))


for(i in 1:length(ri)){
  conti<-med+as.SO3(c(1,0,0),ri[i])
  Rsi<-rbind(Rs,conti)
  class(Rsi)<-'SO3'
  medi<-median(Rsi,type='geometric')
  diffG[i]<-rot.dist(med,medi)
  
  mediE<-median(Rsi)
  diffE[i]<-rot.dist(medE,mediE)
}

plot(ri,diffG,type='l')
lines(ri,diffE,col=2)

##############
#Plot projected IFs
library(ggplot2)
library(plyr)
library(reshape2)
ri<-seq(-pi,pi,length=200)

medianIF<-function(r,kappa){
  a2<-kappa*sqrt(2)*gamma(kappa+2)
  a2<-a2/(3*sqrt(pi)*gamma(kappa+2.5))
  return(sin(r)/(2*a2*sqrt(1-cos(r))))
}

qplot(ri,medianIF(ri,2),geom='line')

meanIF<-function(r,kappa){
  a2<-kappa/(kappa+2)
  return(3*sin(r)/a2)
}

qplot(ri,meanIF(ri,2),geom='line')

IFcomp1<-data.frame(ri=ri,kappa=cayley.kappa(0.25))
IFcomp1$Mean<-meanIF(IFcomp1$ri,IFcomp1$kappa[1])
IFcomp1$Median<-medianIF(IFcomp1$ri,IFcomp1$kappa[1])

IFcomp2<-data.frame(ri=ri,kappa=cayley.kappa(0.5))
IFcomp2$Mean<-meanIF(IFcomp2$ri,IFcomp2$kappa[1])
IFcomp2$Median<-medianIF(IFcomp2$ri,IFcomp2$kappa[1])

IFcomp3<-data.frame(ri=ri,kappa=cayley.kappa(0.75))
IFcomp3$Mean<-meanIF(IFcomp3$ri,IFcomp3$kappa[1])
IFcomp3$Median<-medianIF(IFcomp3$ri,IFcomp3$kappa[1])

#IFcomp<-rbind(IFcomp1,IFcomp2,IFcomp3)
IFcomp<-rbind(IFcomp1,IFcomp3)
mIFcomp<-melt(IFcomp,id=c("ri","kappa"))
mIFcomp$Var<-0.25
mIFcomp[mIFcomp$kappa==2,]$Var<-0.75
mIFcomp$Estimator<-mIFcomp$variable
mIFcomp$Estimator<-factor(mIFcomp$Estimator,labels=c('Proj.Mean','Proj.Median'))  

qplot(ri,value,data=mIFcomp,colour=Estimator,facets=.~Var)+theme_bw()+coord_fixed(1/2)+
  labs(list(x=expression(r),y=expression(IF(R, hat(S)))))+
  geom_hline(yintercept=0,colour='gray50')+geom_vline(xintercept=0,colour='gray50')+geom_line(lwd=I(2))+
  scale_x_continuous(breaks=c(-pi,-pi/2,0,pi/2,pi),
  labels=c(expression(-pi),expression(-pi/2),0,expression(pi/2),expression(pi)))

ggsave("C:/Users/Sta36z/Dropbox/Conferences/FUSION/Submission/Figures/IFComp.pdf",width=9,height=4)

