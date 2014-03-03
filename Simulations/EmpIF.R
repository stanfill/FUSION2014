#Empirical geometric median influence function
library(rotations)
Rs<-ruars(25,rcayley,kappa=1)

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