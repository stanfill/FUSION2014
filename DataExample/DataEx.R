
#####################
#### Find all locations where mean/median differ substantially
#####################

#Randomly select a location with enough spread to be interesting
possibles<-which(loc.stats$dE>.1 & loc.stats$dE<.15)
possibles<-sort(c(possibles,which.max(loc.stats$dE)))
m<-length(possibles)
DataExDF<-data.frame(Location=rep(0,m),MeanNTH=rep(0,m),MedianNTH=rep(0,m),MeanBoot=rep(0,m),MedianBoot=rep(0,m))

for(i in 1:m){
  ex<-loc.stats[possibles[i],]$location 
  DataExDF$Location[i]<-ex
  exRots<-as.SO3(data.matrix(dat.out[dat.out$location==ex,3:11]))
  DataExDF$MeanNTH[i]<-region(exRots,method='moment',type='theory',estimator='mean',alp=.1)*180/pi
  DataExDF$MedianNTH[i]<-region(exRots,method='moment',type='theory',estimator='median',alp=.1)*180/pi
  DataExDF$MeanBoot[i]<-region(exRots,method='moment',type='bootstrap',estimator='mean',alp=.1,m=500)*180/pi
  DataExDF$MedianBoot[i]<-region(exRots,method='moment',type='bootstrap',estimator='median',alp=.1,m=500)*180/pi
}

DataExDF$Location<-as.factor(DataExDF$Location)
xtable(DataExDF,digits=3)