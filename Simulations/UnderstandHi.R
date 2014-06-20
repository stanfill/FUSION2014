library(rotations)
sourceCpp('Source_Code/MeanMedianAsy.cpp')
source("Source_Code/robustFunctions.R")

HnApprox<-function(Qs){
  n<-nrow(Qs)
  Qhat<-mean(Qs)
  Hnia<-rep(0,n)
  SSR<-sum(rot.dist(Qs,Qhat,method='i',p=2))
  for(i in 1:n){
    
    Qsi<-Qs[-i,]
    Qhati<-mean(Qsi)
    Hnia[i] <- SSR - sum(rot.dist(Qsi,Qhati,method='e',p=2))
    
  }
  return(Hnia)
}

#Uncontaminated sample
Qs<-ruars(50,rcayley,space='Q4',kappa=1)
Hn<-HnFun(Qs)
Hnapprox<-HnApprox(Qs)
fit<-lm(Hnapprox~Hn)

plot(Hn,Hnapprox,main=fit$coefficient)
abline(fit)

#Contaminated sample
Scont<-genR(pi/2)
Rs<-ruarsCont(n=25,rangle=rcayley,kappa1=50,p=.1,Scont=Scont,S=id.SO3,kappa2=50)  
Qs<-as.Q4(Rs)
Hn<-HnFun(Qs)
Hnapprox<-HnApprox(Qs)
fit<-lm(Hnapprox~Hn)

plot(Hn,Hnapprox,main=fit$coefficient)
abline(fit)
