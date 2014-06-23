library(rotations)
sourceCpp('Source_Code/MeanMedianAsy.cpp')
source("Source_Code/robustFunctions.R")

HnApprox<-function(Qs){
  n<-nrow(Qs)
  Qhat<-mean(Qs)
  Hnia<-rep(0,n)
  SSR<-sum(rot.dist(Qs,Qhat,method='e',p=2))
  for(i in 1:n){
    
    Qsi<-Qs[-i,]
    Qhati<-mean(Qsi)
    SSRi<-sum(rot.dist(Qsi,Qhati,method='e',p=2))
    Hnia[i] <- (n-2)*(SSR - SSRi)/(SSRi)
    
  }
  return(Hnia)
}

#Uncontaminated sample
n<-1000
Qs<-ruars(n,rfisher,space='Q4',kappa=5)#,S=genR(pi/4,space='Q4'))
Hn<-HnFun(Qs)
Hnapprox<-HnApprox(Qs)
fit<-lm(Hnapprox~Hn)

plot(Hn,Hnapprox,main=fit$coefficient)
abline(fit)

x<-seq(0,5,length=20)
plot(ecdf(Hn),xlim=c(0,max(x)))
lines(x,pf(x,3,3*(n-2)),col=2)
1-pf(max(Hn),1,n-2)

#Contaminated sample
Scont<-genR(pi/2)
n<-100
Rs<-ruarsCont(n=n,rangle=rvmises,kappa1=50,p=1/n,Scont=Scont,S=id.SO3,kappa2=50)  
Qs<-as.Q4(Rs)
Hn<-HnFun(Qs)
Hnapprox<-HnApprox(Qs)
fit<-lm(Hnapprox~Hn)

plot(Hn,Hnapprox,main=fit$coefficient)
abline(fit)

x<-seq(0,5,length=20)
plot(ecdf(Hn),xlim=c(0,max(x)))
lines(x,pf(x,1,n-2))

plot(ecdf(HnFun(Qs[1:99,])),xlim=c(0,max(x)))
lines(x,pf(x,1,n-2))

#Uncontaminated sample
x<-seq(0,5,length=20)
n<-10
Qs<-ruars(n,rcayley,space='Q4',kappa=5)#,S=genR(pi/4,space='Q4'))
Hn<-HnFun(Qs)

plot(ecdf(Hn),xlim=c(0,max(x)))
lines(x,pf(x,3,3*(n-2)),col=2)

1-pf(max(Hn),1,n-2)

################################
##Alternative method, maybe?  Not!  
#You're thinking about this wrong, SSR is the sums of squared errors
#reduced.  Dropping observation i, is the full model and leaving it
#in is the reduced.

HnApproxAlt<-function(Qs){
  n<-nrow(Qs)
  Qhat<-mean(Qs)
  Hnia<-rep(0,n)
  SSR<-sum(rot.dist(Qs,Qhat,method='e',p=2))
  for(i in 1:n){
    
    Qsi<-Qs[-i,]
    Qhati<-mean(Qsi)
    SSRi<-sum(rot.dist(Qsi,Qhati,method='e',p=2))
    Hnia[i] <- (n-2)*(SSR - SSRi)/(SSR)
    
  }
  return(Hnia)
}


#Uncontaminated sample
x<-seq(0,5,length=20)
n<-10
Qs<-ruars(n,rfisher,space='Q4',kappa=2)#,S=genR(pi/4,space='Q4'))
Hn<-HnApproxAlt(Qs)

plot(ecdf(Hn),xlim=c(0,max(x)))
lines(x,pf(x,3,3*(n-1)),col=2)

1-pf(max(Hn),1,n-2)

#compare the two
plot(Hn,HnFun(Qs))
abline(0,1)
