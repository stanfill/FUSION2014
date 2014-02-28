library(rotations)
source("Source_Code/robustFunctions.R")
sourceCpp("Source_Code/robustCpp.cpp")
  
S2 <- genR(pi/2)
Qs<-ruarsCont(n=25,rangle=rcayley,kappa1=50,p=.2,Scont=S2,space='Q4')  

plot(Qs)

SL2<-mean(Qs)
SL1<-median(Qs)
TSL2<-