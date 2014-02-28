ruarsCont<-function(n,rangle,kappa1,kappa2=kappa1,p,S=id.SO3,Scont,space='SO3'){
  
  #n   		- sample size
  #rangle - angular distribution from which to simulate
  #kappa1	- concentration parameter for F data
  #kappa2 - concentration for contaminated data
  #p			- percent of n that will be contaminated
  #S			- central direction of normal data
  #Scont	-	central direction of contaminated data
  #space  - SO3 (default) or quaternions("Q4")
    
  nCont<-floor(p*n)

  nNorm <- n-nCont

  rsNorm<-rangle(nNorm,kappa=kappa1)
  RsNorm<-genR(rsNorm,S)	#Simulate from the normal distribution
  
  if(nCont>0){
    rsCont<-rangle(nCont,kappa=kappa2)
    RsCont<-genR(rsCont,Scont) #Simulated from the contaminated distribution
    Rs<-as.SO3(rbind(RsNorm,RsCont))
  }else{
    Rs<-RsNorm
  }
  
  if(space=='Q4')
    Rs<-as.Q4(Rs)
  
  return(Rs)
  
}

HnFun<-function(Qs,full=TRUE){
  #Compute the statistic proposed by FLW(?) that is a function of the largest eigenvalue
  #when observation i was removed
  #Written for quaternions, so if SO3 is given, make them quaternions
  
  if(class(Qs)=="SO3")
    Qs<-as.Q4(Qs)
  
  if(full){
    Hn<- as.vector(HnCpp(Qs))
  }else{
    Hn<- as.vector(HnCpp(Qs))
  }
  return(Hn)
}

trimMean<-function(Qs,a,anneal=T,smart=F,...){
  #Trim the most extreme a% based on the HnFun results
  #Qs - the sample
  #a - percent of sample to remove
  #anneal - T/F, remove all at once (F) or one at a time (T)
  #... - additional arguements passed to mean function
  if(class(Qs)!="Q4")
    Qs<-as.Q4(Qs)
  
  n<-nrow(Qs)
  nCut<-floor(min(max(0,n*a),n)) #remove at least 0, at most n
  
  #Written for quaternions so change to quaternions if given matrices
  if(class(Qs)=="SO3")
    Qs<-Q4(Qs)
  
  if(nCut==0){
    return(list(Qs=Qs,Shat=mean(Qs)))
  }
  
  if(smart){
    Hn<-HnFun(Qs)
    ord<-order(Hn)
    Qs<-Qs[ord,]
    Hn<-Hn[ord]
    cutAt<-which.max(diff(Hn))+1
    tQs<-Qs[-c(cutAt:length(Hn)),]
    return(list(Qs=tQs,Shat=mean(tQs,...)))
  }
  
  if(anneal){
    
    for(i in 1:nCut){
      Hn<-HnFun(Qs)
      Qs<-Qs[-which.max(Hn),]
    }
    return(list(Qs=Qs,Shat=mean(Qs,...)))
    
  }else{
    Hn<-HnFun(Qs)
    toCut<-which(order(Hn)>(n-nCut))
    tQs<-Qs[-toCut,]
    return(list(Qs=tQs,Shat=mean(tQs,...)))
  }
  
}



trimMeanOld<-function(Qs,a,discordFun,anneal=F,...){
  #Trim the most extreme a% based on the HnFun results
  #Qs - the sample
  #a - percent of sample to remove
  #discordFun - function to identify extreme observations, larger value more extreme obs
  #anneal - T/F, remove all at once (F) or one at a time (T)
  #... - additional arguements passed to mean function
  n<-nrow(Qs)
  nCut<-floor(min(max(0,n*a),n)) #remove at least 0, at most n
  
  #Written for quaternions so change to quaternions if given matrices
  if(class(Qs)=="SO3")
    Qs<-Q4(Qs)
  
  if(nCut==0){
    return(list(Qs=Qs,Shat=mean(Qs)))
  }
  
  if(anneal){
    
    for(i in 1:nCut){
      Hn<-discordFun(Qs)
      Qs<-as.Q4(Qs[-which.max(Hn),])
    }
    return(list(Qs=Qs,Shat=mean(Qs,...)))
    
  }else{
    Hn<-discordFun(Qs)
    toCut<-which(order(Hn)>(n-nCut))
    tQs<-as.Q4(Qs[-toCut,])
    return(list(Qs=tQs,Shat=mean(tQs,...)))
  }
}

