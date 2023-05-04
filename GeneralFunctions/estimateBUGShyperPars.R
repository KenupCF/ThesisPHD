###Get hyper-parameters from MCMC posteriors
estimateBUGShyperPars<-function(input,parDist,
                                start=list(),fix.arg=list(),estMethod=list(),defaultEstMethod="mme",
                                method=c("BUGSoutput","coda.samples","jagsUI")){
  
  ###Importing Packages
  require(dplyr)
  require(tidyr)
  require(stringr)
  require(abind)
  require(fitdistrplus)
  require(reshape2)
  options(stringsAsFactors = FALSE)
  ###Parameters
  #`model` is the fitted JAGS model
  #`parDist` is a named vector of the distributions assigned to each followed parameter
  #`start` is a listing  containg starting values for parameter estimation algorithm 
  ######(optional for most distributions)
  
  if(method=="coda.samples"){
    array<-abind(input,along=3)
    array<-aperm(array,c(1,3,2))
    
  }
  
  if(method=="jagsUI"){
    array<-abind(input$samples,along=3)
    array<-aperm(array,c(1,3,2))
  }
  
  if(method=="BUGSoutput"){
    
    array<-input$BUGSoutput$sims.array
  }
  
  
  #Create a data.frame with   
  parInfo<-data.frame(
    EstPar=dimnames(array)[[3]])%>% #the parameters estimated
    mutate(TruePar=gsub(x=EstPar,"\\[.*$",""))%>%         #the `true parameter (unindexed)
    merge(data.frame(TruePar=names(parDist),Distribution=parDist)) #and the correspoding distribution
  
  #Looping through all parameters which we want to know the underlying hyperparameters
  # hyperPars<-lapply(1:nrow(parInfo),function(p){
  hyperPars<-list()
  for(p in 1:nrow(parInfo)){
    
    distrib<-as.character(parInfo$Distribution[p])
    par<-as.character(parInfo$EstPar[p])
    tpar<-as.character(parInfo$TruePar[p])
    
    
    x=as.numeric(array[,,par])
    
    # Change one value of simulated values so it is not a sum-zero - this will estimate a very certain, very small value, instead of straight zero.
    if(sum(x)==0){x[1]<-10e-3}
    if(var(x)==0 & (!distrib%in%c("bern"))){x[1]<-x[1]*(1+10e-3)}
    
    #use function  `fitdistr` to estimate hyperparameters from MCMC values
    fit<-fitdistrplus::fitdist(
      data = x,
      method=ifelse(is.null(estMethod[[tpar]]),defaultEstMethod,estMethod[[tpar]]),
      distr=distrib,
      order=1,
      fix.arg=fix.arg[[tpar]],
      start=start[[tpar]])
    
    temp<-fit$estimate
    fixed<-unlist(fit$fix.arg)
    temp<-c(temp,fixed,Mean=mean(x),median=as.numeric(quantile(x,0.5)),lcl=as.numeric(quantile(x,0.025)),ucl=as.numeric(quantile(x,0.975)),min=min(x),max=max(x))
    #create a data.frame with all hyperparameters of the current parameter
    hyperPars[[p]]<-
      # return(
      data.frame(
        EstPar=par,TruePar=tpar,Distribution=distrib,
        HyperPar=names(temp),Value=temp)
    # )
    # print(p)
    # Sys.sleep(1)
    # })
  }
  
  # p
  #bind and return all data.frames created during loop
  hyperPars<-rbind.fill(hyperPars)
  resu<-tidyr::spread(hyperPars,key=HyperPar,value=Value)
  return(resu)
  
}