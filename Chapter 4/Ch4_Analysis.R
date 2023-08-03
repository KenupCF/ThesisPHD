### Chapter 4 Analysis

require(plyr)
require(dplyr)
require(magrittr)
require(data.table)
require(ggplot2)
require(fitdistrplus)
require(jagsUI)
require(devtools)
require(abind)
require(reshape2)
require(gtools)
require(lubridate)
require(parallel)
require(truncnorm)

#Running
source_url("https://raw.githubusercontent.com/KenupCF/Kenup-Repo/master/Quick%20Functions.R")
# source_url("https://raw.githubusercontent.com/KenupCF/Kenup-Repo/master/phd_experimental_functions.R")

devtools::source_url("https://raw.githubusercontent.com/KenupCF/ThesisPHD/main/GeneralFunctions/estimateBUGShyperPars.R")
devtools::source_url("https://raw.githubusercontent.com/KenupCF/ThesisPHD/main/GeneralFunctions/ObjectiveScore.R")


#### Cost ###

###Price of stuff

CostOfCameras<-250 #in NZD
CostOfRadios<-500 # in NZD

CostOfSoftSupp<-10e3 #in NZD
CostOfHardSupp<-20e3 #in NZD



# Creating management options
managementOptions<-expand.grid(
  SoftFeed=0:1,
  HardFeed=0:1,
  ReleasedAnimals=c(0,0),
  PoisonDrop=c(0,1),
  ReleaseSchedule=0:2)%>%
  # Removing illogical management alternatives
  dplyr::filter(!(SoftFeed == 1 & HardFeed == 1))%>%
  dplyr::filter(!(ReleasedAnimals == 0 & ReleaseSchedule > 0))%>%
  dplyr::filter(!(ReleasedAnimals > 0 & ReleaseSchedule == 0))

managementOptions<-managementOptions%>%
  dplyr::filter(!duplicated(managementOptions))
  
managementCols<-colnames(managementOptions)

# Dataframe holding possible management options and their cost
managementOptions<-managementOptions%>%
  dplyr::mutate(CostMgmt=(SoftFeed*CostOfSoftSupp)+
                  (HardFeed*CostOfHardSupp)+
                  (ReleasedAnimals*ReleaseSchedule*1000)+
                  (PoisonDrop*30e3),
                mgmt=1:n())
Baseline_index<-which(apply(apply(
  managementOptions[,managementCols],2,function(x){x==0}),
  1,function(x){prod(x)==1}))



# Loading results from VOI analysis for continuation on AM

load(paste0(".\\Chapter 3\\Framework 1\\Chapter 3 - Framework 1 - Results.RData",sep=""))

VOI_backup<-VOI
priors<-VOI_backup$priors


#### Additional effects' priors  on this model 

# Add effect of Mast on Survival (negative)
priors$survMastEffect<-data.frame(mean=-1.7,sd=.25,dist="norm")
# Add effect of a poison drop, which is expected to counter act the effect on Mast
priors$survPoisonDropEffect<-data.frame(mean=+1.7,sd=.45,dist="norm")

# Add cost of release on Survival
priors$survCostRelease<-data.frame(shape1=90,shape2=10,dist="beta")

nSims<-nrow(VOI_backup$TrueParams)
nSims<-1000
nYears<-10
noAgeClasses<-3

# Calculate where the 0 values of a Leslie matrix would go, for a given number of age classes
# (i.e., impossible transitions)
{L<-matrix(0,nrow=noAgeClasses,ncol=noAgeClasses)
  for(i in 1:ncol(L)){
    L[1,i]<-1
    for(j in 1:nrow(L)){
      L[min(i+1,nrow(L)),i]<-1
    }}
  
  # Result is a matrix of n rows (where n is the number of zeroes in the Leslie matrix) an two columns
  # (column 'i' shows the row index of each zero, while column 'j' shows the column index)
  ZeroIndexes<-reshape2::melt(L)%>%
    dplyr::filter(value==0)%>%
    dplyr::mutate(i=Var1,j=Var2)%>%
    dplyr::select(i,j)%>%
    as.matrix
}

# Repeating priors but adding a temporal component to them
priors_timelapse<-lapply(
  # For each prior, run the following
  names(priors),
  function(x){
  
  # Define what is the distribution
  dist<-priors[[x]]$dist
  
  
  # Depending on the distribution, calculate the mean and 95% credible intervals given the priors
  
  if(dist=="beta"){
    Mean<-priors[[x]]$shape1/(priors[[x]]$shape1+priors[[x]]$shape2)
    lcl<-qbeta(p = .025,shape1 = priors[[x]]$shape1,shape2=priors[[x]]$shape2)
    ucl<-qbeta(p = .975,shape1 = priors[[x]]$shape1,shape2=priors[[x]]$shape2)
  }
  
  if(dist=="gamma"){
    Mean<-priors[[x]]$shape/priors[[x]]$rate
    lcl<-qgamma(p = .025,shape = priors[[x]]$shape,rate=priors[[x]]$rate)
    ucl<-qgamma(p = .975,shape = priors[[x]]$shape,rate=priors[[x]]$rate)    
  }
  
  if(dist=="unif"){
    Mean<-(priors[[x]]$min+priors[[x]]$max)/2
    lcl<-qunif(p = .025,min = priors[[x]]$min,max=priors[[x]]$max)
    ucl<-qunif(p = .975,min = priors[[x]]$min,max=priors[[x]]$max)        
  }  
  if(dist=="norm"){
    Mean<-priors[[x]]$mean
    lcl<-qnorm(p = .025,mean = priors[[x]]$mean,sd=priors[[x]]$sd)
    ucl<-qnorm(p = .975,mean = priors[[x]]$mean,sd=priors[[x]]$sd)
  }
  if(dist=="poisson"){
    Mean<-priors[[x]]$lambda
    lcl<-qpois(p = .025,lambda = priors[[x]]$lambda)
    ucl<-qpois(p = .975,lambda = priors[[x]]$lambda)      
  }  
  if(dist=="exp"){
    Mean<-(priors[[x]]$rate^-1)
    lcl<-qexp(p = .025,rate = priors[[x]]$rate)
    ucl<-qexp(p = .975,rate = priors[[x]]$rate)      
  }    
  
  
  # Save those results, copied once for each iteration - those are the values at time 0, without
  resu<-cbind(expand.grid(Year=0,Sim=1:nSims),priors[[x]],EstPar=x,lcl=lcl,Mean=Mean,ucl=ucl)
  
  return(resu)})

# Name this list with the same name as master prior list
names(priors_timelapse)<-names(priors)

# Getting the best monitoring combination from Chapter 3
BestMonitoringOption<-(VOI_backup$finalConsequenceTable%>%pull(m))[1]

# Getting the best monitoring combination that uses all monitoring options
BestCompleteMonitoringOption<-(VOI_backup$finalConsequenceTable%>%
                                 filter(TrackedReleasedBirds>0,NoCameras>0)%>%
                                 pull(m))[1]
# Overrides to the 16th monitoring option
BestCompleteMonitoringOption=16

# Gets the monitoring options from Chapter 3
monitoringOptions<-VOI_backup$monitoringOptions

#------------------------------------------#
#--Run simulations and simulated datasets--#
#------------------------------------------#

i=1;t=1;pt<-0
s=1:nSims

# Data list to be imported into JAGS
dataGen.jags<-list(
  # Management Description
  n_alternatives=nrow(managementOptions),
  SoftFeed=managementOptions$SoftFeed,
  HardFeed=managementOptions$HardFeed,
  PoisonDrop=managementOptions$PoisonDrop,
  ReleasedAnimals=managementOptions$ReleasedAnimals,
  ReleaseSchedule=managementOptions$ReleaseSchedule,
  # Sample size of information collected
  fec.sample.size=monitoringOptions$NoCameras[BestCompleteMonitoringOption],
  surv.sample.size=monitoringOptions$TrackedReleasedBirds[BestCompleteMonitoringOption],
  seen=1,
  # Priors for fecundity 
  FledgNumberIntercept=as.numeric((priors_timelapse$reprodIntercept%>%filter(Year==0))[1,c("mean","sd")]),
  FledgNumberSoftFeedEffect=as.numeric((priors_timelapse$reprodBetaSoftFeed%>%filter(Year==0))[1,c("mean","sd")]),
  FledgNumberHardFeedEffect=as.numeric((priors_timelapse$reprodBetaHardFeed%>%filter(Year==0))[1,c("mean","sd")]),
  # Priors for survival
  survivalIntercept=as.numeric((priors_timelapse$survIntercept%>%filter(Year==0))[1,c("mean","sd")]),
  survivalAgeEffect=as.numeric((priors_timelapse$survAgeEffect%>%filter(Year==0))[1,c("mean","sd")]),
  survivalMastEffect=as.numeric((priors_timelapse$survMastEffect%>%filter(Year==0))[1,c("mean","sd")]),
  survivalPoisonDropEffect=as.numeric((priors_timelapse$survPoisonDropEffect%>%filter(Year==0))[1,c("mean","sd")]),
  survivalSoftFeedEffect=as.numeric((priors_timelapse$survBetaSoftFeed%>%filter(Year==0))[1,c("mean","sd")]),
  survivalHardFeedEffect=as.numeric((priors_timelapse$survBetaHardFeed%>%filter(Year==0))[1,c("mean","sd")]),
  survivalCostOfRelease=as.numeric((priors_timelapse$survCostRelease%>%filter(Year==0))[1,c("shape1","shape2")])
)

# The indexes of where in the Leslie matrix the zeroes would go
dataGen.jags$ZeroIndexes<-ZeroIndexes

# Starting population size, divided by age groups
dataGen.jags$N0<-c(0,25,25)  
# Number of year to simulate management
dataGen.jags$n_years<-nYears
# Number of year to project
dataGen.jags$n_years_proj<-50 ###

parameters.to.save=c(
  # Data generated
  "new.seen"
  ,"new.fledgings"
  #True values generated
  ,"ExpectedFecundity"
  ,"ExpectedSurvival"
  ,"survIntercept","survBetaSoftFeed","survBetaHardFeed"
  ,"survAgeEffect","survMastEffect","survPoisonDropEffect"
  ,"reprodIntercept","reprodBetaSoftFeed","reprodBetaHardFeed"
  ,"isMastYear","PropMastYears","MastProb"
  ,"Lambda"
  ,"TotalN"
  ,"Decline")


# First piece of JAGS code, describing the model (the relationships between variables)
jags.pt1<-readLines(paste0(".\\Chapter 4\\jags_model.txt",sep=""))
# Second piece of JAGS code, describing the operations done (in this case generating data)
jags.pt2<-readLines(paste0(".\\Chapter 4\\jags_data_generation_end.txt",sep=""))


cat(paste0(c(jags.pt1,jags.pt2),collapse="\n"),file=".\\Chapter 4\\data_gen.txt")


# Generated true values, and datasets to maatch, using the JAGS model just assembled
generatedData0<-jagsUI::jags(
  data=dataGen.jags,
  # inits=initial.list,
  model.file=".\\Chapter 4\\data_gen.txt",
  DIC=F,
  parameters.to.save=parameters.to.save,
  n.iter=max(nSims,2000),
  n.adapt=100,
  n.burnin=0,
  n.thin=1,
  n.chains=1,
  # n.chains=1,
  verbose=T,
  # seed="19191919",
  codaOnly = parameters.to.save[length(parameters.to.save)])

# Reshape fecundity data as a long format
generatedData0$sims.list$new.fledgings.melt<-generatedData0$sims.list$new.fledgings%>%
  reshape2::melt()%>%
  dplyr::rename(Sim=Var1,Sample=Var2,mgmt=Var3,Year=Var4)

# Reshape encounter history data as a long format
generatedData0$sims.list$new.seen.melt<-generatedData0$sims.list$new.seen%>%
  reshape2::melt()%>%
  dplyr::rename(Sim=Var1,Sample=Var2,Survey=Var3,mgmt=Var4,Year=Var5)
  

# This function estimates the hyper parameters of the priors defined in 'parDist' 
# In this case, we are just coming up with the with the hyper parameters for each
# parameter of interest in our model on the time step 0 (i.e., before analysing data)
hpars<-estimateBUGShyperPars(input = generatedData0,
                              parDist = c(survIntercept="norm",
                                          survBetaSoftFeed="norm",
                                          survAgeEffect="norm",
                                          survMastEffect="norm",
                                          survBetaHardFeed="norm",
                                          reprodIntercept="norm",
                                          reprodBetaSoftFeed="norm",
                                          reprodBetaHardFeed="norm",
                                          TimeToK="norm",
                                          ExpectedSurvival="beta",
                                          MastProb="beta",
                                          Lambda = "gamma",
                                          Decline = "norm",
                                          ExpectedFecundity="gamma"),
                              #start=start,
                              #fix.arg=fix.arg,
                              method="jagsUI")%>%
  dplyr::mutate(Year=0)%>%
  merge(data.frame(Year=0,Sim=1:nSims),by="Year")

# Add those values to the priors records object
for(p in unique(hpars$TruePar)){
  priors_timelapse[[p]]<-rbind.fill(priors_timelapse[[p]],hpars%>%filter(TruePar==p))%>%
    dplyr::mutate(range=ucl-lcl)%>%
    dplyr::group_by(EstPar)%>%
    # This extracts what's inside the brackets of the parameter name
    dplyr::mutate(ParIndex=gsub(x=empty2na(str_extract_all(EstPar, "\\[[^()]+\\]")[[1]]),"\\[|\\]",""))
}

##############################################################################
# Extracting the true values for each simulation (ie the predicted outcomes)##
##############################################################################

# Probability of a mast year
ProbMastYears.df<-reshape2::melt(generatedData0$sims.list$MastProb)%>%
  dplyr::rename(MastProb=value)%>%
  dplyr::mutate(Sim=1:n())

# Merge it with Lambda true value, and management/monitoring alternatives
# It also calculates the true scores under SMART analysis
trueValues<-reshape2::melt(generatedData0$sims.list$Lambda)%>%
  dplyr::rename(Sim=Var1,mgmt=Var2,Lambda=value)%>%
  merge(ProbMastYears.df,by="Sim")%>%
  merge(managementOptions,by="mgmt")%>%
  # Calculate the cost of management for each management scenario
  dplyr::mutate(CostMgmt=(SoftFeed*CostOfSoftSupp)+
                         (HardFeed*CostOfHardSupp)+
                         (ReleasedAnimals*ReleaseSchedule*1000)+
                         (PoisonDrop*30e3*MastProb))%>%
  dplyr::mutate(CostScore=ObjectiveScore(CostMgmt,minimize = T)*VOI_backup$Weights[2],
                BioScore=ObjectiveScore(Lambda,
                                        Max=VOI_backup$globalLimits$lambda$max,
                                        Min=VOI_backup$globalLimits$lambda$min)*VOI_backup$Weights[1],
                Score=CostScore+BioScore)%>%
  dplyr::group_by(Sim)%>%
  dplyr::mutate(BestScore=Score==max(Score),
                SecondBestScore = Score == sort(Score,decreasing = F)[2])



jags.pt1<-readLines(paste0(".\\Chapter 4\\jags_model.txt",sep=""))
jags.pt2<-readLines(paste0(".\\Chapter 4\\jags_data_interpretation_end.txt",sep=""))
cat(paste0(c(jags.pt1,jags.pt2),collapse="\n"),file=".\\Chapter 4\\data_read.txt")

# Create a table to hold the consequence tables throughout the simulate study
# (so we can follow over each time step)
ConsequenceTableHistory<-trueValues%>%
  dplyr::group_by(mgmt)%>%
  dplyr::summarise(Lambda=mean(Lambda),Score=mean(Score),CostMgmt=mean(CostMgmt))%>%
  dplyr::select(mgmt,Lambda,Score,CostMgmt)%>%
  dplyr::mutate(Year=0)%>%
  dplyr::mutate(BestScore=Score==max(Score),
                SecondBestScore = Score == sort(Score,decreasing = F)[2])%>%
  merge(expand.grid(mgmt=managementOptions$mgmt,Sim=1:nSims),by="mgmt")

# Empty list to hold simulated values  
jags.list<-list()

# Mark when the time consuming part of the analysis starts
start.<-now()

# When to print progress
print.crit<-min(25,nSims)

updatePriors<-T # whether to use updated prior (the posterior from the previous time step) or original ones
useFullData<-T  #whether to use data from previous years or just the just collected data
flushJAGS<-T # Whether to flush the JAGS results to save memory

data.jags.list<-list()

# Marker on how many iterations where conducted
pt<-0


# for each year 
for(t in 1:nYears){
    # and each iteration
    for(i in 1:nSims){
      
      # Add iteration counter
      pt<-pt+1
      
      # Return Decisions historically taken
      decLoop = rbind.fill(ConsequenceTableHistory%>%filter(Sim==i,Year<t-1,BestScore),
                           ConsequenceTableHistory%>%filter(Sim==i,Year==t-1,BestScore))
      
      # Extract the decision deemed as best
      currentDecision = decLoop%>%filter(Year==max(Year))%>%pull(mgmt)%>%as.numeric()
      
      # If we used the full data, the indexes to filter is 1 until current year in the loop
      # Otherwise, just the current year
      if(useFullData){td=1:t}else{td=t}
      
      # Extract and and rearrange fecundity data
      FledgingDataLoop<-lapply(1:nrow(decLoop),function(dd){
        decision=decLoop[dd,"mgmt"]
        year = decLoop[dd,"Year"]+1
        resu<-generatedData0$sims.list$new.fledgings.melt%>%
          dplyr::filter(Sim==i,mgmt==decision,Year==year)
        return(resu)
      })%>%
        plyr::rbind.fill()%>%
        filter(Year%in%td)
      
      # Extract and and rearrange survival data (encounter history)
      SurvDataLoop<-lapply(1:nrow(decLoop),function(dd){
        decision=decLoop[dd,"mgmt"]
        year = decLoop[dd,"Year"]+1
        resu<-generatedData0$sims.list$new.seen.melt%>%
          dplyr::filter(Sim==i,mgmt==decision,Year==year)%>%
        return(resu)
      })%>%
        plyr::rbind.fill()%>%
        dplyr::group_by(Survey)%>%
        dplyr::arrange(Year,Sample)%>%
        dplyr::mutate(Sample2=1:n())%>%
        dplyr::filter(Year%in%td)
      
      # Reshape encounter history as a matrix
      SeenLoop<-SurvDataLoop%>%
        dplyr::select(Sample2,Survey,value)%>%
        reshape2::acast(Sample2~Survey)
      
      # The management that took place (used as an index)
      PhiMgmt=SurvDataLoop%>%filter(Survey==1)%>%pull(mgmt)
      
      # Index of which prior to use. either the first one, or the last generated
      tt<-ifelse(updatePriors,t-1,0)
      
      # Create a data list to be inputed into JAGS
      data.jags<-list(
        ### Management Alternatives ####
        n_alternatives=nrow(managementOptions),
        SoftFeed=managementOptions$SoftFeed,
        HardFeed=managementOptions$HardFeed,
        PoisonDrop=managementOptions$PoisonDrop,
        ReleasedAnimals=managementOptions$ReleasedAnimals,
        ReleaseSchedule=managementOptions$ReleaseSchedule,
        
        ##############
        ### DATA ####
        #############
        fledgings=FledgingDataLoop$value,
        # Binomial identifier of which management was conducted 
        FecSoftFeed = managementOptions$SoftFeed[FledgingDataLoop$mgmt],
        FecHardFeed = managementOptions$HardFeed[FledgingDataLoop$mgmt],
        FecMast=generatedData0$sims.list$isMastYear[i,FledgingDataLoop$Year],
        # Which year was 
        FecTime=FledgingDataLoop$Year,
        ######################
        ### Survival data  ##
        ######################
        seen=SeenLoop,
        # Binomial identifiers of which management was conducted 
        PhiSoftFeed = managementOptions$SoftFeed[PhiMgmt],
        PhiHardFeed = managementOptions$HardFeed[PhiMgmt],
        PhiPoison = managementOptions$PoisonDrop[PhiMgmt],
        PhiMast = generatedData0$sims.list$isMastYear[i,SurvDataLoop%>%filter(Survey==1)%>%
                                                        pull(Year)],
        PhiTime=SurvDataLoop%>%filter(Survey==1)%>%pull(Year),

    
        ##########
        # Priors #
        ##########
        FledgNumberIntercept=as.numeric((priors_timelapse$reprodIntercept%>%filter(Year==tt,Sim==i))[1,c("mean","sd")]),
        FledgNumberSoftFeedEffect=as.numeric((priors_timelapse$reprodBetaSoftFeed%>%filter(Year==tt,Sim==i))[1,c("mean","sd")]),
        FledgNumberHardFeedEffect=as.numeric((priors_timelapse$reprodBetaHardFeed%>%filter(Year==tt,Sim==i))[1,c("mean","sd")]),
        # Priors
        survivalIntercept=as.numeric((priors_timelapse$survIntercept%>%filter(Year==tt,Sim==i))[1,c("mean","sd")]),
        survivalAgeEffect=as.numeric((priors_timelapse$survAgeEffect%>%filter(Year==tt,Sim==i))[1,c("mean","sd")]),
        survivalMastEffect=as.numeric((priors_timelapse$survMastEffect%>%filter(Year==tt,Sim==i))[1,c("mean","sd")]),
        survivalPoisonDropEffect=as.numeric((priors_timelapse$survPoisonDropEffect%>%filter(Year==tt,Sim==i))[1,c("mean","sd")]),
        survivalSoftFeedEffect=as.numeric((priors_timelapse$survBetaSoftFeed%>%filter(Year==tt,Sim==i))[1,c("mean","sd")]),
        survivalHardFeedEffect=as.numeric((priors_timelapse$survBetaHardFeed%>%filter(Year==0))[1,c("mean","sd")]),
        survivalCostOfRelease=as.numeric((priors_timelapse$survCostRelease%>%filter(Year==0))[1,c("shape1","shape2")])
      )
      
      ########################################
      ##Operations on encounter history data##    
      ########################################
      
      # First column of encounter history is always presence (time of banding)
      data.jags$seen[,1]<-1
      # Encounter history of aliveness is the same of detection
      data.jags$alive<-data.jags$seen
      # Except the zeroes are now missing values 
      data.jags$alive[data.jags$alive==0]<-NA
      # Unless they are followed by a 1, in which case we can know the individual was indeed alive 
      for(ii in 1:nrow(data.jags$alive)){
        last.seen<-max(which(data.jags$alive[ii,]==1))
        data.jags$alive[ii,1:last.seen]<-1        
        
      }
      
      data.jags$ZeroIndexes<-ZeroIndexes
      data.jags$N0<-c(0,25,25) ###
      data.jags$n_years<-nYears ###
      data.jags$n_years_proj<-50 ###
      
      jagsLoop<-jagsUI::jags(
        data=data.jags,
        # inits=initial.list,
        model.file=".\\Chapter 4\\data_read.txt",
        parameters.to.save=c(
          "ExpectedFecundity",
          "ExpectedSurvival",
          "survIntercept","survBetaSoftFeed","survBetaHardFeed"
          ,"survAgeEffect","survMastEffect","survPoisonDropEffect"
          ,"reprodIntercept","reprodBetaSoftFeed","reprodBetaHardFeed"
          ,"Lambda"
          ,"Decline","MastProb"
          ,"finalN"
          ,"IsAtCarryingCapacity"
        ),
        n.iter=3e3,
        n.adapt=200,
        n.burnin=1e3,
        n.thin=1,
        n.chains=2,
        # n.chains=1,
        verbose=F,
        codaOnly = T)
        
        # Derived parameters (possible with R only)
        # Time to Reach K 
        jagsLoop$sims.list$TimeToK<-apply(
          
        jagsLoop$sims.list$IsAtCarryingCapacity,
          # For each iteration and management  
          c(1,3),
          # Calculate the first value where N is at Carrying capacity
          function(x){
            y<-min(which(x==1))
            # If it never happens, this returns an infinite value. 
            # Replace it with an arbitraly large one
            y[y==Inf]<-999
            return(y)
          })
        
        # Flag which iteration and year this belongs too
        jagsLoop$Sim=i
        jagsLoop$Time=t
        
        # Save JAGS results to a list
        jags.list[[pt]]<-jagsLoop
        
        # Save the data used in a list as well
        data.jags.list[[pt]]<-data.jags
        data.jags.list[[pt]]$Sim<-i
        data.jags.list[[pt]]$Time<-t
        
        # This optionally gets rid of the JAGS results
		    if(flushJAGS){jags.list[[pt]]<-NULL}
        #jags.list[[pt]]$samples<-MigrateSimsListToSamples(jags.list[[pt]])
        
        jagsLoop$samples<-MigrateSimsListToSamples(jagsLoop)
        
        #fix.arg<-list(survBetaSoftFeed=list(a=0))
        #start<-list(survBetaSoftFeed=list(mean=1,sd=1))
                        
        # start<-list(survBetaSoftFeed=list(a=0,mean=1,sd=1))
        # start$survBetaHardFeed<-
        #   start$reprodBetaSoftFeed<-
        #   start$reprodBetaHardFeed<-start$survBetaSoftFeed
          
        hpars<-estimateBUGShyperPars(input = jagsLoop,
                            parDist = c(survIntercept="norm",
                                        survBetaSoftFeed="norm",
                                        survAgeEffect="norm",
                                        survMastEffect="norm",
                                        survPoisonDropEffect="norm",
                                        survBetaHardFeed="norm",
                                        reprodIntercept="norm",
                                        reprodBetaSoftFeed="norm",
                                        reprodBetaHardFeed="norm",
                                        TimeToK="norm",
                                        ExpectedSurvival="beta",
                                        MastProb="beta",
                                        Lambda = "gamma",
                                        Decline = "norm",
                                        ExpectedFecundity="gamma"),
                            #start=start,
                            #fix.arg=fix.arg,
                            method="jagsUI")%>%
        dplyr::mutate(Year=t,Sim=i)
        
        for(p in unique(hpars$TruePar)){
          priors_timelapse[[p]]<-rbind.fill(priors_timelapse[[p]],hpars%>%filter(TruePar==p))%>%
            dplyr::mutate(range=ucl-lcl)%>%
            dplyr::group_by(EstPar)%>%
            dplyr::mutate(ParIndex=gsub(x=empty2na(str_extract_all(EstPar, "\\[[^()]+\\]")[[1]]),"\\[|\\]",""))
                          }
        
      #  
      if(i%%print.crit==0){print(paste("Sim: ",i,", Time: ",t,sep=""))}
    }
    
    # Making the decision for all simulations, at a given time step
    # Consequence Table
    
    ConseqTableLoop<-priors_timelapse$Lambda%>%filter(Year==t)%>%
      dplyr::group_by(Sim,ParIndex)%>%
      dplyr::summarise(Lambda=mean(Mean))%>%
      dplyr::mutate(mgmt=ParIndex,ParIndex=NULL)%>%
      
      merge(priors_timelapse$MastProb%>%filter(Year==t)%>%
              dplyr::rename(MastProb=Mean)%>%
              dplyr::select(Sim,Year,MastProb),by=c("Sim"))%>%
      
      merge(managementOptions,by="mgmt")%>%
      
      dplyr::mutate(CostMgmt=(SoftFeed*CostOfSoftSupp)+
                      (HardFeed*CostOfHardSupp)+
                      (ReleasedAnimals*ReleaseSchedule*1000)+
                      (PoisonDrop*30e3*MastProb))%>%
      
      dplyr::mutate(CostScore=ObjectiveScore(CostMgmt,minimize = T,Max = max(trueValues$CostMgmt))*VOI_backup$Weights[2],
                    BioScore=ObjectiveScore(Lambda,
                                            Max=VOI_backup$globalLimits$lambda$max,
                                            Min=VOI_backup$globalLimits$lambda$min)*VOI_backup$Weights[1],
                    Score=CostScore+BioScore)%>%
      
      dplyr::group_by(Sim)%>%
      dplyr::mutate(BestScore=Score==max(Score))%>%
      dplyr::mutate(Year=t)

    ConsequenceTableHistory<-rbind.fill(ConsequenceTableHistory,ConseqTableLoop)
    
    
  }
end.<-now()

simsRan=ConsequenceTableHistory%>%filter(Year==1)%>%pull(Sim)%>%max
# Pushbullet Notification
if(1==1){RPushbullet::pbPost(title="JAGS Run Finished",
                             body=paste("Process finished in",
                                        abs(round(difftime(start.,end.,units="mins"),2)),
                                        "minutes. Ran ",
                                        simsRan,
                                        " simulations over ",
                                        ConsequenceTableHistory%>%pull(Year)%>%max,
                                        " time steps."),
                             type = "note",apikey = "o.NfaWdVVPxnLb8DY2Ab27rregtjC2O23M",email ="caio.kenup@gmail.com")
}

priors_timelapse$survMastEffect%>%filter(Sim==2)%>%arrange(Sim,Year)
priors_timelapse$survAgeEffect%>%filter(Sim<=simsRan)%>%arrange(Sim,Year)


alpha.level<-.5
mgmtPalette<-c("1"="#264653","2"="#e9c46a","3"="#e76f51")
scenarioPalette<-c("Perfect Information"="purple","Uncertainty"="gold")
mgmtMap<-c("No Feeding","Low Intensity Feeding","High Intensity Feeding")
SizeAdj=1.1

Fig.lambda.list<-list()
Score.fig.list<-list()

Fig.lambda.list<-list()
Score.fig.list<-list()
Fig.survivalBetas.list<-list()

Ch04_Example<-list(priors_timelapse=priors_timelapse,
                   VOI_backup=VOI_backup,
                   ConsequenceTableHistory=ConsequenceTableHistory,
                   trueValues=trueValues,
                   dataGen.jags=dataGen.jags,
                   managementOptions=managementOptions,
                   generatedData0=generatedData0)
                   # jags.list=jags.list)

save(Ch04_Example,file=paste0(
       "C:\\Users\\",username,"\\Dropbox\\03-Work\\01-Science\\02-PhD\\02-Analysis\\Chapter 4\\AM Example v3.RData",sep=""))


load(file = paste0(
       "C:\\Users\\",username,"\\Dropbox\\03-Work\\01-Science\\02-PhD\\02-Analysis\\Chapter 4\\AM Example v3.RData",sep=""))
# Ch04_Example0<-Ch04_Example


# Pushbullet Notification
if(1==1){RPushbullet::pbPost(title="JAGS Saving Finished",
                             body=paste("Info Saved!"),
                             type = "note",apikey = "o.NfaWdVVPxnLb8DY2Ab27rregtjC2O23M",email ="caio.kenup@gmail.com")
}


prow<-plot_grid(
  AM_fig1 + theme(legend.position="none"),
  AM_fig2 + theme(legend.position="none"))
#   CostMgmt + theme(legend.position="none")+xlab(NULL),
#   ScoreUnderUncertainty + theme(legend.position="none"),
#   align="v",
#   ncol=1,labels = c("A", "B","C"),label_x=0.18,label_y=.98)
legend <- get_legend(
  AM_fig1 +
    guides(fill=guide_legend(title="Alternatives",ncol=1,override.aes = list(alpha=1)))+
    theme(legend.position = "bottom")
)
p<-cowplot::plot_grid(prow, legend, ncol = 1, rel_widths = c(1,.45))



PriorExpect<-trueValues%>%
  dplyr::mutate(PoisonDropLab=ifelse(PoisonDrop,"Poison Drop","No Poison Drop"))%>%
  dplyr::mutate(FeedingLab=ifelse(SoftFeed,"Low Intensity Feeding","No Feeding"))%>%
  dplyr::mutate(FeedingLab=ifelse(HardFeed,"High Intensity Feeding",FeedingLab))
  

ggplot(PriorExpect%>%dplyr::mutate(mgmt=factor(mgmt)),aes(x=Lambda,fill=FeedingLab))+
  geom_density(mapping = aes(group=FeedingLab,fill=FeedingLab),color=NA,alpha=.5)+
  geom_vline(data=PriorExpect%>%
               dplyr::group_by(PoisonDropLab,FeedingLab)%>%
               dplyr::summarise(xintercept=mean(Lambda)),
             aes(xintercept=xintercept,color=FeedingLab),lty="dotted")+
  facet_wrap(~PoisonDropLab,ncol=1)


ggplot(PriorExpect%>%dplyr::mutate(mgmt=factor(mgmt)),aes(x=Score,fill=FeedingLab))+
  geom_density(mapping = aes(group=FeedingLab,fill=FeedingLab),color=NA,alpha=.5)+
  geom_vline(data=PriorExpect%>%
               dplyr::group_by(PoisonDropLab,FeedingLab)%>%
               dplyr::summarise(xintercept=mean(Score)),
             aes(xintercept=xintercept,color=FeedingLab),lty="dotted")+
  facet_wrap(~PoisonDropLab,ncol=1,scales="free_y")

SimExample<-1


plot_data<-rbind.fill(priors_timelapse$Lambda%>%mutate(mgmt=ParIndex))
trueData<-data.frame(TrueLambda=generatedData0$Lambda[SimExample,],mgmt=1:length(generatedData0$Lambda[SimExample,]))
  
AM_fig1<-ggplot(plot_data%>%
                 filter(Sim==SimExample),aes(x=Year,y=Mean,color=mgmt))+
  # geom_smooth(method = "loess",size=.45)+
  geom_line(size=.45)+
  # geom_smooth(method = "loess",aes(x=Year,y=ucl),linetype="dashed",size=.33)+
  # geom_smooth(method = "loess",aes(x=Year,y=lcl),linetype="dashed",size=.33)+
  geom_ribbon(aes(ymin=lcl,ymax=ucl,x=Year,fill=mgmt),alpha=.33)+
  geom_hline(data=trueData,mapping = aes(yintercept=TrueLambda),linetype="dotted",col="red")+
  facet_wrap(~mgmt,scales = "fixed")+
  ylab("")
AM_fig1


score_data<-ConsequenceTableHistory%>%filter(Sim==SimExample)
trueScores<-trueValues%>%filter(Sim==SimExample)%>%mutate(mgmt=as.factor(mgmt))
AM_fig2<-ggplot()+
  # geom_smooth(method = "loess",size=.45)+
  # 
  # Plot the best decisions taking at all points in time
  geom_point(data=score_data%>%
               filter(BestScore),mapping=aes(x=Year,y=Score,color=mgmt))+
  geom_line(data=score_data%>%
              filter(BestScore),mapping=aes(x=Year,y=Score),color="grey75")+
  # Plot the other estimated scores of the decisions
  geom_point(data=score_data,mapping=aes(x=Year,y=Score,color=mgmt))+
  
  # geom_smooth(method = "loess",aes(x=Year,y=ucl),linetype="dashed",size=.33)+
  # geom_smooth(method = "loess",aes(x=Year,y=lcl),linetype="dashed",size=.33)+
  # geom_ribbon(aes(ymin=lcl,ymax=ucl,x=Year,fill=mgmt),alpha=.33)+
  geom_hline(data=trueScores,mapping = aes(yintercept=Score,color=mgmt),linetype="dotted",size=1)+
  # facet_wrap(~mgmt,scales = "fixed")+
  ylab("")
AM_fig2


expand.grid(mgmt=1:managementOptions$m,Year=1:nYear,Sim=1:nSims)
