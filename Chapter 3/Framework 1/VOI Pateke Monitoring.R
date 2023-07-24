###############################################
######Value of Information Pateke Monitoring###
###############################################

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

devtools::source_url("https://raw.githubusercontent.com/KenupCF/ThesisPHD/main/GeneralFunctions/ObjectiveScore.R")
devtools::source_url("https://raw.githubusercontent.com/KenupCF/ThesisPHD/main/GeneralFunctions/estimateBUGShyperPars.R")

devtools::source_url("https://raw.githubusercontent.com/KenupCF/Kenup-Repo/master/Quick%20Functions.R")



# Object with the proper link and ilink funcitons for each distribution
linkFUN<-list("gamma"=log,
              "pois"=log,
              "beta"=logit,
              "normal"=identity)
ilinkFUN<-list("gamma"=exp,
               "pois"=exp,
               "beta"=inv.logit,
               "normal"=identity)


####
#### Part 0
####

# Number of simulations to run
nSims=4
# Large numbers will take quite sometime...
# nSims=1000

# Number of years to simulate for
nYears=1


###########################################
### Part I ################################
### Defining the decision scenario ########
###########################################

#### Cost ###
###Price of stuff

CostOfCameras<-250 #in NZD
CostOfRadios<-500 # in NZD

CostOfSoftSupp<-10e3 #in NZD
CostOfHardSupp<-20e3 #in NZD

# Dataframe holding possible monitoring options
monitoringOptions<-expand.grid(
  
  NoCameras=c(0
              ,5
              ,10
              ,20
  ),
  
  TrackedReleasedBirds=c(0
                         ,5
                         ,10
                         ,20
  ),
  
  RadioSurveyLength=c(0,6))

# Remove illogical monitoring options
monitoringOptions<-monitoringOptions%>%
  dplyr::filter(!(RadioSurveyLength==0 & TrackedReleasedBirds>0))%>%
  dplyr::filter(!(TrackedReleasedBirds==0 & RadioSurveyLength>0))
monitorCols<-colnames(monitoringOptions)


# Calculate  cost of monitoring option
monitoringOptions<-monitoringOptions%>%
  dplyr::mutate(CostMonitoring=(NoCameras*CostOfCameras)+
                  (TrackedReleasedBirds*CostOfRadios),
                # 'm' is the id defining each monitoring scenarion
                m=1:n())


# Extract which of the monitoring options is full uncertainty
# i.e., no new information
Uncertainty_index<-which(apply(
  monitoringOptions[,monitorCols],1,function(x){sum(x)==0}))

# This is an integer
Uncertainty_index

# Creating management options
managementOptions<-expand.grid(
  SoftFeed=0:1,
  HardFeed=0:1,
  # For the  sake of simplicty no further releases were considered
  ReleasedAnimals=c(0,0),
  ReleaseSchedule=0:2)%>%
  # Removing illogical management alternatives
  # Both low and high intensity feeding
  dplyr::filter(!(SoftFeed == 1 & HardFeed == 1))%>%
  # Schedule of releases without animals, and vice versa
  dplyr::filter(!(ReleasedAnimals == 0 & ReleaseSchedule > 0))%>%
  dplyr::filter(!(ReleasedAnimals > 0 & ReleaseSchedule == 0))

# Remove any duplicated management options
managementOptions<-managementOptions%>%
  dplyr::filter(!duplicated(managementOptions))

# Save name of columns on management data.frame
managementCols<-colnames(managementOptions)

# Dataframe holding possible management options and their cost
managementOptions<-managementOptions%>%
  dplyr::mutate(CostMgmt=(SoftFeed*CostOfSoftSupp)+
                  (HardFeed*CostOfHardSupp)+
                  (ReleasedAnimals*ReleaseSchedule*1000),
                mgmt=1:n())
Baseline_index<-which(apply(
  managementOptions[,managementCols],1,function(x){sum(x)==0}))	

############
###Priors###
############

#Defining arbitrary priors
priors<-list(
  reprodIntercept=data.frame(mean=log(1.1),sd=.75,dist="norm"),
  reprodBetaSoftFeed=data.frame(mean=0.27,sd=.45,dist="norm"),
  reprodBetaHardFeed=data.frame(mean=0.54,sd=.55,dist="norm"),
  survIntercept=data.frame(mean=logit(.67),sd=.43,dist="norm"),	
  survAgeEffect=data.frame(mean=0,sd=1e-3,dist="norm"),	
  survBetaSoftFeed=data.frame(mean=0.50,sd=.54,dist="norm"),	
  survBetaHardFeed=data.frame(mean=0.75,sd=.54,dist="norm"),
  survCostRelease=data.frame(shape1=90,shape2=10,dist="beta")
)
# Define number of age classes for species
noAgeClasses<-3


# Calculate where the 0 values of a Leslie matrix would go, for a given number of age classes
# (i.e., impossible transitions)
{
  L<-matrix(0,nrow=noAgeClasses,ncol=noAgeClasses)
  
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

# Objects to hold results
monitoringResults<-list()
jagsResults<-list()

# Start iterations
start..<-now()

# Empty objects to hold results
no_survivors.df.General<-data.frame()
n_fledglings.General<-data.frame()
generatedData<-list()
trueValues<-list()

# Two JAGS text files
jags.pt1<-readLines(".\\jags_model.txt") #one describing the relationships between variables (i.e. the model)
jags.pt2<-readLines(".\\jags_data_generation_end.txt") # one describing the data that will be generated from model
# Bind two files together
cat(paste0(c(jags.pt1,jags.pt2),collapse="\n"),file="data_gen.txt")


### For each monitoring scenario
for(m in 1:nrow(monitoringOptions)){
  
  # Object with data to input into JAGS
  dataGen.bugs<-list(
    # Management Description
    n_alternatives=nrow(managementOptions),
    SoftFeed=managementOptions$SoftFeed,
    HardFeed=managementOptions$HardFeed,
    ReleasedAnimals=managementOptions$ReleasedAnimals,
    ReleaseSchedule=managementOptions$ReleaseSchedule,
    
    fec.sample.size=monitoringOptions$NoCameras[m],
    surv.sample.size=monitoringOptions$TrackedReleasedBirds[m],
    seen=1,
    # Priors
    FledgNumberIntercept=as.numeric((priors$reprodIntercept)[1,c("mean","sd")]),
    FledgNumberSoftFeedEffect=as.numeric((priors$reprodBetaSoftFeed)[1,c("mean","sd")]),
    FledgNumberHardFeedEffect=as.numeric((priors$reprodBetaHardFeed)[1,c("mean","sd")]),
    # Priors
    survivalIntercept=as.numeric((priors$survIntercept)[1,c("mean","sd")]),
    survivalAgeEffect=as.numeric((priors$survAgeEffec)[1,c("mean","sd")]),
    survivalSoftFeedEffect=as.numeric((priors$survBetaSoftFeed)[1,c("mean","sd")]),
    survivalHardFeedEffect=as.numeric((priors$survBetaHardFeed)[1,c("mean","sd")]),
    survivalCostOfRelease=as.numeric((priors$survCostRelease)[1,c("shape1","shape2")]) )
  
  
  dataGen.bugs$ZeroIndexes<-ZeroIndexes
  # Starting population (split by age class)
  dataGen.bugs$N0<-c(0,25,25)  
  
  dataGen.bugs$n_years<-nYears
  
  # Number of years to project population growth
  dataGen.bugs$n_years_proj<-50 ###
  
  runToyMCMC=T
  
  # Run JAGS to generate simulated monitoring data
  generatedData[[m]]<-jagsUI::jags(
    data=dataGen.bugs,
    # inits=initial.list,
    model.file="data_gen.txt",
    parameters.to.save=c(
      "new.fledgings"
      ,"new.seen"
      ,"isMastYear"
      ,"Lambda"
      ,"Decline"
      ,"ExpectedSurvival"
      ,"ExpectedFecundity"
    ),
    n.iter=nSims,
    n.adapt=100,
    n.burnin=0,
    n.thin=1,
    n.chains=1,
    # seed=26051991,
    verbose=F,
    codaOnly = T)
  
  # Fixing instances of no monitoring (all data values are NA if there is no monitoring of a given quantity)
  if(dataGen.bugs$fec.sample.size==0){generatedData[[m]]$sims.list$new.fledgings = generatedData[[m]]$sims.list$new.fledgings*NA}  
  if(dataGen.bugs$surv.sample.size==0){generatedData[[m]]$sims.list$new.seen = generatedData[[m]]$sims.list$new.seen*NA}  
  
  
  #Reshape fecundity 
  generatedData[[m]]$sims.list$new.fledgings.melt<-generatedData[[m]]$sims.list$new.fledgings%>%
    reshape2::melt()%>%
    dplyr::rename(Sim=Var1,Sample=Var2,mgmt=Var3,Year=Var4)
  
  generatedData[[m]]$sims.list$new.seen.melt<-generatedData[[m]]$sims.list$new.seen%>%
    reshape2::melt()%>%
    dplyr::rename(Sim=Var1,Sample=Var2,Survey=Var3,mgmt=Var4,Year=Var5)
  
  #df containing true values of lambda for each simulation
  trLambda<-generatedData[[m]]$sims.list$Lambda%>%
    reshape2::melt()%>%
    dplyr::rename(Sim=Var1,mgmt=Var2,Lambda=value)
  
  # Same for values of whether population was declining
  trDecline<-generatedData[[m]]$sims.list$Decline%>%
    reshape2::melt()%>%
    dplyr::rename(Sim=Var1,mgmt=Var2,Decline=value)
  
  # Same for survival
  trSurvival<-generatedData[[m]]$sims.list$ExpectedSurvival%>%
    reshape2::melt()%>%
    dplyr::rename(Sim=Var1,Age=Var2,mgmt=Var3,Survival=value)%>%
    dplyr::mutate(Age = paste("Age_",Age,"_Survival",sep=""))%>%
    tidyr::pivot_wider(names_from=Age,values_from=Survival)
  
  # Same for fecundity
  trFecundity<-generatedData[[m]]$sims.list$ExpectedFecundity%>%
    reshape2::melt()%>%
    dplyr::rename(Sim=Var1,Age=Var2,mgmt=Var3,Fecundity=value)%>%
    dplyr::mutate(Age = paste("Age_",Age,"_Fecundity",sep=""))%>%
    tidyr::pivot_wider(names_from=Age,values_from=Fecundity)
  
  # Merge all those values into a single data.frame, and save to the corresponding monitoring scenario
  trueValues[[m]]<-merge(trLambda,trDecline,by=c("Sim","mgmt"))%>%
    merge(trSurvival,by=c("Sim","mgmt"))%>%
    merge(trFecundity,by=c("Sim","mgmt"))%>%
    dplyr::mutate(m=m)
}




## Calculate predicted outcomes for management scenarios

PredictedOutcomesFull<-plyr::rbind.fill(trueValues)

str(PredictedOutcomesFull) # full data.frame showing predicted outcome for each iteration of each scenario

# Calculate limits on predicted outcomes as the 95% range
globalLimits<-list(lambda=list(min=quantile(PredictedOutcomesFull$Lambda,.025),
                               max=quantile(PredictedOutcomesFull$Lambda,0.975)))


### Objective weights ###
Weights<-c("biological"=.82,"monetary"=.18)

PredictedOutcomesSummary<-PredictedOutcomesFull%>%
  dplyr::group_by(mgmt)%>%
  dplyr::summarise(Lambda=mean(Lambda),Decline=mean(Decline))%>%
  merge(managementOptions)%>%
  dplyr::mutate(CostScore=ObjectiveScore(CostMgmt,minimize = T)*Weights[2],
                BioScore=ObjectiveScore(Lambda,
                                        Max=globalLimits$lambda$max,
                                        Min=globalLimits$lambda$min)*Weights[1],
                BioScore2=ObjectiveScore(Decline,Min = 0,Max=0.5,minimize = T)*Weights[1],
                Score=CostScore+BioScore,
                Score2=CostScore+BioScore2)

decisionUnderUncertainty<-PredictedOutcomesSummary%>%
  filter(Score==max(Score))%>%
  pull(mgmt)


# jags.list<-list()
# jags.pred<-list()

totalSims<-expand.grid(m=monitoringOptions$m,
                       i=1:nrow(generatedData[[1]]$samples[[1]]))



jags.list<-list()
hPars<-list()
# For each index of the parameter space 
ss<-1:nrow(totalSims)
pt<-0


jags.pt1<-readLines(".\\jags_model.txt")
jags.pt2<-readLines(".\\jags_data_interpretation_end.txt")
cat(paste0(c(jags.pt1,jags.pt2),collapse="\n"),file="data_read.txt")


print.crit<-50

pb <- winProgressBar(title = "Running", # Window title
                     label = "Percentage completed", # Window label
                     min = 0,      # Minimum value of the bar
                     max = length(ss), # Maximum value of the bar
                     initial = 0,  # Initial value of the bar
                     width = 300L) # Width of the window 
flushJags=FALSE

for(s in ss){
  
  pt<-pt+1
  i=totalSims$i[s]
  m=totalSims$m[s]
  
  FledgingDataLoop<-generatedData[[m]]$sims.list$new.fledgings.melt%>%
    dplyr::filter(Sim==i,mgmt==decisionUnderUncertainty)
  
  SurvDataLoop<-generatedData[[m]]$sims.list$new.seen.melt%>%
    dplyr::filter(Sim==i,mgmt==decisionUnderUncertainty)%>%
    dplyr::group_by(Survey)%>%
    dplyr::arrange(Year,Sample)%>%
    dplyr::mutate(Sample2=1:n())
  
  SeenLoop<-SurvDataLoop%>%
    dplyr::select(Sample2,Survey,value)%>%
    acast(Sample2~Survey)
  
  PhiMgmt<-SurvDataLoop%>%filter(Survey==1)%>%pull(mgmt)
  # Extracting Data from Previous Collected Dataset
  data.bugs<-list(
    ### Management Alternatives ####
    n_alternatives=nrow(managementOptions),
    SoftFeed=managementOptions$SoftFeed,
    HardFeed=managementOptions$HardFeed,
    ReleasedAnimals=managementOptions$ReleasedAnimals,
    ReleaseSchedule=managementOptions$ReleaseSchedule,
    
    ##############
    ### DATA ####
    #############
    fledgings=FledgingDataLoop$value,
    FecSoftFeed = managementOptions$SoftFeed[FledgingDataLoop$mgmt],
    FecHardFeed = managementOptions$HardFeed[FledgingDataLoop$mgmt],
    FecTime=FledgingDataLoop$Year,
    
    ###############################
    ### Survival data and priors ##
    ###############################
    seen=SeenLoop,
    PhiSoftFeed = managementOptions$SoftFeed[PhiMgmt],
    PhiHardFeed = managementOptions$HardFeed[PhiMgmt],
    PhiTime=SurvDataLoop%>%filter(Survey==1)%>%pull(Year),

    # Priors
    FledgNumberIntercept=as.numeric((priors$reprodIntercept)[1,c("mean","sd")]),
    FledgNumberSoftFeedEffect=as.numeric((priors$reprodBetaSoftFeed)[1,c("mean","sd")]),
    FledgNumberHardFeedEffect=as.numeric((priors$reprodBetaHardFeed)[1,c("mean","sd")]),
    # Priors
    survivalIntercept=as.numeric((priors$survIntercept)[1,c("mean","sd")]),
    survivalAgeEffect=as.numeric((priors$survAgeEffec)[1,c("mean","sd")]),
    survivalMastEffect=as.numeric((priors$survMastEffect)[1,c("mean","sd")]),
    survivalSoftFeedEffect=as.numeric((priors$survBetaSoftFeed)[1,c("mean","sd")]),
    survivalHardFeedEffect=as.numeric((priors$survBetaHardFeed)[1,c("mean","sd")]),
    survivalCostOfRelease=as.numeric((priors$survCostRelease)[1,c("shape1","shape2")]) )
  
  data.bugs$seen[,1]<-1
  data.bugs$alive<-data.bugs$seen
  data.bugs$alive[data.bugs$alive==0]<-NA
  data.bugs$ZeroIndexes<-ZeroIndexes
  data.bugs$N0<-c(0,25,25) ###
  data.bugs$n_years<-nYears ###
  data.bugs$n_years_proj<-50 ###
  
  # Run model (not a lot of chains)
  jags.list[[pt]]<-jagsUI::jags(
    data=data.bugs,
    # inits=initial.list,
    model.file="data_read.txt",
    parameters.to.save=c(
      "Lambda"
      ,"Decline"
    ),
    n.iter=ifelse(runToyMCMC,100,3e3),
    n.adapt=ifelse(runToyMCMC,5,200),
    n.burnin=ifelse(runToyMCMC,5,1e3),
    n.thin=ifelse(runToyMCMC,1,1),
    n.chains=ifelse(runToyMCMC,1,2),
    # n.chains=1,
    verbose=F,
    codaOnly = T)
  
  jags.list[[pt]]$m<-m
  jags.list[[pt]]$i<-i
  jags.list[[pt]]$model <- NULL
  
  hPars[[pt]]<-estimateBUGShyperPars(input = jags.list[[pt]],
                                     parDist = c(survIntercept="norm",
                                                 survBetaSoftFeed="norm",
                                                 survAgeEffect="norm",
                                                 survBetaHardFeed="norm",
                                                 reprodIntercept="norm",
                                                 reprodBetaSoftFeed="norm",
                                                 reprodBetaHardFeed="norm",
                                                 ExpectedSurvival="beta",
                                                 Lambda = "gamma",
                                                 Decline = "norm",
                                                 ExpectedFecundity="gamma"),
                                     # start=list(Decline=list(prob=.001)),
                                     method="jagsUI")%>%
    dplyr::mutate(m=m,i=i)
  
  if(flushJags){jags.list[[pt]]<-NULL}
  
  if(s%%print.crit==0){
    print(paste(s,"out of",max(ss)))
    
    # This line of code prints the progress as a text file - not necessary
    # cat(paste(s,object.size(jags.list)),file="progress.txt")}
  
  
  
}



# Create a database
fullHyperPars<-rbind.fill(hPars)%>%
  # For each estimated parameter
  dplyr::group_by(EstPar)%>%
  # Create a variable that is explicitly which is the index of said parameter (extracted from between brackets)
  dplyr::mutate(ParIndex=gsub(x=empty2na(str_extract_all(EstPar, "\\[[^()]+\\]")[[1]]),"\\[|\\]",""))%>%
  dplyr::ungroup()



fullTrueValues<-rbind.fill(trueValues)


Lambda<-fullHyperPars%>%
  dplyr::filter(TruePar=="Lambda")%>%
  dplyr::select(ParIndex,Mean,m,i)%>%
  dplyr::rename(mgmt=ParIndex,lambda=Mean)

Decline<-fullHyperPars%>%
  dplyr::filter(TruePar=="Decline")%>%
  dplyr::select(ParIndex,Mean,m,i)%>%
  dplyr::rename(mgmt=ParIndex,probDecline=Mean)

Objectives<-merge(Lambda,Decline,by=c("mgmt","m","i"))%>%
  dplyr::mutate(mgmt=as.numeric(mgmt))%>%
  data.table(key="mgmt")%>%
  merge(data.table(managementOptions%>%
                     dplyr::select(mgmt,CostMgmt),key="mgmt"))%>%
  dplyr::mutate(CostScore=ObjectiveScore(CostMgmt,minimize = T)*Weights[2],
                BioScore=ObjectiveScore(lambda,
                                        Max=globalLimits$lambda$max,
                                        Min=globalLimits$lambda$min)*Weights[1],
                BioScore2=ObjectiveScore(probDecline,Min = 0,Max=0.5,minimize = T)*Weights[1],
                Score=CostScore+BioScore)


system.time({resu3<-Objectives%>%
  # For each monitoring m and paramater space index i
  dplyr::group_by(m,i,mgmt)%>%
  dplyr::summarise(
    lambda=mean(lambda),
    probDecline=mean(probDecline),
    Score=mean(Score),	  
    BioScore=mean(BioScore),
    BioScore2=mean(BioScore2),
    Score2=mean(BioScore2+CostScore),
    CostScore=mean(CostScore))%>%
  # Only keep the management that maximizes score
  # filter(Lambda==max(Lambda))%>%
  dplyr::group_by(m,i)%>%
  dplyr::mutate(MaxScore=Score==max(Score))%>%
  dplyr::mutate(MaxScore2=Score2==max(Score2))%>%
  dplyr::ungroup()%>%
  # Append true values of parameter space index i 
  merge(fullTrueValues%>%
          dplyr::rename(i=Sim,
                        trueLambda=Lambda,
                        trueDecline=Decline),
        by=c("m","i","mgmt"))%>%
  dplyr::group_by(m,i,mgmt)%>%
  dplyr::ungroup()%>%
  merge(managementOptions,by="mgmt")%>%
  merge(monitoringOptions,by="m")%>%
  dplyr::mutate(Cost=CostMonitoring+CostMgmt)%>%
  dplyr::ungroup()%>%
  dplyr::mutate(trueBioScore=ObjectiveScore(trueLambda,
                                            Min = globalLimits$lambda$min,
                                            Max = globalLimits$lambda$max)*Weights[1],
                trueCostScore=ObjectiveScore(Cost,Min = 0,minimize = T)*Weights[2],
                trueScore=trueBioScore+trueCostScore)%>%
  dplyr::group_by(m,i)%>%
  dplyr::mutate(MaxTrueScore=Score==max(trueScore))%>%	 
  ungroup()})

### Consequence Table for Monitoring Options
system.time({finalConsequenceTable<-resu3%>%
  # For each monitoring scenario
  dplyr::filter(MaxScore)%>%
  # dplyr::filter(MaxScore2)%>%
  dplyr::group_by(m)%>%
  # Get the expected outcomes (average of parameter space)
  dplyr::summarise(b=mean(lambda),DeclineProb=mean(probDecline),
                   c=mean(Cost),c2=mean(CostMgmt))%>%
  # Get the expected outcomes (average of parameter space)
  # dplyr::summarise(b=mean(trueLambda),DeclineProb=mean(decline),
  # c=mean(Cost),c2=mean(CostMgmt))%>%
  # Create scores for each outcome
  dplyr::mutate(BioScore=ObjectiveScore(b,
                                        Min = globalLimits$lambda$min,
                                        Max = globalLimits$lambda$max)*Weights[1],
                
                # This might need to change - is the Cost scored based on management or on monitoring?						 
                CostScore=ObjectiveScore(c,Min = 0,Max = ,minimize = T)*Weights[2],
                CostScoreMgmt=ObjectiveScore(c2,Min = 0,Max = ,minimize = T)*Weights[2],
                
                Score=BioScore+CostScore,
                nSims=nrow(generatedData[[1]]$samples[[1]]))%>%
  # Append monitoring information
  merge(monitoringOptions)%>%
  dplyr::arrange(desc(Score))%>%
  dplyr::ungroup()%>%
  ### Calculate expected value for each monitoring scenario
  dplyr::mutate(EVSI_Lambda=b-b[m==Uncertainty_index],
                EVSI_DeclineProb=DeclineProb-DeclineProb[m==Uncertainty_index],
                EVSI_Cost = c-c[m==Uncertainty_index],
                EVSI_CostMgmt = c2-c2[m==Uncertainty_index],
                EVSI_Score = Score - Score[m==Uncertainty_index])})


# Values to save from analysis
VOI<-list()
# The final consequence table sa
VOI$finalConsequenceTable<-finalConsequenceTable
# Results summarised by iteration, monitoring scenario and management alternative
VOI$results.summ<-list(resu3=resu3)
VOI$monitoringOptions<- monitoringOptions
VOI$managementOptions<-managementOptions
# A list of the priors used to generate the results
VOI$priors<-priors
# The global limits of the fundamental objectives (in this case, only Lambda)
VOI$globalLimits<-globalLimits  
# A dataset of the actual values (Lambda, survival and fecundity of each class under each management alternative)
VOI$fullTrueValues<-fullTrueValues
# Full hyper parameters of posterior distributios of each variable of interest. For each iteration under each monitoring alternative
VOI$fullHyperPars<-fullHyperPars
# A list with the full predicted outcomes, and a summary of those
VOI$PredictedOutcomes<-list(Full=PredictedOutcomesFull,
                            Summary=PredictedOutcomesSummary)
# The weights used in the decision making
VOI$Weights<-Weights

save(VOI,file="Chapter 3 - Framework 1 - Results.RData")


