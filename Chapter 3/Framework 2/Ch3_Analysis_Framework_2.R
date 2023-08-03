###---------------------------------------------###
###Check for buggy lines of code, tag = '@buggy'###
###---------------------------------------------###

#Bayesian libraries
require('rjags') # Note: OpenBUGS may work but didn't see it
require(R2jags) # Note: OpenBUGS may work but didn't see it
require(jagsUI)
require(MCMCvis) #visualization of MCMC chains
require(dclone) ##what does it do?
require(mc2d)
require(lhs)
require(parallel)
require(gtools)
# general purpose packages
require(plyr)
require(dplyr)
require(tidyr)
require(tidyverse)
require(magrittr)
require(reshape2)
require(lubridate)
require(gtools)
require(devtools)
require(ggplot2)
require(reshape2)

# Custom functions
# source_url("https://raw.githubusercontent.com/KenupCF/Kenup-Repo/master/phd_experimental_functions.R")

devtools::source_url("https://raw.githubusercontent.com/KenupCF/Kenup-Repo/master/Quick%20Functions.R")
devtools::source_url("https://raw.githubusercontent.com/KenupCF/ThesisPHD/main/GeneralFunctions/PriorSampling.R")

####How can we detect declines in Kaka (due to Stoat predation or otherwise) through mist-nets counts and differential sex ratios?

###HYPOTHETICAL SITUATION
#we get a starting population of 10 female and 10 male Kaka
#sex ratio is 1:1 - p(f) = p(m) = 0.5
#survival is at first equal during the first 15 years, than stoat density increases to the point where female survival = 50% of male survival. This goes on for 35 more years


#Monitoring is kept constant throughout the 50 years. There are three ways to detect something went wrong with the population:
#1) a significant decline was detected (what is `significant?`)
#2) we realize female survival dropped in relation to male survival
#3) the sex ratio of individuals caught on mistnets becomes male-biased, signalling a change in mortality


#All three ways can be detected using the same monitoring scheme (mist-netting kakas, measuring and banding them.)

#Assumptions:
#1) Sex is determined with 0% error rate
#2) Detection is constant over time and sex


#Simulating what is happening with the population

start.<-now() # UNIX format time that analysis started

noSims<-9 # number of simulations to run

totalYears<-30 # number of years to simulate

sampling.method<-"oalhs"

usePriors<-TRUE


###Deterministic values
if(usePriors==FALSE){
trueParams<-expand.grid(
d=.5,      			#probability of catching a given individual in a net 
phi.fem=0.92,    	#yearly survival of females
phi.mal=0.92, 		#yearly survival of males
fec=.7,				#mean number of surviving fledgings per female
K=100,			 	#maximum number (of females)
MaxAge=24,			#maximum age of reproduction
p.female=.5, 		#expected proportion of females on each offspring
yearOfTheStoat=Inf,	#year where the stoats start having an effect of female kakas
ParSim=1:noSims)
}
###Values sample from priors

# loading priors from Expert Elicitation

load("D:\\03-Work\\01-Science\\01-PhD\\01-Expert Elicitation\\Expert Results.RData") #Expert Elicitation Stuff


#################################################
## Using priors to generate possible scenarios ##
################################################


if(usePriors==TRUE){
priors<-list(
  #Probability of nest detection: beta  (conjugate with binomial)
  # d=ExpertResults$Consensus$KakaCaptureProbability_Consensus,
  
  #@note_caio quick fix while R data does not arrive
  d=data.frame(shape1=10,shape2=90,dist="beta"),
  
  
  phi.alpha=ExpertResults$ElicitedTerms$logit.alpha,
  phi.beta.sex=ExpertResults$ElicitedTerms$logit.beta.sex,
  phi.beta.mast=ExpertResults$ElicitedTerms$logit.beta.mast,
  phi.beta.mast.sex=ExpertResults$ElicitedTerms$logit.beta.mast.sex,
  fec=ExpertResults$Consensus$KakaFecundityWithMast_Consensus,
  
  MaxAge=data.frame(min=24,max=24,dist="unif"),
  yearOfTheStoat=data.frame(min=16,max=21,dist="unif"),
  ### Probability of a year being a mast year is not properly defined
  mast.prob=data.frame(shape1=20,shape2=40,dist="beta"),
	# K=data.frame(shape=8,rate=3,dist="gamma"),
	K=data.frame(min=100,max=100,dist="unif"),
	# test=data.frame(mean=8,sd=2,dist="norm"))
  p.female=data.frame(shape1=3e3,shape2=3e3,dist="beta"))



#Creat a data.frame with that holds the 'true value' of each parameter on each simulation
trueParams<-priorSampling(priors,size = noSims,
                          seed=26051991,method="oalhs")%>%
  dplyr::mutate(ParSim=sample(1:n(),replace=F))%>%
  dplyr::arrange(ParSim)
}


#Scenarios are things you want to explicitly compare
scenarios<-list(
	# yearOfTheStoat=c(17:21),
  phi.StoatSexEffect=c(-2,-1),
	dummy=1)

trueParams<-expand.grid.mixed(data=trueParams,
                              new.index="BioSim",
                              old.index="ParSim",L=scenarios)

trueParams$yearOfTheStoat%<>%ceiling

nrow(trueParams)



# trueParams$d<-.4
trueParams$mast.prob<-0

#for the first survey, all individuals are alive, unbanded and unseen

#No of intially marked
no.mark.0=50

EncounterHistoryBegin<-expand.grid(
	Year=0,Sex=c(0,1),
	Ind=1:(no.mark.0/2),
	Alive=1,
	Cohort=1,
	Age=4,
	Seen=1)%>%
	dplyr::mutate(Ind=randomStrings(n()))


simList<-list()

paralleling=TRUE
nchains=1

#Dictionary of variables for survival
var.dict<-data.frame(
  Variable=c("Intercept","Sex",
             "Mast","Sex:Mast"),
  JAGS_coef=c("phi.alpha","phi.beta.sex",
              "phi.beta.mast","phi.beta.mast.sex"))


# Code to set up parallel processing 
if(paralleling){
  
noCores<-min(detectCores()-1, nrow(trueParams))

cl<-makeCluster(noCores) #make cluster of processing

split.sims<-split(1:nrow(trueParams),
                  sort(1: nrow(trueParams)%%noCores)) #split simulation into (roughly) equal chunks
clusterEvalQ(cl, library(plyr))#load R2OpenBUGS on each cluster
clusterEvalQ(cl, library(dplyr))#load R2OpenBUGS on each cluster
clusterEvalQ(cl, library(gtools))#load R2OpenBUGS on each cluster
# clusterEvalQ(cl, load.module("glm"))#load R2OpenBUGS on each cluster
# clusterEvalQ(cl, library(jagsUI))#load R2OpenBUGS on each cluster
clusterExport(cl,c("trueParams",
                   "EncounterHistoryBegin",
                   "totalYears",
                   "var.dict",
                   "jags.predict",
                   "randomStrings",
                   "simList"))                #Export all objects from main enviroment to each cluster enviroment
}

system.time({simResu<-parLapply(cl, split.sims, function(s){


# system.time({
for(i in 
    s){
    # 1:nrow(trueParams)){

###Generating encounter history for the population

EncounterHistoryLoop<-EncounterHistoryBegin

# Female survival 
# phi.femLoop <-trueParams$phi.fem[i]
# phi.alphaLoop<-trueParams$PhiAlpha[i]

# Male survival
# phi.malLoop<-trueParams$phi.mal[i]
# phi.betaSexLoop<-trueParams$PhiBetaSex[i]
}
phi.StoatSexLoopEffect<- trueParams$phi.StoatSexEffect[i]  

# phi.map<-data.frame(Sex=0:1,Phi=c(phi.femLoop,phi.malLoop))

mast.track<-numeric()
for(y in 1:totalYears){

StoatIncrease<-0  

if(y>=trueParams$yearOfTheStoat[i]){StoatIncrease<-1}

# Stochastic sampling whether a given year is a mast year
# (not biologically accurate)
isMast<-rbinom(n = 1,size=1,prob = trueParams$mast.prob[i])
mast.track[y]<-isMast


##########################################  
## Survival sampling for the population ##
##########################################

  EH.prev.year<-EncounterHistoryLoop%>%
	  dplyr::filter(Year==y-1)%>%
    dplyr::ungroup()%>%
    dplyr::mutate(Mast=rep(isMast,n()),Stoat=StoatIncrease)



  
  if(nrow(EH.prev.year)>0){
  
  # Predicting the survival for each individual of the population, 
  # according 
  phi.pred<-unlist(jags.predict(jags=list(sims.list=as.list(trueParams[i,])),
                            var.dict=var.dict,
                            new.data=EH.prev.year)$summary$median)
  
  EH.prev.year$Phi<-phi.pred
    
  EncounterHistoryAlreadyBanded<-EH.prev.year%>%
    
	dplyr::group_by(Ind)%>%
	dplyr::mutate(
	  
	    Year=y,
	    
		  Alive=
		    as.numeric(Age<trueParams$MaxAge[i]) * 
		    Alive * 
		    
		    rbinom(n=1,size=1,prob=
		             
		             #### This line defines the model 
		             #estimating the survival probability for each individual, 
		             #depending on its condition
		             inv.logit(logit(Phi+
		             #@buggy, add effect of stoats on the linear model                  
		             (phi.StoatSexLoopEffect*((Sex==0)*StoatIncrease))))
		           
		           ),
		  
	    Seen=Alive*rbinom(n=1,size=1,prob=trueParams$d[i]),
		  
		  Age=Age+1)

	
  EncounterHistoryLoop<-plyr::rbind.fill(EncounterHistoryLoop,
                                   EncounterHistoryAlreadyBanded)%>%
    #Remove dead individuals from data.frame 
    #(dead and not seen)
    dplyr::filter(!(Alive==0 & Seen==0))%>%
    dplyr::ungroup()
  
############################################  
##Reproduction sampling for the population##
############################################
  
aliveall<-EncounterHistoryLoop%>%
	dplyr::filter(Year==y)%>%
	dplyr::pull(Alive)%>%
    sum
  
alivefemales<-EncounterHistoryLoop%>%
	dplyr::filter(Year==y,Sex==0)%>%
	dplyr::pull(Alive)%>%
  sum

# How many new individuals were generated

newindivsMax<-sum(rpois(n=alivefemales,lambda=trueParams$fec[i]*isMast))

#Capping new individuals 
#This is assumes a ceiling model for the population
newindivs<-min(trueParams$K[i]-aliveall,newindivsMax)

if(newindivs>0){
NewIndiv<-data.frame(Ind=randomStrings(newindivs))%>%
  dplyr::mutate(Mast=isMast,Stoat=StoatIncrease)%>%
	dplyr::group_by(Ind)%>%
	dplyr::mutate(
	    Sex=rbinom(n=1,size=1,prob=1-trueParams$p.female[i]),
		  Alive=1,Age=1,Year=y,Cohort=y,
		  Seen=rbinom(n=1,size=1,prob=trueParams$d[i]))

		  
EncounterHistoryLoop<-plyr::rbind.fill(EncounterHistoryLoop,NewIndiv)%>%
  mutate(Mast=rep(isMast,n()),Stoat=StoatIncrease)
}

}

}	
simList[[i]]<-EncounterHistoryLoop
print(i)

}
i
return(simList)  
})})  
  
if(paralleling){stopCluster(cl)}  

simList<-reconstruct.parallel(simResu)


data.jags.list<-lapply(simList,function(X){
  
  
  
  })


sapply(simList,nrow)


######Simulation runs####
monitorTrueParams<-expand.grid(list(
	BioSim=trueParams$BioSim,
	SurveyFrequency=c(1,2,5)))%>%    ####@adjusting
	dplyr::arrange(BioSim)%>%
	dplyr::mutate(MonSim=1:n())

#@adjusting
timeSeriesBinYear<-1+(max(monitorTrueParams$SurveyFrequency)*2) #setting the time series bin year so as to always use at least 3 surveys

####Generating monitoring data from  population simulation
useMNKA=0
timeSeriesBin=5 #arbitrary

data.jags.list<-list()
truth.list<-list()
start.<-now()
options(warn=-1)

for(i in 1:nrow(monitorTrueParams)){

# The frequency on which Kaka will be monitored
surveyFrequency=monitorTrueParams$SurveyFrequency[i]

# The years that are actually sampled
yearsSampled<-seq(from=1,to=totalYears,by=surveyFrequency)

# Data.frame mapping each survey for a given year
index.dict<-data.frame(Survey=1:length(yearsSampled),Year=yearsSampled)

# Just the encounter history of the individuals that were actually
# caught once or more
# ...................
# We do not `see` any individuals that avoided capture throughout 
# the whole study
EncounterHistorySeenOnce<-simList[[monitorTrueParams$BioSim[i]]]%>%
  dplyr::group_by(Ind)%>%
  dplyr::arrange(Cohort)%>%
	merge(index.dict,by="Year",all.y=T)%>%
  tidyr::complete(
    nesting(Ind,Sex),nesting(Year,Survey),
    fill=list(Seen=0,Alive=0))%>%
  dplyr::filter(!is.na(Ind))%>%
  dplyr::ungroup()%>%
  # Check only the encounter history on the years that were sampled
  dplyr::filter(Year>=1,Year%in%yearsSampled)%>%
  dplyr::group_by(Ind)%>%
  # Keep the individuals that were seen at least once on those years
  # dplyr::filter(sum(Seen)>0)%>%
  dplyr::ungroup()
  
#Cumulative time elapsed between each survey
int<-EncounterHistorySeenOnce%>%
  dplyr::arrange(Year)%>%
  dplyr::pull(Year)%>%
  unique
int<-int[2:length(int)]  -   int[1:(length(int)-1)]
cumTime<-cumsum(c(0,int))

# Individuals are not known to be alive on surveys after their last sighting  	
EncounterHistorySeenOnce<-EncounterHistorySeenOnce%>%
  dplyr::group_by(Ind)%>%
	dplyr::mutate(LastSeen=max(Survey[Seen==1]),
	              FirstSeen=min(Survey[Seen==1]))%>%
	# dplyr::ungroup()%>%
	# tidyr::complete(Ind,Survey,fill=list(Alive=NA,Seen=0))%>%
	# dplyr::mutate(AliveQ=replace(Alive ,Survey<FirstSeen ,0))%>%
	dplyr::mutate(AliveQ=replace(Alive,Survey>LastSeen ,NA))

	
#Individuals ordered by cohort
cohort.order<-EncounterHistorySeenOnce%>%
	ungroup%>%
	filter(!duplicated(Ind))%>%
	arrange(Cohort)%>%
	pull(Ind)
	
# Create matrix of sightings
seen<-reshape2::acast(EncounterHistorySeenOnce,Ind~Survey,value.var="Seen",fill=0)
seen<-seen[sort(rownames(seen)),]
if(is.null(dim(seen))){seen<-t(t(seen))}

# Create matrix of knowledge of "aliveness" (good wording! /s)
alive<-reshape2::acast(EncounterHistorySeenOnce,Ind~Survey,value.var="AliveQ")
alive<-alive[rownames(seen),]
if(is.null(dim(alive))){alive<-t(t(alive))}


#fill alive matrixes with zeroes before first detection
# alive<-t(apply(alive,1,function(x){
	# first.detection<-min(which(x==1))
	# x[0:(first.detection-1)]<-0
	# return(x)
	# }))

# Vector declaring where each individual was detected for the first time
first.detection<-EncounterHistorySeenOnce%>%
  dplyr::group_by(Ind)%>%
  dplyr::summarise(f=min(Survey[Seen==1]))%>%
  dplyr::arrange(Ind)%>%
	pull(f)

# Which individuals are female
females.ind<-EncounterHistorySeenOnce%>%
  dplyr::filter(Sex==0)%>%
  dplyr::pull(Ind)%>%
	unique

# Which individuals are male	
males.ind<-EncounterHistorySeenOnce%>%
  dplyr::filter(Sex==1)%>%
  dplyr::pull(Ind)%>%
	unique
		
# The sex of each individual
sex<-EncounterHistorySeenOnce%>%
  dplyr::group_by(Ind)%>%
  dplyr::summarise(s=unique(Sex))%>%
  dplyr::arrange(Ind)%>%
	pull(s)
	
# Counts of individuals on each survey
Counts<-EncounterHistorySeenOnce%>%
		dplyr::group_by(Survey,Sex)%>%
    dplyr::summarise(Count=sum(Seen))%>%
    dplyr::ungroup()%>%
    tidyr::complete(Survey, Sex, fill = list(Count = 0))
		
# Unbanded individuals on each survey
UB<-EncounterHistorySeenOnce%>%
		dplyr::group_by(Survey,Sex)%>%
		dplyr::summarise(UB=sum(Survey==FirstSeen))%>%
		dplyr::ungroup()%>%
		dplyr::mutate(Survey=as.numeric(Survey),Sex=factor(Sex))%>%
		tidyr::complete(Survey,Sex,fill=list(UB=0))
# Unbanded females on each survey
fem.UB<-UB%>%
  filter(Sex==0)%>%
  arrange(Survey)%>%
  pull(UB)	
# Unbanded males on each survey
mal.UB<-UB%>%
  filter(Sex==1)%>%
  arrange(Survey)%>%
  pull(UB)	
# Counts of females on each survey
fem.counts<-Counts%>%
  filter(Sex==0)%>%
  arrange(Survey)%>%
  pull(Count)	
# Counts of males on each survey
mal.counts<-Counts%>%
  filter(Sex==1)%>%
  arrange(Survey)%>%
  pull(Count)	

# Minimum number known alive
mnka.f<-apply(matrix(alive[females.ind,],ncol=ncol(alive)),2,sum,na.rm=T)
mnka.m<-apply(matrix(alive[males.ind,],ncol=ncol(alive)),  2,sum,na.rm=T)

# Total number of marks added to the population on each survey
total.marks<-fem.UB+mal.UB
total.marks[1]<-total.marks[1]+no.mark.0

# Transform the total into a cumulative sum
total.marks<-cumsum(total.marks)

##########
###BUGS###
##########
# {		
data.jags <- list(
			#Observed values
      seen=seen,
			alive=alive,
			mast = EncounterHistorySeenOnce%>%group_by(Survey)%>%
			  summarise(mast=mean(Mast,na.rm=T))%>%pull(mast),
			stoat =  EncounterHistorySeenOnce%>%group_by(Survey)%>%
			  summarise(stoat=mean(Stoat,na.rm=T))%>%pull(stoat),
			fem.counts=fem.counts,mal.counts=mal.counts,
			fem.UB=fem.UB,mal.UB=mal.UB,total.UB=fem.UB+mal.UB,
			total.counts=fem.counts+mal.counts,
			mnka.f=mnka.f,mnka.m=mnka.m,mnka=mnka.f+mnka.m,
			intervals=int,
			total.marks=total.marks,
			# TRUTH=list(N.true,Trend.true),
			RE.year.detectability=as.numeric(F), #0 or 1
			RE.year.survival=as.numeric(F), #0 or 1
			timeSeriesBinSurvey=timeSeriesBin,
			timeSeriesBinYear=timeSeriesBinYear,
			useMNKA=as.numeric(useMNKA),
            sex=sex,
            first=first.detection,
			index.dict=index.dict,
			years=index.dict$Year,surveys=index.dict$Survey,
			n.indiv=nrow(seen),n.surveys=ncol(seen))
data.jags$observed.marks<-apply(data.jags$seen,2,sum)
data.jags$recaptures<-data.jags$total.counts-data.jags$total.UB

# Measured over two+ surveys
data.jags$n.surveys


TrendIndexEnd<-data.jags$surveys[-1]-1


##THIS IS WRONG, REPAIR @buggy
TrendIndexBegin<-sapply(TrendIndexEnd, #for each interval between surveys
	function(x){

	Y<-data.jags$years[x+1] #year in which interval was calculated
	timeE<-Y-data.jags$years[-length(data.jags$years)] #time elapsed between 

	min(which(timeE%in%(0:(data.jags$timeSeriesBinYear-1))))
	
	})

data.jags$TrendIndexes<-cbind(TrendIndexBegin,TrendIndexEnd)

data.jags.list[[i]]<-data.jags

###Calculating ABSOLUTE truths from each simulation

# Number of individuals on each surveyed year  
N.true<-simList[[monitorTrueParams$BioSim[i]]]%>%
	group_by(Year)%>%
	summarise(N=sum(Alive))%>%
	merge(index.dict,all.x=T,all.y=T)%>%
  tidyr::complete(nesting(Year,Survey),fill=list(N=0))%>%
  mutate(BioSim=monitorTrueParams$BioSim[i],MonSim=i)

TrendSurvey.true<-rbind.fill(
  # For each survey
    lapply(na.rm(unique(N.true$Survey)),function(x){
    
    #The years over which we are calculating a trend  
  	minYear<-N.true%>%
  	  dplyr::filter(Survey==x)%>%
  	  dplyr::pull(Year)
  	
  	maxYear<-N.true%>%
  	  dplyr::filter(Survey==(x+timeSeriesBin-1))%>%
  	  dplyr::pull(Year)
  	
  	
  	if(length(maxYear)>0){
  	df<-N.true%>%
  	  # Get the values between each year
  	  dplyr::filter(Year%in%(minYear:maxYear))%>%
  	  # Arrange them by year
  	  dplyr::arrange(Year)%>%
  	  # Calculate a lambda for each timestep
  		dplyr::mutate(Lambda=c(NA,N[2:n()]/N[1:(n()-1)]))%>%
  	  # Calculate an overall trend over these years (geometric average)
  		dplyr::summarise(
  		  Trend=prod(Lambda,na.rm=T)^(1/diff(range(Year))),
  			Survey=max(Survey,na.rm=T),
  			Year=max(Year),minYear=minYear)}else{
  		df<-data.frame(Trend=0,Survey=0,Year=0)[0,]}
  	return(df)
  	}))%>%
	mutate(BioSim=rep(monitorTrueParams$BioSim[i],n()),MonSim=rep(i,n()))



Lambda.true<-N.true%>%
	dplyr::filter(!is.na(Survey))%>%
	dplyr::arrange(Year)%>%
  dplyr::ungroup()%>%
	dplyr::mutate(Lambda=c(NA,N[2:n()]/N[1:(n()-1)]),
		   TimeElapsed=c(NA,Year[2:n()]-Year[1:(n()-1)]),
		   minYear=c(NA,Year[1:(n()-1)]),
		   N=NULL)%>%
	dplyr::filter(!is.na(Lambda))	  

TrendYear.true<-lapply(			#for each:
	1:nrow(Lambda.true),		#occasion that we have a finite rate of increase 
								#(i.e., each interval between surveys)
	function(x){

	minLambda<-TrendIndexBegin[x]
	maxLambda<-TrendIndexEnd[x]

	ldf<-Lambda.true[minLambda:maxLambda,]%>%
		summarise(Trend=prod(Lambda)^(1/sum(TimeElapsed)),minYear=min(minYear),Year=max(Year))
	return(ldf)
	})%>%
	rbind.fill()%>%
	dplyr::mutate(BioSim=monitorTrueParams$BioSim[i],MonSim=i)

fullTruth<-simList[[monitorTrueParams$BioSim[i]]]
	
truth.list[[i]]<-list(N=N.true,Lambda=Lambda.true,
                      TrendY=TrendYear.true,TrendS=TrendSurvey.true,
                      index.dict=index.dict,fullTruth=fullTruth)





#Adding priors information to data.jags object
data.jags.list[[i]]$priorProbDetection<-as.numeric(priors$d[1,1:2])

# Dealing with never seen individuals

seenOnce<-which(rowSums(data.jags.list[[i]]$seen)!=0)

data.jags.list[[i]]$seen<-data.jags.list[[i]]$seen[seenOnce,]
data.jags.list[[i]]$alive<-data.jags.list[[i]]$alive[seenOnce,]
data.jags.list[[i]]$first<-data.jags.list[[i]]$first[seenOnce]
data.jags.list[[i]]$sex<-data.jags.list[[i]]$sex[seenOnce]
data.jags.list[[i]]$n.indiv<-length(seenOnce)
# Adding a dummy survey

data.jags.list[[i]]$seen<-cbind(data.jags.list[[i]]$seen,NA)
data.jags.list[[i]]$alive<-cbind(data.jags.list[[i]]$alive,NA)
data.jags.list[[i]]$mnka%<>%c(0)
data.jags.list[[i]]$mast%<>%c(0)
data.jags.list[[i]]$stoat%<>%c(1)
data.jags.list[[i]]$n.surveys<-data.jags.list[[i]]$n.surveys+1
data.jags.list[[i]]$fem.counts%<>%c(NA)
data.jags.list[[i]]$mal.counts%<>%c(NA)
data.jags.list[[i]]$total.counts%<>%c(NA)
data.jags.list[[i]]$intervals%<>%c(1)



## Add elicited terms to Phi to bugs.list

for(e in 1:length(ExpertResults$ElicitedTerms)){
  data.jags.list[[i]][[length(data.jags.list[[i]])+1]]<-as.numeric(ExpertResults$ElicitedTerms[[e]][,1:2])
  names(data.jags.list[[i]])[length(data.jags.list[[i]])]<-names(ExpertResults$ElicitedTerms)[e]
  
}

if(i%%100==0){print(i)}

}

options(warn=0)
timeTaken<-difftime(now(),start.)
timeTaken/nrow(monitorTrueParams)
timeTaken

save.image("Temp.RData")
getwd()


jags.file<-".\\Chapter 3\Framework 2\\jags_model.txt"

# Parameters to save
params<-c("ProbDetection","AlfaPhi",
          "BetaSexPhi","BetaMastPhi","BetaSexMastPhi","BetaStoatPhi",
          "Beta.sex.phi.switch",
		  "Lambda","N","re.t.p",
          "TrendSurvey","isDeclineTrendSurvey","Z","first.Trend",
          "TrendYear","isDeclineTrendYear")
          # "Lambda","N","int",
			# "TrueSexRatio","RE.y.p.switch","AlphaFecundity","total.UBs",
			# "live.marks",
			# "B","fec")
	
# system.time({JM <- rjags::jags.model(data=data.jags,
           # file =jags.file , n.chains=1);
		   # JM_out <- rjags::coda.samples(JM,
          # variable.names = params,n.iter = 3e3,thin = 1)		   
		   # })

paralleling=TRUE
runToyMCMC=FALSE
nchains=2


if(paralleling){
noCores<-min(detectCores(),length(data.jags.list))
cl<-makeCluster(noCores) #make cluster of processing
split.sims<-split(1:length(data.jags.list),sort(1:length(data.jags.list)%%noCores)) #split simulation into (roughly) equal chunks
clusterEvalQ(cl, library(rjags))#load R2OpenBUGS on each cluster
clusterEvalQ(cl, load.module("glm"))#load R2OpenBUGS on each cluster
clusterEvalQ(cl, library(jagsUI))#load R2OpenBUGS on each cluster
clusterExport(cl,c("split.sims","params","nchains","data.jags.list","truth.list","jags.file","paralleling","runToyMCMC","extractBugsParsSummary","qsummary"))                #Export all objects from main enviroment to each cluster enviroment
}


if(!paralleling){s<-1:length(data.jags.list)}

####THIS IS THE COMMAND THAT RUNS ALL THE JAGS ITERATIONS

if(paralleling){
system.time({simResu<-parLapply(cl, split.sims, function(s){


JAGS.list<-list()
# JAGS.list2<-list()

for(i in s){
#Command running the JAGS script of interpreting monitoring data

# initial.list<-list()

# for(j in 1:nchains){initial.list[[j]]<-list(N.c=data.jags.list[[i]]$mnka)}

jagsLoop<- jagsUI::jags(
		   data=data.jags.list[[i]],
		   # inits=initial.list,
           model.file=jags.file,
		   parameters.to.save=params,
		   n.iter=ifelse(runToyMCMC,100,10e3),
		   n.adapt=ifelse(runToyMCMC,5,200),
		   n.burnin=ifelse(runToyMCMC,10,2e3),
		   n.thin=ifelse(runToyMCMC,1,2),
		   n.chains=ifelse(runToyMCMC,1,nchains))

JAGS.list[[i]]<-list(summary=extractBugsParsSummary(jagsLoop))
  
#JAGS.list[[i]]$full<-jagsLoop
#JAGS.list[[i]]$data<-data.jags.list[[i]]	   
#JAGS.list[[i]]$params<-monitorTrueParams[i,]	   
#JAGS.list[[i]]$truth<-truth.list[[i]]	
		  
gc()		   

}



# print(which(s == i))


return(JAGS.list)
})})
} else{


start.<-now()
JAGS.list<-list()
# JAGS.list2<-list()
initial.list<-list()

for(j in 1:nchains){initial.list[[j]]<-list(N.c=data.jags.list[[i]]$mnka)}

for(i in s){
#Command running the JAGS script of interpreting monitoring data

jagsLoop<- jagsUI::jags(
		   data=data.jags.list[[i]],
		   # inits=initial.list,
           model.file=jags.file,
		   parameters.to.save=params,
		   n.iter=ifelse(runToyMCMC,100,2e3),
		   n.adapt=ifelse(runToyMCMC,5,200),
		   n.burnin=ifelse(runToyMCMC,10,5e2),
		   n.thin=ifelse(runToyMCMC,1,2),
		   n.chains=ifelse(runToyMCMC,1,nchains))

JAGS.list[[i]]<-list(summary=extractBugsParsSummary(jagsLoop))
  
#JAGS.list[[i]]$full<-jagsLoop
#JAGS.list[[i]]$data<-data.jags.list[[i]]	   
#JAGS.list[[i]]$params<-monitorTrueParams[i,]	   
#JAGS.list[[i]]$truth<-truth.list[[i]]	
		  
print(i)
}


end.<-now()

}

if(paralleling) {stopCluster(cl)}
if(!paralleling){difftime(start.,end.)}

beep(sound=wc2sound)
Sys.sleep(2)


#Reconstructing parallelized object
if(paralleling){
JAGS.list<-list()
for(i in 1:length(simResu)){

x<-simResu[[i]]

valid_models<-!sapply(x,is.null)

indexes <- 1:sum(valid_models)  

indexes <- indexes+length(JAGS.list)

JAGS.list[indexes]<-simResu[[i]][valid_models]


}
# rm(simResu)
}


time.elapsed<-round(abs(difftime(start.,now(),units="hours"),2))

RPushbullet::pbPost(type = "note",title = paste("JAGS Analysis Finished!",time.elapsed,sep=";"),apikey = "o.8v2VQ88kZXsMdUrk86PtaxXfbsaPOwMm")

#Parsing results from all simulations at once
popResults<-lapply(1:length(JAGS.list),function(i){

#Getting information from each simulation
x<-JAGS.list[[i]] #object containing JAGS results
data<-data.jags.list[[i]]	  #object containing JAGS input
BioSim<-monitorTrueParams$BioSim[i]
MonSim<-monitorTrueParams$MonSim[i]
resu<-JAGS.list[[i]]$summary

#Create mapping for each Trend Index and its corresponding survey
#(in this case, the index of the last survey used to calculate trend)

trendYear.df<-data$TrendIndexes%>%
	as.data.frame%>%
	dplyr::mutate(Index1=1:n(),Survey=TrendIndexEnd+1)%>%
	dplyr::select(Index1,Survey)

###Time-dependent variables
trendYear.vars<-c("TrendYear","isDeclineTrendYear")
trendSurvey.vars<-c("TrendSurvey","isDeclineTrendSurvey")
survey.vars<-c("N")
survey.int.vars<-c("Lambda")

temp<-list()

#Arranging variables that are related to each survey
temp[[length(temp)+1]]<-resu%>%
	dplyr::filter(Par%in%survey.vars)%>%
	dplyr::mutate(Survey=Index1)

temp[[length(temp)+1]]<-resu%>%
	dplyr::filter(Par%in%survey.int.vars)%>%
	dplyr::mutate(Survey=Index1+1)

temp[[length(temp)+1]]<-resu%>%filter(Par%in%trendYear.vars)%>%
	merge(trendYear.df,by="Index1",all=T)
	
temp[[length(temp)+1]]<-resu%>%
	filter(Par%in%trendSurvey.vars)%>%
	dplyr::mutate(Survey=Index1+data$timeSeriesBinSurvey-1)	
	
	
resu2<-plyr::rbind.fill(temp)%>%
	merge(data$index.dict%>%select(Survey,Year),by="Survey",all.x=TRUE)%>%
	dplyr::arrange(Par,Survey)%>%
	dplyr::mutate(BioSim=BioSim,MonSim=MonSim)
	
return(resu2)
})%>%plyr::rbind.fill()

popResults<-merge(popResults,monitorTrueParams,by=c("BioSim","MonSim"))

#Split by parameter of interest
popResults<-split(popResults,popResults$Par)

save(popResults,truth.list,trueParams,monitorTrueParams,priors,file="wd.RData")

save(JAGS.list,file="JAGS.list.RData")


#Dealing with trends
trendCompare<-merge(popResults$isDeclineTrendSurvey,trueParams%>%
	dplyr::select(BioSim,yearOfTheStoat))%>%
	merge(monitorTrueParams)
	
trendTruth<-rbind.fill(lapply(truth.list,function(x){x$TrendS}))

PopSizeTruth<-rbind.fill(lapply(truth.list,function(x){x$N}))%>%
	dplyr::group_by(BioSim,MonSim)%>%
	dplyr::arrange(Year)%>%
	dplyr::mutate(LostIndividuals=cumsum(clip(min=0,c(0,N[1:(n()-1)]-N[2:n()]))),
				 LostIndividualsPercent=LostIndividuals/max(N))%>%
	dplyr::ungroup()%>%
	dplyr::arrange(BioSim,MonSim,Year)


trendCompare<-merge(trendCompare,
					trendTruth%>%dplyr::select(Trend,MonSim,Year),
					by=c("MonSim","Year"),all.x=T,all.y=F)%>%
					merge(PopSizeTruth%>%
						dplyr::select(N,LostIndividuals,LostIndividualsPercent,MonSim,Year),
					by=c("MonSim","Year"),all.x=T,all.y=F)%>%
					dplyr::arrange(BioSim,MonSim,Year)
					
trendCompare2<-trendCompare%>%
	dplyr::group_by(yearOfTheStoat,SurveyFrequency,BioSim,MonSim)%>%
	dplyr::arrange(Year)%>%
	dplyr::summarise(YearOfTrueDetection=
		min(first(Year[mean>0.95 & Trend<1]),last(Year),na.rm=T),
		YearsUntilDetection=YearOfTrueDetection-first(Year[Trend<1]),
		LostIndividualsBeforeDetection=
		min(first(LostIndividuals[mean>0.95 & Trend<1]),last(LostIndividuals),na.rm=T))%>%
	dplyr::arrange(yearOfTheStoat,SurveyFrequency)
	
trendCompare3<-trendCompare2%>%	
	dplyr::group_by(SurveyFrequency)%>%
	dplyr::summarise(
		YearOfTrueDetection=mean(YearOfTrueDetection),
		YearsUntilDetection=mean(YearsUntilDetection),
		ChanceLosing50Indiv=sum(LostIndividualsBeforeDetection>50)/n(),
		LostIndividualsBeforeDetection=mean(LostIndividualsBeforeDetection))
						

	
S=1


ggplot(data=popResults$N%>%filter(BioSim%in%S),aes(y=median,x=Year))+
	geom_line()+
	geom_point()+
	geom_errorbar(aes(ymin=lcl,ymax=ucl))+
	geom_vline(data=trueParams%>%filter(BioSim%in%S),
		mapping=aes(xintercept=yearOfTheStoat),linetype="dashed",color="red")+ #for vertical lines
	geom_vline(data=trendCompare2%>%filter(BioSim%in%S),
		mapping=aes(xintercept=YearOfTrueDetection),linetype="dashed")+
	geom_point(data=PopSizeTruth%>%filter(BioSim%in%S),aes(x=Year,y=N),color="red")+
	xlim(0,max(yearsSampled)+.5)+
	ylab("Population Size")+
	facet_grid(~ SurveyFrequency)+
	# ylim(0,200)+
	theme_bw()	
	
S=S+1	

	
	





###Draft - This bit does not run
if(1==2){

####Graphical ways of showing influence on survey frequency
lostIndivDens<-lapply(unique(trendCompare2$SurveyFrequency),function(x){
lostIndivs<-trendCompare2%>%filter(SurveyFrequency==x)%>%pull(LostIndividualsBeforeDetection)
lostIndivsDens<-density(lostIndivs)
return(lostIndivsDens)
})


#Heat map?
temp<-list()
for(i in 1:max(trendCompare2$LostIndividualsBeforeDetection)){

temp[[i]]<-trendCompare2%>%	
	dplyr::group_by(SurveyFrequency)%>%
	dplyr::summarise(
		ChanceLosingIndiv=sum(LostIndividualsBeforeDetection>i)/n(),
		LostIndiv=i)


}
temp<-rbind.fill(temp)

ggplot(temp,aes(x=as.factor(SurveyFrequency),y=LostIndiv,fill=ChanceLosingIndiv))+
geom_tile(colour=NA)+
scale_y_reverse()+
scale_fill_gradientn(
    			name='Probability',
    			values=c(0,.2,.5,1),
				limits=c(0,1),
    			# limits=c(0,max(daily.coffee$Intake,na.rm=T)),
    			# colours=c('#FFFFFF',"#898989","#000000"),na.value = "grey75",
    			colours=c('#72B3B3',"#F3C98B","#F41000","#F41000"),na.value = "grey75",
    			breaks=c(0,.2,.5,.8,1))+
ylab("Number of Individuals Lost")+xlab("Survey Frequency")+
theme(
				panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    			panel.background = element_blank())
				#, 
				a#xis.line = element_blank())
				
				# ,
    			)


}




