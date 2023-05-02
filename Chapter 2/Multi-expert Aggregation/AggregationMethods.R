### Aggregation Methods

require(devtools)
require(dplyr)

#Import custom functions
devtools::source_url("https://raw.githubusercontent.com/KenupCF/ThesisPHD/main/GeneralFunctions/PriorSampling.R")
devtools::source_url("https://raw.githubusercontent.com/KenupCF/ThesisPHD/main/GeneralFunctions/QuickFunctions.R")



# Load example data
load(file=".\\Example Data v2.RData")



# Calculate average for each question, using quantile aggregation (i.e. averaging across experts)

QuantileAggregation<-Ch2$Answers%>%
  dplyr::group_by(question)%>%
  dplyr::summarise(min=mean(MinQ),mode=mean(mode),max=mean(MaxQ))%>%
  dplyr::mutate(alias="Quantile Aggregation")


# Now calculate the fitted aggregated value for each question

# Create a varaible that is the combination of expert and question
Ch2$Answers<-Ch2$Answers%>%
  mutate(aliasQuestion=paste(alias,question,sep="_"))

# Split data.frame by such variable
splitted<-split(Ch2$Answers,Ch2$Answers$aliasQuestion)
# Convert to simple data.frame class
splitted<-lapply(splitted,function(x){class(x)<-"data.frame";return(x)})
# Check results
str(splitted)


# Sample 1000 values for each expert distribution and aggregate them by question
samples<-priorSampling(splitted,method="lhs",size=1e4)%>%
  reshape2::melt()%>%
  dplyr::mutate(aliasQuestion=variable,variable=NULL)%>%
  merge(Ch2$Answers%>%
          dplyr::select(aliasQuestion,question,alias,
                        originalDistribution),
        all.x=T,all.y = F)%>%
  dplyr::mutate(dist=originalDistribution)


# Nudging values to not be exactly 0 or 1 (bad for beta estimation)
samples$value[samples$value==0]<-0+1e-3
samples$value[samples$value==1]<-1-1e-3

# Now split by question
samples.split<-split(samples,samples$question)




# For each set of sampled values, fit a distribution

fitted<-lapply(
    names(samples.split),
    function(n){
  
  x<-samples.split[[n]]
  
  y<-fitdistrplus::fitdist(
      data=x$value,
      distr = unique(x$dist))
  
  z<-as.data.frame(t(y$estimate))%>%
    dplyr::mutate(dist=y$distname,question=n,alias="Fitted aggregation",
                  aliasQuestion = paste(alias,question,sep="_"))
  
  return(z)
  
})

str(fitted)

# Saving to example file 

Ch2$Aggregation<-list(quantile=QuantileAggregation,fitted=fitted)


save(Ch2,file=".\\Example Data v3.RData")





