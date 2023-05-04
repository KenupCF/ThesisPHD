##########################################################################
### Transforming elicited distributions into linear model coefficients ###
##########################################################################

require(plyr)
require(dplyr)
require(tidyr)
require(devtools)
require(ggplot2)
require(gtools)
require(jagsUI)
require(pracma)


#Named list of possible mean functions
meanFUN<-list("beta"=function(x){
  a<-x$shape1
  b<-x$shape2
  return(a/(a+b))},
  "gamma"=function(x){x$shape/x$rate},
  "poisson"=function(x){x$lambda},
  "norm"=function(x){x$mean})

#Named list of possible variance functions
varFUN<-list("beta"=function(x){
  a<-x$shape1
  b<-x$shape2
  return((a*b)/(((a+b)^2)*(a+b+1)))},
  "gamma"=function(x){x$shape/(x$rate^2)},
  "poisson"=function(x){x$lambda},
  "norm"=function(x){x$sd^2})

#Named list of possible link functions
linkFUN=list("beta"=gtools::logit,
             "gamma"=log,
             # identity function
             "norm"=function(x){x})

#Named list of possible quantile functions
qFUN=list(beta=qbeta,gamma=qgamma,norm=qnorm)


# List of elicited distributions, where M0 is status quo and M1-3 are an applied management action
# Note that for simplicity, all management actions have the same central estimate, but 
elicited<-list(M0=data.frame(mean=0,sd=1.0,dist="norm"),
               M1=data.frame(mean=1,sd=1.0,dist="norm"),
               M2=data.frame(mean=1,sd=0.5,dist="norm"),
               M3=data.frame(mean=1,sd=2.0,dist="norm"))



###
# Create JAGS script to indepedently sample from the distributions to generate betas
# (Assumption: No covariance) between values of all three distributions
###

{
  sink("model_no_covariance.txt")
  cat(paste("
model{

   SQ~dnorm(",elicited$M0$mean,",",(elicited$M0$sd^-2),")
   Alt1~dnorm(",elicited$M1$mean,",",(elicited$M1$sd^-2),")
   Alt2~dnorm(",elicited$M2$mean,",",(elicited$M2$sd^-2),")
   Alt3~dnorm(",elicited$M3$mean,",",(elicited$M3$sd^-2),")
 
 # Generate Coefficients
 
 alpha<-(SQ)
 beta1<-(Alt1) - alpha
 beta2<-(Alt2) - alpha
 beta3<-(Alt3) - alpha
}"),fill=T)
  sink()
}


# Run the script
no_cov_coef<-jagsUI::jags(
  # Dummy data object - no data needed since this is a prediction based on priors only
  data=list(x=1),
  model.file="model_no_covariance.txt",
  parameters.to.save= c("alpha","beta1","beta2","beta3"),
  n.iter=10e3,
  n.adapt=500,
  n.burnin=1e3,
  n.thin=2,
  n.chains=2)


###
# Generate coefficiennts assuming maximum covariance
###


# Vector of arbitraly large number of quantiles to calculate

q.vec<-seq(from=10e-5,to=1-10e-5,length.out=10e3)


# Quantile extractions
# Note: the quantile function is the inverse function of the CDF
qformed=lapply(elicited,function(x){
  qFUN[[x$dist]](
    q.vec,x[1,1],x[1,2])
})

# Convert each quantile value calculated with with a link function
qformed.link<-lapply(names(elicited),function(n){
  x<-elicited[[n]]
  return(linkFUN[[x$dist]](qformed[[n]]))
})
# rename object
names(qformed.link)<-names(elicited)

# Use these (ordered) values to generate a distribution of coefficients 
perf_cov_coef.sample<-list(
                    alpha=qformed.link$M0,
                    beta1=qformed.link$M1 - qformed.link$M0,
                    beta2=qformed.link$M2 - qformed.link$M0,
                    beta3=qformed.link$M3 - qformed.link$M0)

# Integrate those values over the quantiles to get mean and SD
# (trapz function returns a trapezoid approximation of the integral)
perf_cov_coef<-lapply(perf_cov_coef.sample,function(y){
  # m is mean, s is standard deviation
  m <- pracma::trapz(x=q.vec,y=y)
  s <- sqrt(pracma::trapz(x=q.vec,y=(y-m)^2))
  return(c(mean=m,sd=s))
})

# Same operation, but returning results as a data.frame
perf_cov_coef.df<-lapply(perf_cov_coef.sample,function(y){
  # m is mean, s is standard deviation
  m <- pracma::trapz(x=q.vec,y=y)
  s <- sqrt(pracma::trapz(x=q.vec,y=(y-m)^2))
  return(data.frame(mean=m,sd=s,dist="norm"))
})


###
# Comparing the two methods 
# on how they perform in predicting the original elicited distributions
###

sink("model_terms_to_outcomes.txt")
cat(paste(
  "model{

alpha~dnorm(alpha.input[1],1/pow(alpha.input[2],2))
beta1~dnorm(b1.input[1],1/pow(b1.input[2],2))
beta2~dnorm(b2.input[1],1/pow(b2.input[2],2))
beta3~dnorm(b3.input[1],1/pow(b3.input[2],2))

M0<-alpha
M1<-alpha + beta1 
M2<-alpha + beta2
M3<-alpha + beta3
# Aggregated alternative extrapolating from the generated betas
M1_2<-alpha + beta1 + beta2

}"),fill=T)
sink()



### Generating prediction using the MCMC method (no covariance)

# Input object for the generated coefficients assuming no covariance
data.list.no_cov=list(
  alpha.input=c(
    mean(no_cov_coef$sims.list$alpha),
    sd(no_cov_coef$sims.list$alpha)),
  b1.input=c(
    mean(no_cov_coef$sims.list$beta1),
    sd(no_cov_coef$sims.list$beta1)),
  b2.input=c(
    mean(no_cov_coef$sims.list$beta2),
    sd(no_cov_coef$sims.list$beta2)),
  b3.input=c(
    mean(no_cov_coef$sims.list$beta3),
    sd(no_cov_coef$sims.list$beta3)))



