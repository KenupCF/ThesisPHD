###Sample values from priors

priorSampling<-function(L,size,seed=NULL,method="random"){
  
  require(lhs)
  require(mc2d)
  
  if(!is.null(seed)){set.seed(seed)}
  
  
  #generate a list of sampling functions for each distribution
  rFUN<-list(beta=rbeta,norm=rnorm,unif=runif,binom=rbinom,
                gamma=rgamma,poisson=rpois,lnorm=rlnorm,exp=rexp,pert=rpert)
				
  dFUN<-list(beta=dbeta,norm=dnorm,unif=dunif,binom=dbinom,
                gamma=dgamma,poisson=dpois,lnorm=dlnorm,exp=dexp,pert=dpert)
  qFUN<-list(beta=qbeta,norm=qnorm,unif=qunif,binom=qbinom,
                gamma=qgamma,poisson=qpois,lnorm=qlnorm,exp=qexp,pert=qpert)
 
  npars<-list(beta=2,norm=2,unif=2,binom=2,gamma=2,poisson=1,lnorm=2,exp=1,pert=4)
 
  if(method%in%c("lhs","oalhs")){
	  require(lhs)}
	  
  if(method==c("lhs"))   {M<-randomLHS(size,length(L))}
  if(method==c("oalhs")) {M<-create_oalhs(size,length(L),TRUE,FALSE)}
  if(method==c("random")){M<-matrix(runif(size*length(L)),nrow=size,ncol=length(L))}
  
  resu<-sapply(1:length(L),function(i){
	  
		  y<-names(L)[i]
		  values<-M[,i]  
		  
		  #get prior data.frame object for each parameter
		  x<-L[[y]]
		  
		  dist<-as.character(x$dist)
		  
		  discrete<-na2false(x$integer)
		  if(length(discrete)==0){discrete<-FALSE}
		  
		  x<-x[,
			   colnames(x)%in%c("min","mode","max","shape","rate","mean","sd","lambda","shape1","shape2")]
		  
		  x$dist<-dist
		
		  #if distribution is not poisson nor exponetial (they need only 1 parameter)
		  if(npars[[x$dist]]==2){
		  z<-qFUN[[x$dist]](p=values,
			x[1,1]-(1*(dist=="unif")*discrete),
			x[1,2]) #generate n values using meta parameters
		  #if distribution IS poisson nor exponetial (they need only 1 parameter)
		  }
		  
		  if(npars[[x$dist]]==1){
		  z<-qFUN[[x$dist]](p=values,x[1,1])        #generate n values using meta parameters
		  }
		  
		  if(npars[[x$dist]]==4){
		  z<-qFUN[[x$dist]](p=values,x[1,1],x[1,2],x[1,3],x[1,4])        #generate n values using meta parameters
		  }
		  
		  if(discrete & dist=="unif"){z<-ceiling(z)}
		  
		  z<-matrix(z,nrow=1) #convert values to a matrix of 1 row
		  
		  rownames(z)<-y #rename the matrix with the parameter name
		  
		  return(z)
		  })
		  
	if(is.null(dim(resu))){resu<-t(resu)}
	resu<-as.data.frame(resu)
	colnames(resu)<-names(L)
	# if(!is.null(seed)){}
	return(resu) 
	
	}
