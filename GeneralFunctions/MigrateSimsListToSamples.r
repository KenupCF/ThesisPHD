###Function to create analogues to variables in $sims.list to $samples
###Useful if creating derived parameters from the MCMCchains in sims.list directly through R
MigrateSimsListToSamples<-function(jags){
	require(dplyr)
	require(stringr)
	require(reshape2)
	
	jags<-list(sims.list=jags$sims.list,samples=jags$samples)
	
	n1<-names(jags$sims.list)
	n2<-unique(gsub(x=dimnames(jags$samples[[1]])[[2]],"\\[.*$",""))
	migrate<-n1[!n1%in%n2]
	
	if(length(migrate)>0){
	
	migrated<-lapply(migrate,function(n){
		x<-jags$sims.list[[n]]
		molten<-reshape2::melt(x)
		VarCols<-str_subset(colnames(molten),"Var")
		VarCols<-VarCols[VarCols!="Var1"]
		if(length(VarCols)>0){
			molten$Index <-paste0("[",apply(as.matrix(molten[,VarCols]),1,paste0,collapse=","),"]")
		}else{molten$Index<-""}
		molten[,VarCols]<-NULL
		if(!"Var1"%in%colnames(molten)){molten$Var1<-1:nrow(molten)}
		molten$Chain<-rep(1:length(jags$samples),each=nrow(jags$samples[[1]]))
		molten$Par<-paste0(n,molten$Index,sep="")
        molten$Index<-NULL
		splitted<-split(molten,molten$Chain)
		resu<-lapply(splitted,function(df){reshape2::acast(df,Var1~Par)})
		return(resu)
	})
	names(migrated)<-migrate
	
	for(ch in 1:length(jags$samples)){
	 for(m in 1:length(migrated)){
		jags$samples[[ch]]<-cbind(jags$samples[[ch]],migrated[[m]][[ch]]) 
	 }
	}
	
   }
   return(jags$samples)
}
