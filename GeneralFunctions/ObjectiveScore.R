#Function to calculate an objective's score
#(simple standardization to the 0-1 scale)
ObjectiveScore<-function(x,minimize=FALSE,Max=NULL,Min=NULL,clip=TRUE){
  if(is.null(Max)){Max=max(x)}
  if(is.null(Min)){Min=min(x)}
  
  if(clip){
	  x[x<Min]<-Min
          x[x>Max]<-Max}	
  
  y<- (x - Min) / (Max-Min)
  
  if(minimize){y<- 1 - y}
  
  return(y)}
