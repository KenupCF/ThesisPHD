
zero_pad<-function(x,w){
require(stringr)
y<-str_pad(x, w, side = "left", pad = "0")
return(y)}	


na2false<-function(x){
  x[is.na(x)]<-FALSE
  return(x)
}


null2false<-function(x){if(length(x)==0){x<-FALSE};return(x)}	
	
	step<-function(x){
		y<-rep(0,length(x))
		y[x>=0]<-1
		return(y)}
		
	
	pow<-function(x,y){
		x^y		
		}
	
    inv.logit<-function (x, min = 0, max = 1){
          p <- exp(x)/(1 + exp(x))
          p <- ifelse(is.na(p) & !is.na(x), 1, p)
          p <- p * (max - min) + min
          return(p)
      }
  	na.rm<-function(x){x<-x[!is.na(x)];return(x)}
  	empty2na<-function(x){if(length(x)==0){x<-NA};return(x)}
  	inf2nan<-function(x){x[x%in%c(Inf,-Inf)]<-NaN;return(x)}
  	nan2zero<-function(x){x[is.nan(x)]<-0;return(x)}
  	na2zero<-function(x){x[is.na(x)]<-0;return(x)}
  	na2false<-function(x){x[is.na(x)]<-FALSE;return(x)}
  	cap<-function(x,max=Inf,min=-Inf){x[x>max]<-max;x[x<min]<-min;return(x)}
  	clip<-function(x,min=-Inf,max=Inf){
  		x[x<min]<-min
  		x[x>max]<-max
  		return(x)}
    randomStrings<-function(n = 5000,no_letters=6,no_numbers=4) {
         a <- do.call(paste0, replicate(no_letters, sample(LETTERS, n, TRUE), FALSE))
          paste0(a, 
                 sprintf(paste("%0",no_numbers,"d",sep=""), 
                         sample(as.numeric(paste0(rep(9,no_numbers),collapse="")), n, TRUE)),
                 sample(LETTERS, n, TRUE))
    }
