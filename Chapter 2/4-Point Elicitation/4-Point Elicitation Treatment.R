
# Set directory to downloaded repository folder
# setwd()

require(devtools)
require(dplyr)

# Import 4-point convertion function
devtools::source_url("https://raw.githubusercontent.com/KenupCF/ThesisPHD/main/Chapter%202/4-Point%20Elicitation/SolvePERT.R")

load(file=".\\Example Data.RData")

# Calculate a' and c' through linear interpolation
Ch2$Answers<-Ch2$Answers%>%
  mutate(MinLinExt=mode - ((mode-min)*(100/confidence)),
         MaxLinExt=mode + ((max-mode)*(100/confidence)))

# Now calculating it through quantile fitting

# Empty vectors that will hold values
MinQ<-numeric()
MaxQ<-numeric()

for (i in 1:nrow(Ch2$Answers)){
  
  # Run solvePERT function over each answer
  quantileFitting<-solvePERT(Mode=Ch2$Answers$mode[i],
                             Min=Ch2$Answers$min[i],
                             Max=Ch2$Answers$max[i],
                             Conf = Ch2$Answers$confidence[i]/100)
  
  # Save answers to vector
  MinQ[i]<-quantileFitting$trueMin
  MaxQ[i]<-quantileFitting$trueMax
  
  print(paste0(i," out of ",nrow(Ch2$Answers)))  
}

Ch2$Answers$MinQ<-MinQ
Ch2$Answers$MaxQ<-MaxQ

save(Ch2,file=".\\Example Data v2.RData")
