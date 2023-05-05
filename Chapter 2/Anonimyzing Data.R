### Creating anonymized data
require(dplyr)
require(stringr)

load("C:\\Users\\caiok\\Dropbox\\03-Work\\02-R_Freelancing\\02-Science_Jobs\\DOC_Expert_Elicitation\\Data Handling\\Biological Responses\\Kuaka Elicited v2.RData")
devtools::source_url("https://raw.githubusercontent.com/KenupCF/ThesisPHD/main/GeneralFunctions/QuickFunctions.R")

Ch2_Repo<-list()


###
# Anonimyzed data

Ch2_Repo$Answers<-Ch2$KuakaExampleAnswers%>%
  dplyr::select(email,question,min,mode,max,confidence)%>%
  dplyr::mutate(alias=email,email=NULL)%>%
  dplyr::mutate(originalDistribution="beta",dist="pert",shape=4)

Ch2_Repo$Answers$originalDistribution[str_detect(Ch2_Repo$Answers$question,"MaximumDensity")]<-"gamma"

Ch2_Repo$Answers$alias<-LETTERS[as.numeric(as.factor(Ch2_Repo$Answers$alias))]


Ch2_Repo$Answers$question<-paste0("Q",zero_pad(as.numeric(as.factor(Ch2_Repo$Answers$question)),2),sep="")


Ch2<-Ch2_Repo

save(Ch2,file=".\\Example Data.RData")





