#Spectrum parameters for disease progression, etc.
spectrumPars <- read.csv(paste(dirData,country,'ART.txt',sep=''),header=FALSE,row.names=1)
cd4LowerLimits <- as.numeric(spectrumPars['CD4LOWLIMITS',])
numCD4compartments <- length(cd4LowerLimits)
lambdaCD4 <- as.numeric(spectrumPars['LAMBDA',1:(length(cd4LowerLimits)-1)])   #vector of parameter controlling movement between HIV+ not on ART compartments [CD4 to CD4] 
lambdaCD4length <- length(lambdaCD4)
muCD4 <- as.numeric(spectrumPars['MU',])             #AIDS-related mortality in HIV+ not on ART compartments [CD4] (stratified by sex in Spectrum but one matrix for EPP(?))
alpha1CD4 <- as.numeric(spectrumPars['ALPHA1',])        #AIDS-related mortality in HIV+ on ART for first 6 months [CD4]
alpha2CD4 <- as.numeric(spectrumPars['ALPHA2',])      #AIDS-related mortality in HIV+ on ART for second 6 months [CD4]
alpha3CD4 <- as.numeric(spectrumPars['ALPHA3',])    #AIDS-related mortality in HIV+ on ART for more than 1 yr [CD4]
tempBeta1 <- tempBeta2 <- 0  #set to zero for now    
infectLevelART <- 1 - spectrumPars['INFECTREDUC',1]
propNewInfectsCD4_350 <- spectrumPars['NEWINFECTSCD4_350',1]


#Population sizes on ART and percent/numbers for ART coverage
if (is.element(country, c("Botswana", "Kenya", "Uganda", "Lesotho", "Liberia", "Namibia", "Nigeria", "Tanzania", "Ethiopia")))
spectrumARTpops <- read.csv(paste(dirData,country,'ARTpops.csv',sep=''),header=FALSE)
if (!is.element(country, c("Botswana", "Kenya", "Uganda", "Lesotho", "Liberia", "Namibia", "Nigeria", "Tanzania", "Ethiopia")))
spectrumARTpops <- read.csv(paste(dirData,country,'ARTPops.csv',sep=''),header=FALSE)

colnames(spectrumARTpops) <- c('Year','Percent1stLine','Num1stLine','Percent2ndLine','Num2ndLine','Eligibility')
spectrumARTpops <- spectrumARTpops[spectrumARTpops$Year<=endYear,]
temp1stLine <- spectrumARTpops$Num1stLine
temp2ndLine <- spectrumARTpops$Num2ndLine
temp1stLineIsPercent <- temp2ndLineIsPercent <- rep(NA,length(temp1stLine))
temp1stLineIsPercent[spectrumARTpops$Percent1stLine=="P"] <- 1
temp1stLineIsPercent[spectrumARTpops$Percent1stLine=="N"] <- 0
temp2ndLineIsPercent[spectrumARTpops$Percent2ndLine =="P"] <- 1
temp2ndLineIsPercent[spectrumARTpops$Percent2ndLine =="N"] <- 0
tempEligibilityThreshold <- spectrumARTpops$Eligibility


CD4_GT_500in <- 1
CD4_350_499in <- 2
CD4_LT_50in <- 7
EPPConstantsDEFAULT_CD4_THRESHOLD <- 200


#Calculate scaling parameter for apportioning ART between Urban and Rural populations
# EPP bases this on relative population sizes in 2010 [DRH modification: or earlier if projection is halted before 2010]
fnCalcScaling <- function(){
 return(scaling=1)
}



fnGenARTarrays <- function(){

 #Indices and objects
 artEarliestYr <- min(spectrumARTpops$Year)
 postARTStartYrs <- endYear - artEarliestYr
 preARTStartYrs <- artEarliestYr - beginYear - 1
 numART1in <- numART2in <- rep(NA,numYears+1)
 numART1inLength <- length(numART1in)
 art1stIsPercentIn <- art2ndIsPercentIn <- rep(FALSE,numYears+1)
 eligibilityThresholdCD4 <- rep(NA,numYears+1)

 numART1 <- rep(NA,numsteps)
 numART2 <- rep(NA,numsteps)
 art1stIsPercent <- rep(FALSE,numsteps)
 art2ndIsPercent <- rep(FALSE,numsteps)
 eligibilityThreshold <- rep(NA,numsteps)

 scaling <- fnCalcScaling()

 #instantaneous transition rates
 lambdaCD4ts <- rep(NA,length(lambdaCD4))
 for(i in 1:length(lambdaCD4)){
  lambdaCD4ts[i] <- 1-exp(-1/lambdaCD4[i] * timestep)
  }
 beta1CD4 <- beta2CD4 <- rep(NA,numCD4compartments-1)         
 for(i in 1:length(beta1CD4)){
  beta1CD4[i] <- tempBeta1;
  beta2CD4[i] <- tempBeta2;
 }

 #Input arrays for number on ART are not for full projection, so need to expand them
 numART1in <- numART2in <- rep(NA,numYears+1)
 art1stIsPercentIn <- art2ndIsPercentIn <- rep(NA,numYears+1)
 eligibilityThresholdCD4 <- rep(NA,numYears+1)
 for(i in 1:(preARTStartYrs+1)){
  numART1in[i] = numART2in[i] = 0.0                     #No ART prior to the earliest ART start year (1995 usually)
  art1stIsPercentIn[i] = art2ndIsPercentIn[i] = TRUE    #Assume 0 is percent until user entered data starts
  eligibilityThresholdCD4[i] = EPPConstantsDEFAULT_CD4_THRESHOLD
 }

 for(i in (preARTStartYrs + 2):length(numART1in)){ #need to shift by 1
  numART1in[i] <- temp1stLine[i-preARTStartYrs-1]
  numART2in[i] <- temp2ndLine[i-preARTStartYrs-1]
  art1stIsPercentIn[i] <- temp1stLineIsPercent[i-preARTStartYrs-1]
  if(!art1stIsPercentIn[i]) numART1in[i] <- numART1in[i] * scaling
  art2ndIsPercentIn[i] <- temp2ndLineIsPercent[i-preARTStartYrs-1]
  if(!art2ndIsPercentIn[i]) numART2in[i] <- numART2in[i] * scaling
  eligibilityThresholdCD4[i] <-tempEligibilityThreshold[i-preARTStartYrs-1]
 }

 #Now set up the number on 1st and 2nd line ART arrays and the 1st year survival array for all timesteps
 #Set first year on ART to zero noting that ART is specified on Dec 31, i.e. at our first timestep
 # in the next year and we have a half year offset
 for(i in 1:5){
  numART1[i] <- numART2[i] <- 0.0
  art1stIsPercent[i] <- art2ndIsPercent[i] <- TRUE
 }
 #Now fill in the rest of them from the values given
 for(i in 1:(numsteps-6)){
  j <- ceiling(i/timestepsInYear)
  numART1[i+5] = numART1in[j] + (numART1in[(j+1)]-numART1in[j])* ((i-1)%%timestepsInYear)/timestepsInYear
  numART2[i+5] = numART2in[j] + (numART2in[(j+1)]-numART2in[j])* ((i-1)%%timestepsInYear)/timestepsInYear
 }
 numART1[numsteps] <- (numART1in[numYears+1]+numART1in[numYears])/2
 numART2[numsteps] <- (numART2in[numYears+1]+numART2in[numYears])/2
   
 for(i in 1:5){
  eligibilityThreshold[i] <- EPPConstantsDEFAULT_CD4_THRESHOLD
  art1stIsPercent[i] <- art2ndIsPercent[i] <- TRUE
 }
 for(i in 1:(numsteps-6)){
  j <- ceiling(i/timestepsInYear)
  eligibilityThreshold[i+5] <- eligibilityThresholdCD4[j]
  art1stIsPercent[i+5] <- art1stIsPercentIn[(j+1)]
  art2ndIsPercent[i+5] <- art2ndIsPercentIn[(j+1)]
 }
 eligibilityThreshold[numsteps] <- eligibilityThresholdCD4[numYears]
 art1stIsPercent[numsteps] <- art1stIsPercentIn[numYears]
 art2ndIsPercent[numsteps] <- art2ndIsPercentIn[numYears]      


return(list(numART1=numART1,numART2=numART2,art1stIsPercent=art1stIsPercent,art2ndIsPercent=art2ndIsPercent,eligibilityThreshold=eligibilityThreshold,lambdaCD4ts=lambdaCD4ts,beta1CD4=beta1CD4,beta2CD4=beta2CD4))
}

ARTarrays <- fnGenARTarrays()
