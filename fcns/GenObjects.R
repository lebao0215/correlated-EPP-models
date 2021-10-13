
#Generate population arrays
source(paste(dirFcns,'fnSpectrumPops.R',sep=''))
p15to49 <- spectrumObjects$p15to49    #adults 15 to 49
p15dt <- spectrumObjects$p15dt        #entrants (new adults)
p50dt <- spectrumObjects$p50dt        #exits (aging out of 15 to 49)
migr <- spectrumObjects$migr          #migrants
surv <- spectrumObjects$surv
bckgndMort <- spectrumObjects$bckgndMort  #background mortality

#Generate treatment arrays
source(paste(dirFcns,'fnSpectrumART.R',sep=''))
numART1 <- ARTarrays$numART1
numART2 <- ARTarrays$numART2
art1stIsPercent <- ARTarrays$art1stIsPercent       #Size of pop receiving ART can be specified as number of individuals or percent
art2ndIsPercent <- ARTarrays$art2ndIsPercent
eligibilityThreshold <- ARTarrays$eligibilityThreshold    #CD4 count for eligibility
lambdaCD4ts <- ARTarrays$lambdaCD4ts    #transition rates between CD4 compartments
beta1CD4 <- ARTarrays$beta1CD4          #"beta progression" that allows people to improve CD4 while on ART. Not in use.
beta2CD4 <- ARTarrays$beta2CD4


####################################################################################################################################
####################################################################################################################################

#Objects used in the disease model (the "integrator" function). 
#To reduce computation time in IMIS, they are created here once, instead of each time the model integrator function is called.
#Note: includes objects for turnover, which are not included in main integrator yet

# dur <- 0                 #Flow of turnover population
# turnoverOff <- 1         #hard code no turnover for now
turnoverOff <- 0         #hard code turnover for now

#Population sizes at different time steps
Ya <- rep(NA,numsteps); aidsDeathRate=rep(NA,numsteps)
YnoART <- YARTlt1yr <- YARTgt1yr <- rep(NA,numsteps)       #Instantaneous totals in CD4 compartments
Za <- rep(NA,numsteps)
N <- rep(NA,numsteps)
YEx <- rep(NA,numsteps)
YExNew <- rep(NA,numsteps)

#Turnover related
exNoART <-  exARTlt1yr <- exARTgt1yr <- rep(NA,numsteps)

#Arrays for storing ART relevant information
NumOnART <- rep(NA,numsteps)           #Total number on ART at different time steps
EligibleForART <- rep(NA,numsteps)     #Total number eligible for ART at different time steps
ARTNeed <- rep(NA,numsteps)
NumNewART <- rep(NA,numsteps)

#Now create the arrays used to extract the outputs needed
prevs <- rep(NA,numsteps)   #prevalence
incid <- rep(NA,numsteps)   #incidence
times <- rep(NA,numsteps)   #times
newInfections <- rep(NA,numsteps)      #new infections in a timestep
infectionsToDate <- rep(NA,numsteps)   #infections to date (used for annual infections)
aidsDeaths <- rep(NA,numsteps)         #aidsDeaths in a timestep
aidsDeathsToDate <- rep(NA,numsteps)   #cumulative AIDS deaths (used for annual deaths to AIDS)
nonAIDSDeaths <- rep(NA,numsteps)      #nonAIDSDeaths in a timestep
nonAIDSDeathsToDate <- rep(NA,numsteps)   #cumulative nonAIDS deaths (used for annual non-AIDS deaths)

#Turnover related
cumExHIV <- rep(NA,numsteps)

#Annual results (as calculated by Tim)
annualPrevs <- rep(0,numYears+1)
annualIncid <- rep(0,numYears+1)
annualMu <- rep(0,numYears+1)

#Vector indices for CD4 compartments
CD4_GT_500 <- CD4_GT_500in
CD4_350_499 <- CD4_350_499in
CD4_LT_50 <- CD4_LT_50in

#ART compartment arrays needed here
#First the queue for holding 1st year ART entries, then the with and without ART CD4 compartments
popCD4onARTYr1 <- matrix(0,nrow=numCD4compartments,ncol=timestepsInYear)      
popCD4woART <- rep(0,numCD4compartments)
popCD4onART1Plus <- rep(0,numCD4compartments)

#Turnover related
exCD4onARTYr1 <- matrix(0,nrow=numCD4compartments,ncol=timestepsInYear)
exCD4woART <- rep(0,numCD4compartments)
exCD4onART1Plus <- rep(0,numCD4compartments)
 
#Arrays for AIDS-related death calculations in compartments
ynoARTAIDSDeaths <- rep(NA,numCD4compartments)
yARTlt1yrAIDSDeaths <- matrix(NA,nrow=numCD4compartments,ncol=timestepsInYear)
yARTgt1yrAIDSDeaths <- rep(NA,numCD4compartments)

#Turnover related arrays for AIDS-death calculations
exNoARTAIDSDeaths <- rep(NA,numCD4compartments)
exARTlt1yrAIDSDeaths <- matrix(NA,nrow=numCD4compartments,ncol=timestepsInYear)
exARTgt1yrAIDSDeaths <- rep(NA,numCD4compartments)

#Turnover related
cumTurnoverEntries <- rep(NA,numsteps)
cumTurnoverExits <- rep(NA,numsteps)

newART1CD4 <- rep(NA,numCD4compartments)       #Number in a CD4 compartment newly receiving ART, index [CD4]
exitsNoARTCD4 <- rep(NA,numCD4compartments)    #Number exiting this timestep from a CD4 compartment, index [CD4]
entriesNoARTCD4 <- rep(NA,numCD4compartments)  #Number entering this timestep to a CD4 compartment, index [CD4]
exitsART1yrCD4 <- matrix(NA,nrow=numCD4compartments,ncol=timestepsInYear)   #Number of exits this timestep from CD4 compartment <1yr, index [CD4][timestepInYr]
entriesART1yrCD4 <- matrix(NA,nrow=numCD4compartments,ncol=timestepsInYear) #Number of entries this timestep to CD4 compartment <1yr, index [CD4][timestepInYr]
 
#Turnover related
yEx = 0.0   #Instantaneous number of surviving ex-group members
newExART1CD4 <- rep(NA,numCD4compartments)  #Number in a CD4 compartments newly receiving ART among ex-group members [CD4]
exitsExNoARTCD4 <- rep(NA,numCD4compartments)
entriesExNoARTCD4 <- rep(NA,numCD4compartments)
exitsExART1yrCD4 <- matrix(NA,nrow=numCD4compartments,ncol=timestepsInYear)
entriesExART1yrCD4 <- matrix(NA,nrow=numCD4compartments,ncol=timestepsInYear)
totalTurnoverEntrants <- 0.0
totalTurnoverExits <- 0.0
      
#Scalars for tracking movement, infections and deaths
cumInfectionsToDate <- 0.0           #Cumulative HIV infections during projection
cumAIDSDeathsToDate <- 0.0           #Cumulative deaths due to HIV during projection
cumNonAIDSDeathsToDate <- 0.0        #Cumulative deaths to non-HIV related causes
cumulativeExHIVPositives <- 0.0      #Cumulative number of HIV positives leaving by turnover
totalEntrants <- 0.0
totalMigrants <- 0.0
totalAgeouts <- 0.0
cumInternalMigrationToDate <- 0.0

#Annual summaries
annualPopWithHIV <- annualNumWithHIV <- annualNumNewHIV <- annualCumInfections <- annualAIDSDeaths <- annualNonAIDSDeaths <- annualCumAIDSDeaths <- annualCumNonAIDSDeaths <- years <- rep(NA,numYears+1)
annualNewExHIV <- annualCurExHIV <- annualCumExHIV <- rep(NA,numYears+1)


