#Tim Brown's function for approximating differential equations of the disease model for EPP 2011 R16, as implemented in R by Dan Hogan.
#Tim's code contains additional components that allow for turnover in concentrated epidemics, but the following is only for generalized epidemics.

#Notes:
#(1) In an effort to speed up the code a bit, I've replaced some loops with vector notation.
#    Rather than remove them, the orginal loops are commented out, as they may be more interpretable.
#(2) See the "GenObjects.R" and "GetParameters.R" files for many of the objects used in this function.
#    (It is more efficient to create them once, rather than each time the integrator is called from IMIS).

# 7-parameter Mean-Shift model by Le Bao is as follows
## 7 lines of changes marked by ## new ##, Le Bao, 09/14/2012


fnRMS2011 <- function(theta,turnoverOff=T){
# theta[1] is t0, theta[2] is t1, t0+t1 is the starting year of stabilization
# theta[3] is the initial log r(t) at t0; theta[4:7] are beta's

#rVar: vector of force of infection parameter values over time
#iota: pulse that initiates the epidemic
#t0: year epidemic starts
#turnoverOff: if TRUE, then no turnover occurs   			
#Set population values before iterating

 t0 <- round(theta[1]);		t1 <- round(theta[2])							## new ##
 y <- 0.0               			# infected pop this timestep  
 z <- pop15to49[t0-beginYear+1]       	# uninfected pop this timestep
 rVec = rVar = rep(NA, numsteps)
 iota = 0.0025											## new ##
 timestepEpidemicStart <- (t0-tProjStart)*timestepsInYear + 1			## new ##
 timeEquilibrium = (t0+t1-tProjStart)*timestepsInYear + 1				## new ##

 # The following code is the key part
 for(t in timestepEpidemicStart:numsteps){
  # value for r at this timestep
  if (t==timestepEpidemicStart)	rVec[t] = theta[3]				## new ##
  if (t>timestepEpidemicStart){
	r.diff = theta[4] * (theta[5] - rVar[t-1]) + theta[6]*prevs[t-1]		## new ##
	if (t>timeEquilibrium)
		r.diff = r.diff + theta[7]*(t-timeEquilibrium)*timestep*(rVar[t-1]*(1-prevs[t-1])-aidsDeathRate[t-1])		## new ##
	rVec[t] = rVec[t-1]+r.diff
	if (rVec[t]<=-300)	{rVec[t]=-300;	r.diff=rVec[t]-rVec[t-1]}		# to make sure the program will not crash
	if (rVec[t]>=300)		{rVec[t]=300;	r.diff=rVec[t]-rVec[t-1]}
  }
  rVar[t] = exp(rVec[t])

  # The following code are the same as what Dan has
  time <- (t-1)*timestep
  times[t] <- time+tProjStart+0.5   

  #Calculate y as an aggregate variable of the various CD4 compartment populations
  ynoART <- sum(popCD4woART)
  yARTgt1yr <- sum(popCD4onART1Plus)
  yARTlt1yr <- sum(popCD4onARTYr1)
  y <- ynoART + yARTgt1yr + yARTlt1yr
  
  #Store population sizes at each time step
  Za[t] <- z                   #Unifected
  Ya[t] <- y                   #Total infected
  YnoART[t] <- ynoART          #Infected, not on ART
  YARTgt1yr[t] <- yARTgt1yr    #Infected, on ART for more than 1 year
  YARTlt1yr[t] <- yARTlt1yr    #Infected, on ART for less than 1 year
  N[t] <- y + z                #Total population
  prevs[t] <- y / N[t]         #Prevalence
  NumOnART[t] <- yARTgt1yr + yARTlt1yr

  #New infections
  if(t==timestepEpidemicStart){
    newinfects <- (rVar[t] * (ynoART + infectLevelART * (yARTlt1yr + yARTgt1yr))/N[t] + iota) * z * timestep
   }else{
    newinfects <- (rVar[t] * (ynoART + infectLevelART * (yARTlt1yr + yARTgt1yr))/N[t]) * z * timestep
   }
  if (newinfects>0)
  incid[t] <- newinfects/max(N[t] - Ya[t], newinfects)    #Incidence (defined as NewInfections/Susceptibles)
  if (newinfects==0)		incid[t] <- 0

  newInfections[t] <- newinfects
  cumInfectionsToDate <- cumInfectionsToDate + newinfects
  infectionsToDate[t] <- cumInfectionsToDate

  #Calculate non-AIDS deaths for the entire population
  #Note: This assumes that entries, exits and migration all occur at end of individual timestep
  newNonAIDSDeathsThisTimestep <- N[t] * bckgndMort[t]
  cumInternalMigrationToDate <- cumInternalMigrationToDate - N[t]*((1.0 - surv[t]) - bckgndMort[t])
  cumNonAIDSDeathsToDate <- cumNonAIDSDeathsToDate + newNonAIDSDeathsThisTimestep
  nonAIDSDeaths[t] <- newNonAIDSDeathsThisTimestep
  nonAIDSDeathsToDate[t] <- cumNonAIDSDeathsToDate

  #Calculate AIDS mortality in the compartments here before applying migration, ageout, etc.
  #Background mortality is applied before calculating HIV related deaths to deal with competing risks.
  #Multiply by timestep to adjust annual rates to values for a single shorter timestep. These will
  #be subtracted off later in the calculation.
  newAIDSDeathsThisTimestep <- 0.0

  ynoARTAIDSDeaths <- popCD4woART * surv[t] * muCD4 * timestep
  yARTgt1yrAIDSDeaths <- popCD4onART1Plus * surv[t] * alpha3CD4 * timestep  
  newAIDSDeathsThisTimestep <- newAIDSDeathsThisTimestep + sum(ynoARTAIDSDeaths) + sum(yARTgt1yrAIDSDeaths)

  yARTlt1yrAIDSDeaths[,1:5] <- popCD4onARTYr1[,1:5] * surv[t] * alpha1CD4 * timestep
  yARTlt1yrAIDSDeaths[,6:timestepsInYear] <- popCD4onARTYr1[,6:timestepsInYear] * surv[t] * alpha2CD4 * timestep
  newAIDSDeathsThisTimestep <- newAIDSDeathsThisTimestep + sum(yARTlt1yrAIDSDeaths)

  cumAIDSDeathsToDate <- cumAIDSDeathsToDate + newAIDSDeathsThisTimestep + bckgndMort[t]*Ya[t]		## this is not the aids-speific death, but with the background death
  aidsDeaths[t] <- newAIDSDeathsThisTimestep + bckgndMort[t]*Ya[t]
  aidsDeathsToDate[t] <- cumAIDSDeathsToDate
  if (Ya[t]>0)  aidsDeathRate[t] = aidsDeaths[t]/Ya[t] #+bckgndMort[t]
  if (Ya[t]==0)  aidsDeathRate[t] = 0 #+bckgndMort[t]

  #Apply the new infections and demographic changes to the compartments and track desired quantities
  #Calculate the aging out and migration rates
  #Note: this is based on non-HIV pop to give proportion of total pop aging out
  if(t < numsteps){
   a50Rate = p50dt[t+1]/p15to49[t]
  }else{   #Use last calculated value as an approximation
   a50Rate = p50dt[t]/p15to49[t]
  }
  #Note: migration is held fixed at values provided by Spectrum, regardless of HIV impacts
  migRate = migr[t]/N[t]
  totalMigrants <- totalMigrants + migr[t]
  if(t < numsteps){
   exitsAgeouts = p50dt[t+1] * N[t]/p15to49[t]
   totalAgeouts <- totalAgeouts + exitsAgeouts
  }else{
   exitsAgeouts <- p50dt[t] * N[t]/p15to49[t]
   totalAgeouts <- totalAgeouts + exitsAgeouts
  }

  #Now apply timestep survival, aging out, migration to all compartments
  combinedChangeRate <- surv[t] + migRate - a50Rate
  z <- z * combinedChangeRate
  popCD4woART <- popCD4woART * combinedChangeRate
  popCD4onART1Plus <- popCD4onART1Plus * combinedChangeRate 
  popCD4onARTYr1 <- popCD4onARTYr1 * combinedChangeRate                     
 
  popCD4woART <- popCD4woART - ynoARTAIDSDeaths
  popCD4onART1Plus <- popCD4onART1Plus - yARTgt1yrAIDSDeaths
  popCD4onARTYr1 <- popCD4onARTYr1 - yARTlt1yrAIDSDeaths                    

  #Calculate eligibility and newART based on those who will carry over to the next timestep in
  #i.e. after various mortality, ageout and migration factors are taken into account.
  #This is necessary because we will set newART equal to the eligible if more slots
  #are available than eligible people, and that is the number who will survive to the
  #next timestep. This will keep the flows consistent.
  eligibleForART <- 0.0
  numCD4eligible <- 0
  for(j in 1:numCD4compartments){
   if(cd4LowerLimits[j] < eligibilityThreshold[t]){     
     eligibleForART <- eligibleForART + popCD4woART[j]
     numCD4eligible <- numCD4eligible+1
    }
   }
  EligibleForART[t] <- eligibleForART
  survivingOnART <- 0.0
  #Now calculate the need for ART as eligibility plus those on ART surviving to next timestep         
  survivingOnART <- sum(c(popCD4onART1Plus,popCD4onARTYr1))
  needForART <- eligibleForART + survivingOnART
  ARTNeed[t] <- needForART
  
  #Correction for a number to percent transition after the first data (right now
  #interpolated value below is incorrect downward interpolation of the last numerical value
  #and the percentage value. This corrects for that by converting the next 10 entries
  #in the art1stIsPercent array to correct percentages based on current coverage)
  currentCoverage1st <- survivingOnART/needForART * 100.0
  if(t < numsteps){                
   if((!art1stIsPercent[t]) && art1stIsPercent[t+1]){                    
    numART1[t+1] = currentCoverage1st
    for(i in 2:timestepsInYear){
      numART1[t+i] = currentCoverage1st + (numART1[t+11] - currentCoverage1st) * (i-1)/timestepsInYear
    }
   }                
  }

  #Finally calculate the total number going onto 1st line ART in this timestep
  if(t < numsteps){
   if(art1stIsPercent[t+1]){    				
    newART1 <- needForART * numART1[t+1]/100.0 - survivingOnART 
   }else{
    newART1 <- numART1[t+1] - survivingOnART
   }
  }            

  #If new ART exceeds those eligible for ART, downward adjust the new ART
  if(newART1 > eligibleForART){
   newART1 <- eligibleForART
  }
  NumNewART[t] <- newART1

  #Apportion new ART among the HIV+ group not on ART by population proportion
  for(j in 1:numCD4compartments){
   if(cd4LowerLimits[j] < eligibilityThreshold[t] && eligibleForART != 0.0){
    newART1CD4[j] <- newART1 * popCD4woART[j]/eligibleForART
   }else{
    newART1CD4[j] <- 0.0
   }
  }

  #Next do the movements among the various compartments
  #Adjust the X and Z compartments to add in the new entrants
  if(t < numsteps){
   entrants <- p15dt[t+1]
  }else{
   entrants <- p15dt[t]
  }
  totalEntrants <- totalEntrants + entrants
  z <- z + entrants 

  #Move new infections from z to the CD4>500 not on ART compartment. If fewer susceptibles than
  #new infections downward adjust the new infections to deplete the susceptibles
  #CONDITION on z: downward adjust new infections if they deplete susceptible population
  if(newinfects > z){
   newinfects <- z
   newInfections[t] <- z
  }
  #Move new infections from z to the CD4 > 500 no ART compartment
  z <- z - newinfects

  #Migrate HIV+ people not on ART between the CD4 compartments with lambda and subtract off new ART
  #First calculate the entries and exits in each compartment
  entriesNoARTCD4[CD4_GT_500] <- (1.0 - propNewInfectsCD4_350)*newinfects    #This compartment absorbs largest proportion of new infections
  for(j in 1:(numCD4compartments-1)){
   exitsNoARTCD4[j] <- lambdaCD4ts[j]*popCD4woART[j] + newART1CD4[j]            
   entriesNoARTCD4[j+1] <- lambdaCD4ts[j]*popCD4woART[j]
  }
  #Add in entries from new infections to CD4_350_499
  entriesNoARTCD4[CD4_350_499] <- entriesNoARTCD4[CD4_350_499] + propNewInfectsCD4_350 * newinfects
  exitsNoARTCD4[CD4_LT_50] <- newART1CD4[CD4_LT_50]     #This compartment only has exits to ART

  popCD4woART <- popCD4woART + entriesNoARTCD4 - exitsNoARTCD4
 
   
  ########################################################################################################################
  ########################################################################################################################
  
  #"Beta progressions"
  
  #Add feedback effects of beta in moving people up the CD4 chain and then migrate HIV+ people on ART < 1 yr
  #between compartments. First calculate the beta progressions, then move everybody to the next compartment
  #and finally add new infections to the 1st time compartment. Note: beta1 and beta2 now refer to 1st and 2nd
  #6 months of infections per John's latest model. Final compartment moves into the ART > 1 yr group                      

  #DRH: These beta progressions have not been used, so be sure to re-visit if choose to use them

  if(sum(c(beta1CD4,beta2CD4))!=0){         #Save time and do not calculate these progressions if beta is turned off

     print('Beta progressions are on')

  #For general population
  #First calculate the beta feedback
  for(j in 1:numCD4compartments){
   for(k in 1:5){
    if(j < numCD4compartments){
     entriesART1yrCD4[j,k] = popCD4onARTYr1[j+1,k]*beta1CD4[j]*timestep
    }else{
     entriesART1yrCD4[j,k] = 0.0
    }
    if(j > 1){
     exitsART1yrCD4[j,k] = popCD4onARTYr1[j,k]*beta1CD4[j-1]*timestep
    }else{
     exitsART1yrCD4[j,k] = 0.0     #No beta feedback from the CD4 > 500 compartment
    }
   }
   for(k in 6:timestepsInYear){
    if(j < numCD4compartments){
     entriesART1yrCD4[j,k] = popCD4onARTYr1[j+1,k]*beta2CD4[j]*timestep
    }else{
     entriesART1yrCD4[j,k] = 0.0
    }
    if(j > 1){
     exitsART1yrCD4[j,k] = popCD4onARTYr1[j,k]*beta2CD4[j-1]*timestep
    }else{
     exitsART1yrCD4[j,k] = 0.0
    }
   }
  }
  #Now apply the beta feedback
  for(j in 1:numCD4compartments){
   for(k in 1:timestepsInYear){
    popCD4onARTYr1[j,k] <- popCD4onARTYr1[j,k] + entriesART1yrCD4[j,k] - exitsART1yrCD4[j,k]
   }
  }
 
  } #close if statement for calculating beta feedback
  
    #End of Beta progressions

  ########################################################################################################################
  ########################################################################################################################
    
  #Move populations through ART compartments
  
  #General population
  popCD4onART1Plus <- popCD4onART1Plus + popCD4onARTYr1[,timestepsInYear]  
  popCD4onARTYr1[,timestepsInYear:2] <- popCD4onARTYr1[,(timestepsInYear:2)-1]
  popCD4onARTYr1[,1] <- 0.0   #Have now moved this to second step
  popCD4onARTYr1[,1] <- newART1CD4

 } #End looping over time


 #Calculate annual summaries (as generated by EPP)
 numberOfYears <- numYears + 1

 #Annual prevalence, r(t), incidence, and mortality in July 1 of each year
 yrIndex <- c(1,(1:numYears)*timestepsInYear+1)
 annualPrevs <- prevs[yrIndex]	
 annualR <- rVar[yrIndex]
 annualIncid <- incid[yrIndex]/timestep
 annualMu <- aidsDeathRate[yrIndex]/timestep

 return(list(year=times,prev=prevs,annualFittingPrevs=annualPrevs,annualIncid=annualIncid,annualR=annualR,annualMu=annualMu))
}



