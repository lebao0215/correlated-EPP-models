#These functions pull in Spectrum population projections and ready them for use in EPP

#Population data from Spectrum
popData <-read.csv(paste(dirData,country,'Pops.txt',sep=''),header=FALSE)
colnames(popData) <- c('Year','Pop15to49','Pop15','Pop50','TotalImmigration')
popData <- popData[popData$Year<=endYear,]

#National populations
pop15to49national <- popData[,'Pop15to49']
pop15national <- popData[,'Pop15']
pop50national <- popData[,'Pop50']
netMigrationNational <- popData[,'TotalImmigration']

#Create projection sub-population
percentUrban <- read.csv(paste('EPPUrbanDB09.csv',sep=''),header=TRUE)     #percent of country residing in urban areas
pUrban <- cbind(seq(1950,2050),as.numeric(percentUrban[percentUrban[,1]==country,3:(ncol(percentUrban)-1)]))
pUrban <- pUrban[pUrban[,1]>=beginYear & pUrban[,1]<=endYear,2]/100

# if(region=='Urban') urbanFactor <- pUrban	#new
# if(region=='Rural') urbanFactor <- 1-pUrban	#new
urbanFactor = 1	#new

pop15to49 <- pop15to49national * urbanFactor
pop15  <- pop15national * urbanFactor
pop50  <- pop50national * urbanFactor
netMigration  <- netMigrationNational * urbanFactor

#Interpolated populations (i.e., above are annual values but need values for each time step)
p15to49 <- p15dt <- p50dt <- migr <- surv <- rep(NA,numsteps)
p15to49national <- p15dtNational <- p50dtNational <- migrNational <- survBckgnd <- backgroundMort <- rep(NA,numsteps)


fnSpectrumPops <- function(){
 
 for(i in 1:(numsteps-1)){
  j <- ceiling(i/timestepsInYear)
  p15to49[i] <- pop15to49[j] + (pop15to49[j+1] - pop15to49[j]) * ((i-1)%%timestepsInYear)/timestepsInYear
  p15dt[i] <- (pop15[j]+(pop15[j+1] - pop15[j])*((i-1)%%timestepsInYear)/timestepsInYear) * timestep
  p50dt[i] <- (pop50[j]+(pop50[j+1] - pop50[j])*((i-1)%%timestepsInYear)/timestepsInYear) * timestep
  migr[i]  <- (netMigration[j] + (netMigration[j+1]-netMigration[j])*((i-1)%%timestepsInYear)/timestepsInYear) * timestep
 }
  
 p15to49[numsteps] <- pop15to49[numYears+1]
 p15dt[numsteps] <- pop15[numYears+1]*timestep
 p50dt[numsteps] <- pop50[numYears+1]*timestep
 migr[numsteps] <- netMigration[numYears+1]*timestep

 #Normalize migration so that it produces same number of migrants out each year that Spectrum provides
 scalingFactor <- 1
 for(i in 0:(numYears-1)){
  
  sum <- 0
  for(j in 1:timestepsInYear){
  	sum <- sum + migr[i*timestepsInYear + j]
  }
  
  if(sum!=0){
  	scalingFactor=netMigration[i+1]/sum
  }else{
  	scalingFactor <- 1
  }
  
  for(j in 1:timestepsInYear){
  	migr[i*timestepsInYear+j] <- migr[i*timestepsInYear+j]*scalingFactor
  }
 }
 migr[length(migr)] <-  migr[length(migr)]*scalingFactor

 #surv is the intstantaneous survival rate at each timestep, corrected for external and internal migration, background mortality, etc
 surv <- rep(NA,numsteps)
 for(i in 1:length(p15to49)){
  surv[i] <- (p15to49[i+1] - p15dt[i+1] + p50dt[i+1] - migr[i]) / p15to49[i]
 }
 #For the last year use scale according to trend since no (i+1) data
 surv[length(surv)] <- surv[length(surv)-1] * surv[length(surv)-1] / surv[length(surv)-2]
  
 #bckgndMort is background mortality at the national level, calculated from national population figures
 #NOTE: this is used to calculate non-AIDs Deaths. If excess mortality for IDUs is needed, this routine
 #also adjusts for that in the surv array
   
   #NEED TO ADD FUNCTIONALITY FOR IDU EXCESS MORTALITY HERE IF RUNNING CONCENTRATED EPIDEMICS
  
  bckgndMort <- fnBackgroundMortality(surv=surv)
  
  return(list(bckgndMort=bckgndMort,surv=surv,p15to49=p15to49,p15dt=p15dt,p50dt=p50dt,migr=migr))

}  #end function


#Background mortality based on national population data
fnBackgroundMortality <- function(surv){   
	
 #Interpolate population sizes from Spectrum annual output
 p15to49national <- rep(NA,numsteps)
 p15dtNational <- rep(NA,numsteps)
 p50dtNational <- rep(NA,numsteps)
 migrNational <- rep(NA,numsteps)
 
 for(i in 1:(numsteps-1)){
  j <- ceiling(i/timestepsInYear)
  p15to49national[i] <- pop15to49national[j] + (pop15to49national[j+1] - pop15to49national[j]) * ((i-1)%%timestepsInYear)/timestepsInYear
  p15dtNational[i] <- (pop15national[j]+(pop15national[j+1] - pop15national[j])*((i-1)%%timestepsInYear)/timestepsInYear) * timestep
  p50dtNational[i] <- (pop50national[j]+(pop50national[j+1] - pop50national[j])*((i-1)%%timestepsInYear)/timestepsInYear) * timestep
  migrNational[i]  <- (netMigrationNational[j] + (netMigrationNational[j+1]-netMigrationNational[j])*((i-1)%%timestepsInYear)/timestepsInYear) * timestep
 }
  
 p15to49national[numsteps] <- pop15to49national[numYears+1]
 p15dtNational[numsteps] <- pop15national[numYears+1]*timestep
 p50dtNational[numsteps] <- pop50national[numYears+1]*timestep
 migrNational[numsteps] <- netMigrationNational[numYears+1]*timestep


#Normalize migration so that it produces same number of migrants out each year that Spectrum provides
 scalingFactor <- 1
 for(i in 0:(numYears-1)){
  
  sum <- 0
  for(j in 1:timestepsInYear){
  	sum <- sum + migrNational[i*timestepsInYear + j]
  }
  
  if(sum!=0){
  	scalingFactor= netMigrationNational[i+1]/sum
  }else{
  	scalingFactor <- 1
  }
  
  for(j in 1:timestepsInYear){
  	migrNational[i*timestepsInYear+j] <- migrNational[i*timestepsInYear+j]*scalingFactor
  }
 }
 migrNational[length(migrNational)] <-  migrNational[length(migrNational)]*scalingFactor
 sum(migrNational[11:20]) == netMigrationNational[2]

 #SurvBackground is the intstantaneous survival rate at each timestep at the national level, corrected for external and internal migration, background mortality, etc
 survBckgnd <- rep(NA,numsteps)
 backgroundMort <- rep(NA,numsteps)
 for(i in 1:length(p15to49national)){
  survBckgnd[i] <- (p15to49national[i+1] - p15dtNational[i+1] + p50dtNational[i+1] - migrNational[i]) / p15to49national[i]
  backgroundMort[i] <- 1.0 - survBckgnd[i]
 }
 #For the last year use scale according to trend since no (i+1) data
 survBckgnd[length(survBckgnd)] <- survBckgnd[length(survBckgnd)-1] * survBckgnd[length(survBckgnd)-1] / survBckgnd[length(survBckgnd)-2]
 backgroundMort[length(backgroundMort)] = 1.0 - survBckgnd[length(survBckgnd )]
 
 #NEED TO ADD IDU ADJUSTMENT HERE
  
 return(backgroundMort)

}

spectrumObjects <- fnSpectrumPops()

