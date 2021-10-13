#Function to load and reorganize ANC surveillance data from csv file
fnGetData <- function(country,u_r){ 

 #Reorganize data
 data <- read.csv(paste(dirData,country,u_r,"ANC.csv",sep=""), header=TRUE, fill=TRUE, colClasses="character")
 data <- data[,-1]           #cut off site names
 data[is.na(data)] <- "-1"
 data[data==""] = "-1"
 for (j in 1:ncol(data))
	data[,j] = as.numeric(gsub(",","", data[,j]))
 check <- c()
 for(i in 1:ncol(data)){
 	check[i] <- length(unique(data[,i]))==1 && data[1,i]==-1
 }
 firstyear <- as.numeric(substr(colnames(data)[1],start=2,stop=5))
 years <- seq(firstyear,firstyear+ncol(data)-1)
 colnames(data) <- years
 data_start_yr <- years[check==FALSE][1]
 data_end_yr <- years[check==FALSE][length(years[check==FALSE])]
 data <- data[,years %in% data_start_yr:data_end_yr]
 data[data==-1] <- NA
 num_sites <- nrow(data)/2
 num_yrs <- ncol(data)

 #Create matrices of Denominator (N), Prevalance (P), and Numerator (O) (0-1)
 N_table <- data[seq(from=2,by=2,to=nrow(data)),]         #denominator
 P_table <- data[seq(from=1,by=2,to=nrow(data)),]/100     #prevalence
 I_table <- N_table*P_table
 colnames(N_table) <- colnames(P_table) <- colnames(I_table) <- seq(data_start_yr,(data_start_yr+ncol(N_table)-1))

 #Check that all that given observation has numerator and denominator
 count <- length(N_table[is.na(N_table)==FALSE & is.na(I_table)==TRUE])   #number of observations that were inconsistent
 N_table[is.na(N_table)==FALSE & is.na(I_table)==TRUE] <- NA
 #if(count>0) { cat("Data for country",country,"region",u_r,"have",count," observations with missing prevalence but defined denominator","\n"); flush.console() }

 N_table <- as.matrix(N_table)
 I_table <- as.matrix(I_table)
 P_table <- as.matrix(P_table)
  
 out <- list(N_table,I_table,P_table)
 names(out) <- c("N_table","I_table","P_table")
 return(list(data=out,data_start_yr=data_start_yr,data_end_yr=data_end_yr))
}

#Create data objects
dataOut <- fnGetData(country,region)
data <- dataOut$data
data_start_yr <- dataOut$data_start_yr
data_end_yr <- dataOut$data_end_yr

