#Function to load and reorganize ANC surveillance data from csv file
fnGetNPBS <- function(country,u_r){ 

 #Reorganize data
 data <- read.csv(paste(dirData,country,u_r,"NPBS.csv",sep=""), header=TRUE, fill=TRUE, colClasses="character")
 year = as.numeric(data[,1])
 rate = as.numeric(data[,2])
 std = as.numeric(data[,3])
 size = as.numeric(data[,4])

 out <- list(year, rate/100, std/100, size)
 names(out) <- c("year", "rate", "std", "size")
 return(out)
}

NPBS <- fnGetNPBS(country,region)
NPBS.order = sort(NPBS$year, index.return = TRUE)$ix
NPBS$year = NPBS$year[NPBS.order]
NPBS$rate = NPBS$rate[NPBS.order]
NPBS$size = NPBS$size[NPBS.order]

var_NPBS = 2*pi*NPBS$rate*(1-NPBS$rate)/NPBS$size*exp(qnorm(NPBS$rate)^2)
