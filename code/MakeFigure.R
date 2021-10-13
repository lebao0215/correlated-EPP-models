par(mex=0.7, mar = c(4.5, 4.5, 5.5, 2.5))

# source(file = "code/plot_prevNew.R")
if (region=="Urban"){
  prev.old = pnorm(sweep(qnorm(prev1.old), 1, X1.old[,8], "+"))
  prev.new = pnorm(sweep(qnorm(prev1.new), 1, X1.new[,8], "+"))
}
if (region=="Rural"){
  prev.old = pnorm(sweep(qnorm(prev2.old), 1, X2.old[,8], "+"))
  prev.new = pnorm(sweep(qnorm(prev2.new), 1, X2.new[,8], "+"))
}

max_prev = max(dataOut$data$P_table*100, na.rm=T)
max_prev = max(c(max_prev, apply(prev.old,2,quantile,0.975,na.rm=T)*100, apply(prev.new,2,quantile,0.975,na.rm=T)*100), na.rm=T)

country.print = country
if (country=="RDC")
  country.print = "Democratic Republic of the Congo"

na_sites = dataOut$data$P_table;  isna_sites = is.na(na_sites)
nr_sites = dim(na_sites)[1];    nr_years = dim(na_sites)[2]
plot(prev.old[1,]*100~seq(1970,2015), 
     xlab = "Year", lty=2, ylab = "Prevalence(%)", cex.axis=1.2, cex.lab=1.2, cex.main=1.5, 
     main=paste(country.print, region), xlim = c(1970, 2015), ylim = c(0,max_prev), col="red", lwd=2, type = "n")


for (i in 1:nrow(prev.new))
  lines(prev.old[i,]*100~seq(1970,2015), lwd=1, lty=2, col="gray")
# plot.clinic()

for (site in 1:nr_sites)    points(na_sites[site,]*100~seq(data_start_yr,data_end_yr), col = "black", lwd = 2, pch = 1, cex=0.5)

lines(apply(prev.old,2,mean,na.rm=T)*100~seq(1970,2015), lwd=4, col="black")
lines(apply(prev.old,2,quantile,0.025,na.rm=T)*100~seq(1970,2015), lwd=2, lty=2, col="black")
lines(apply(prev.old,2,quantile,0.975,na.rm=T)*100~seq(1970,2015), lwd=2, lty=2, col="black")

lines(apply(prev.new,2,mean,na.rm=T)*100~seq(1970,2015), lwd=4, col="blue")
lines(apply(prev.new,2,quantile,0.025,na.rm=T)*100~seq(1970,2015), lwd=2, lty=2, col="blue")
lines(apply(prev.new,2,quantile,0.975,na.rm=T)*100~seq(1970,2015), lwd=2, lty=2, col="blue")

legend(1970, max_prev, c("Independent Model", "Correlated Model", "Data Average"), text.col = c("black", "blue", "red"), cex=1.2,
       col = c("black", "blue", "red"), lty = c(1, 1, 1)) 


check_1 = apply(!is.na(na_sites[,(ncol(na_sites)-2):ncol(na_sites)]), 1, sum)>0
check_2 = apply(!is.na(na_sites[,1:(ncol(na_sites)-3)]), 1, sum)>0
na_sites = na_sites[which(check_1 & check_2),]
if (!is.vector(na_sites))
  data.avg = apply(na_sites, 2, mean, na.rm=T)
if (is.vector(na_sites))
  data.avg = na_sites
index = which(data.avg>=0)
lines(data.avg[index]*100~seq(dataOut$data_start_yr, dataOut$data_end_yr)[index], lwd=2, col="red")
# points(NPBS$rate[1]*100~NPBS$year[1], col = "red", lwd = 2, pch = 16, cex=1.7)
# round(cbind(1970:2015, apply(prev.old, 2, mean)*100, apply(prev.new, 2, mean)*100),1)
