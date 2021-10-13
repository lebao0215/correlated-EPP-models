plot.clinic <- function(){

for (site in 1:nr_sites) {
	color = "black"
	type = 1
	if (site ==2) color = "blue"
	if (site ==3) color = "purple"
	if (site ==4) color = "green"
	if (site ==5) color = "red"
	if (site ==6) color = "pink"
	if (site ==7) color = "yellow"
	if (site ==8) color = "maroon"
	if (site ==9) color = "orange"
	if (site ==10) color = "azure"
	if (site ==11) color = "beige"
	if (site ==12) color = "bisque"
	if (site ==13) color = "brown"
	if (site ==14) color = "chocolate"
	if (site ==15) color = "coral"
	if (site ==16) color = "cyan"
	if (site ==17) color = "gold"
	if (site ==18) color = "magenta"
	if (site ==19) color = "navy"
	if (site ==20) color = "darkcyan"
	if (site ==21) color = "darkgreen"
	if (site ==22) color = "darkred"
	if (site ==23) color = "darkorange"
	if (site ==24) color = "darkmagenta"
	if (site ==25) color = "darkblue"
	if (site ==26) color = "deeppink"
	if (site ==27) color = "deepskyblue"
	if (site ==28) color = "darksalmon"
	if (site ==29) color = "darkorchid"
	points(na_sites[site,]*100~seq(data_start_yr,data_end_yr), col = site, lwd = 2, pch = type, cex=0.5)
}
#index = which(apply(na_sites, 2, median, na.rm=T)>=0)
#lines(apply(na_sites, 2, mean, na.rm=T)[index]*100~seq(data_start_yr,data_end_yr)[index], lwd=2, col="red")
}

# pairs plot
panel.hist <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
#   r <- abs(cor(x, y))
	r <- cor(x, y)
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
#    text(0.5, 0.5, txt, cex = cex.cor * r)
	text(0.5, 0.5, txt, cex = cex.cor*1.2)
}
