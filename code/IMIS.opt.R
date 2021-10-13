################################################################################################
# IMIS algorithm draws samples from the posterior distribution. 
# The user has to define the following R functions in advance: 
# 1. prior(x) calculates prior density of x, 
# 2. likelihood(x) calculates the likelihood of x, 
# 3. sample.prior(n) draws n samples from the prior distribution. 
# B = incremental sample size, B.re = resample size, 
# number_k = max number of iterations
################################################################################################
## 1 line of changes marked by ## new ##, Le Bao, 09/14/2012

IMIS.opt <- function(B, B.re, number_k){

	B0 = B*10								# initial sample size
	X_all = X_k = sample.prior(B0)				# Draw initial samples from the prior distribution
	Sig2_global = cov(X_all)					# the prior covariance
	stat_all = matrix(NA, 6, number_k)				# 6 diagnostic statistics at each iteration
	center_all = prior_all = prior_new_all = like_all = c()			# centers of Gaussian components, prior densities, and likelihoods
	sigma_all = list()						# covariance matrices of Gaussian components
	D = 1									# the number of optimizers
	
	iter.start.time = proc.time()
	for (k in 1:number_k ){

		prior_all = c(prior_all, prior(X_k))		# Calculate the prior densities
		# prior_new_all = c(prior_new_all, prior.new(X_k))  	# Calculate the prior densities
		if (k==1){
			ptm.opt = proc.time()
			loglike_all = NULL
			for (i in 1:dim(X_k)[1]){
			  Xi1 = X_k[i,seq(1,16,2)]
			  Xi2 = X_k[i,seq(2,16,2)]
			  loglike_all = c(loglike_all, likelihood.log(Xi1, dataOut1, NPBS1) + likelihood.log(Xi2, dataOut2, NPBS2))			# Calculate the log likelihoods
			}
			ptm.use = (proc.time() - ptm.opt)[3]
			loglike_max = max(loglike_all)
			like_all = exp(loglike_all-loglike_max)
			envelop_all = prior_all				# envelop stores the sampling densities
			print(paste(length(which(like_all>min(like_all))), "non-zero likelihoods, time used=", 
			round(ptm.use/60,2), "minutes"))
		}
		if (k>=2){
			ptm.opt = proc.time()
			loglike_new = NULL
			for (i in 1:dim(X_k)[1])
				loglike_new = c(loglike_new, likelihood.log(X_k[i,]))		# Calculate the likelihoods
			loglike_all = c(loglike_all, loglike_new)
			like_all = c(like_all, exp(loglike_new-loglike_max))
			ptm.use = (proc.time() - ptm.opt)[3]
			envelop_all = apply( rbind(prior_all*B0/B, gaussian_all), 2, sum) / (B0/B+D+(k-2))
		}

		Weights = prior_all * like_all / envelop_all	# importance weight is determined by the posterior density divided by the sampling density
		stat_all[1,k] = log(mean(Weights))			# the raw marginal likelihood
		Weights = Weights / sum(Weights)
		stat_all[2,k] = sum(1-(1-Weights)^B.re)             # the expected number of unique points
		stat_all[3,k] = max(Weights)                                # the maximum weight
		stat_all[4,k] = 1/sum(Weights^2)                    # the effictive sample size
		stat_all[5,k] = -sum(Weights*log(Weights), na.rm = TRUE) / log(length(Weights))     # the entropy relative to uniform
		stat_all[6,k] = var(Weights/mean(Weights))  # the variance of scaled weights
		if (k==1)   print("Stage   MargLike   UniquePoint   MaxWeight   ESS   IterTime")
		iter.stop.time = proc.time()
		print(c(k, round(stat_all[1:4,k], 3), as.numeric(iter.stop.time - iter.start.time)[3]/60))
		
		if (k>=1){
			important = which(Weights == max(Weights))
			if (length(important)>1)	important = important[1]
			X_imp = X_all[important,]									# X_imp is the maximum weight input
			distance_all = mahalanobis(X_all, X_imp, diag(diag(Sig2_global)) )		# Calculate the distances to X_imp
			label_nr = sort(distance_all, decreasing = FALSE, index=TRUE)			# Sort the distances
			which_var = label_nr$ix[1:B]									# Pick B inputs for covariance calculation
			Sig2 = cov.wt(X_all[which_var,], wt = Weights[which_var]+1/length(Weights), cor = FALSE, center = X_imp, method = "unbias")$cov
			center_all = rbind(center_all, X_imp)
			sigma_all[[D+k-1]] = Sig2+diag(sd0^2/10^4)
			X_k = rmvnorm(B, X_imp, sigma_all[[D+k-1]])									# Draw new samples
			X_k[,1:4] = round(X_k[,1:4])			                          ## new ##
			X_all = rbind(X_all, X_k)
		}

		# Calculate the Gaussian densities
		if (k==1)
			gaussian_all = dmvnorm(X_all, X_imp, sigma_all[[1]])
		if (k>1){
			gaussian_new = matrix(0, D+k-1, dim(X_all)[1] )
			gaussian_new[1:(D+k-2), 1:(dim(X_all)[1]-B)] = gaussian_all
			gaussian_new[D+k-1, ] = dmvnorm(X_all, X_imp, sigma_all[[D+k-1]])
			for (j in 1:(D+k-2))	gaussian_new[j, (dim(X_all)[1]-B+1):dim(X_all)[1] ] = dmvnorm(X_k, center_all[j,], sigma_all[[j]])
			gaussian_all = gaussian_new
		}

		if (stat_all[2,k] > (1-exp(-1))*B.re)	break							# Check convergence
	} # end of k
	
  
	Weights.sort = sort(Weights, decreasing=T)
  if (length(Weights)>5000)
    cutoff = max(1e-10, Weights.sort[5000])
	nonzero = which(Weights>=cutoff)
	resample_X = X_all[nonzero,]
	resample_w = Weights[nonzero]
	
	X1 = resample_X[,seq(1,16,2)]
	prev = inc = rt = mu = NULL
	for (i in 1:length(resample_w)){
	  xx = fnRMS2011(X1[i,], turnoverOff=TRUE)
	  prev = rbind(prev, xx$annualFittingPrevs)
	  inc = rbind(inc, xx$annualIncid)
	  rt = rbind(rt, xx$annualR)
	  mu = rbind(mu, xx$annualMu)
	}
	prev1 = prev; inc1 = inc; rt1 = rt; mu1 = mu
	
	X2 = resample_X[,seq(2,16,2)]
	prev = inc = rt = mu = NULL
	for (i in 1:length(resample_w)){
	  xx = fnRMS2011(X2[i,], turnoverOff=TRUE)
	  prev = rbind(prev, xx$annualFittingPrevs)
	  inc = rbind(inc, xx$annualIncid)
	  rt = rbind(rt, xx$annualR)
	  mu = rbind(mu, xx$annualMu)
	}
	prev2 = prev; inc2 = inc; rt2 = rt; mu2 = mu
	
	return(list(stat_all=stat_all, resample=resample_X, weight=resample_w, 
	            prev1 = prev1, inc1=inc1, mu1=mu1, 
	            prev2 = prev2, inc2=inc2, mu2=mu2, 
	            n.eval = sum(!is.na(X_all[,1])),
	            time_used=as.numeric(iter.stop.time - iter.start.time)[3]))

	# Return diagnostic statistics, resamples of parameters, pravelnce, incidence, r(t), and mortality rates
} # end of IMIS

