## Modified from IMIS package by Le Bao (http://cran.r-project.org/web/packages/IMIS/)

IMIS <- function(B0, B, B.re, number_k, D, opt_iter=0){
  
  X_k <- sample.prior(B0)  # Draw initial samples from the prior distribution
  X_all <- matrix(0, B0 + B*((D-1)*opt_iter+number_k), dim(X_k)[2])
  n_all <- 0
  X_all[1:B0,] <- X_k
  
  Sig2_global = cov(X_all[1:B0,])        # the prior covariance
  stat_all = matrix(NA, 6, number_k)                            # 6 diagnostic statistics at each iteration
  center_all = NULL                                             # centers of Gaussian components
  prior_all = like_all = loglike_all = gaussian_sum = vector("numeric", B0 + B*((D-1)*opt_iter+number_k))
  sigma_all = list()                                            # covariance matrices of Gaussian components
  
  iter.start.time = proc.time()
  for (k in 1:number_k ){
    ##Rprof(paste("IMIS-k", k, ".out", sep=""))
    ptm.like = proc.time()
    prior_all[n_all + 1:dim(X_k)[1]] <-  prior(X_k)
    for (i in 1:dim(X_k)[1]){
      Xi1 = X_k[i,seq(1,16,2)]
      Xi2 = X_k[i,seq(2,16,2)]
      loglike_all[n_all + i] =  likelihood.log(Xi1, dataOut1, NPBS1, data_start_yr1, data_end_yr1) + 
                                likelihood.log(Xi2, dataOut2, NPBS2, data_start_yr2, data_end_yr2)			
      # Calculate the log likelihoods
      like_all[n_all + i] = exp(loglike_all[n_all + i])
    }
    ptm.use = (proc.time() - ptm.like)[3]
    if (k==1)   print(paste(B0, "likelihoods are evaluated in", round(ptm.use/60,2), "minutes"))
    which_pos <- which(like_all[1:(n_all + dim(X_k)[1])] > 0)
    
    if(k>=2){
      if(k <= (opt_iter+1)){
        if(k > 2 && length(which_pos) > n_pos)
          gaussian_sum[(n_pos+1):length(which_pos)] <- rowSums(matrix(sapply(1:(D*(k-2)), function(j)dmvnorm(X_all[which_pos[(n_pos+1):length(which_pos)],,drop=FALSE], center_all[j,], sigma_all[[j]])), nrow=length(which_pos)-n_pos))
        n_pos <- length(which_pos)
        gaussian_sum[1:n_pos] <- gaussian_sum[1:n_pos] + rowSums(matrix(sapply(D*(k-2)+1:D, function(j) dmvnorm(X_all[which_pos,,drop=FALSE], center_all[j,], sigma_all[[j]])),nrow=n_pos))
      } else {
        if(k > 2 && length(which_pos) > n_pos)
          gaussian_sum[(n_pos+1):length(which_pos)] <- rowSums(matrix(sapply(1:((D-1)*opt_iter + k-2), function(j)dmvnorm(X_all[which_pos[(n_pos+1):length(which_pos)],,drop=FALSE], center_all[j,], sigma_all[[j]])), nrow=length(which_pos)-n_pos))
        n_pos <- length(which_pos)
        gaussian_sum[1:n_pos] <- gaussian_sum[1:n_pos] + dmvnorm(X_all[which_pos,], center_all[(D-1)*opt_iter + k-1,], sigma_all[[(D-1)*opt_iter + k-1]])
      }
    }

    if (k==1)   envelop_pos = prior_all[which_pos]        # envelop stores the sampling densities
    if (k>1)    envelop_pos = (prior_all[which_pos]*B0/B + gaussian_sum[1:n_pos]) / (B0/B+(k-1))
    Weights = prior_all[which_pos]*like_all[which_pos] / envelop_pos  # importance weight is determined by the posterior density divided by the sampling density

    stat_all[1,k] = log(mean(Weights)*length(which_pos)/(n_all+dim(X_k)[1]))                  # the raw marginal likelihood
    Weights = Weights / sum(Weights)
    stat_all[2,k] = sum(1-(1-Weights)^B.re)             # the expected number of unique points
    stat_all[3,k] = max(Weights)                                # the maximum weight
    stat_all[4,k] = 1/sum(Weights^2)                    # the effictive sample size
    stat_all[5,k] = -sum(Weights*log(Weights), na.rm = TRUE) / log(length(Weights))     # the entropy relative to uniform
    stat_all[6,k] = var(Weights/mean(Weights))  # the variance of scaled weights
    if (k==1)   print("Stage   MargLike   UniquePoint   MaxWeight   ESS   IterTime")
    iter.stop.time = proc.time()
    print(c(k, round(stat_all[1:4,k], 3), as.numeric(iter.stop.time - iter.start.time)[3]/60))
    #iter.start.time = iter.stop.time
    
    ## choose the important point
    important = which(Weights == max(Weights))
    if (length(important)>1)  important = important[1]
    X_imp <- X_all[which_pos[important],]            # X_imp is the maximum weight input
    center_all <- rbind(center_all, X_imp)
    distance_all <- mahalanobis(X_all[1:(n_all+dim(X_k)[1]),], X_imp, diag(diag(Sig2_global)) )
    label_nr = sort(distance_all, decreasing = FALSE, index=TRUE, method="quick")             # Sort the distances
    which_var = label_nr$ix[1:B]                                                              # Pick B inputs for covariance calculation
    
    ## #########
    weight_close <- Weights[match(which_var, which_pos)]
    weight_close[!which_var %in% which_pos] <- 0
    
    n_all <- n_all + dim(X_k)[1]
    Sig2 <- cov.wt(X_all[which_var,], wt = weight_close+1/n_all, cor = FALSE, center = X_imp, method = "unbias")$cov
    sigma_all[[(D-1)*opt_iter+k]] <- Sig2
    X_k <- rmvnorm(B, X_imp, Sig2)                           # Draw new samples
    X_all[n_all + 1:B,] <- X_k
    
    if (stat_all[2,k] > (1-exp(-1))*B.re)       break
  } # end of k

  if (length(which_pos)<=5000){
    resample_X = X_all[which_pos,]
    resample_w = Weights
  }
    
  if (length(Weights)>5000){
    Weights.sort = sort(Weights, decreasing=T)
    cutoff = max(1e-10, Weights.sort[5000], na.rm=T)
    nonzero = which(Weights>=cutoff)
    resample_X = X_all[which_pos[nonzero],]
    resample_w = Weights[nonzero]
  }
  
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
              time_used=as.numeric(iter.stop.time - iter.start.time)[3],
              n.eval = sum(!is.na(X_all[,1])) ) )

} # end of IMIS
