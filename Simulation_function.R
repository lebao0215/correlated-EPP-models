library(mvtnorm)
library(MASS)
library("mosaic")
library(wordspace)


Evaluation <- function(p, n1, n2, rho, variance, n.k){
  result = NULL
  # p = 2;    n1 = 5;   n2 = 1000;  rho = 0.8;    n.k = 10; variance=0.01
  # define the covariance
  Sigma0 = diag(rep(1,2*p))
  for (i in 1:p) Sigma0[i,i+p] = Sigma0[i+p,i] = rho
  Sigma0 = Sigma0*variance
  Sigma0.inv = solve(Sigma0)
  SigmaX = SigmaX.inv = diag(rep(1,2*p))
  sample.size = c(rep(n1,p),rep(n2,p))
  n.mc = 5000
  # w1 = (n1*(n2+1) - n1*n2*rho^2)/((n1+1)*(n2+1) - n1*n2*rho^2)
  # w2 = (n2*rho)/((n1+1)*(n2+1) - n1*n2*rho^2)
  
  # repeat the analysis n.k times to reduce the noise
  MAE.cor = MAE.ind = MAE.new = n.unique = NULL
  for (k in 1:n.k){
    # generate data
    set.seed(k*1000)
    mu = mvrnorm(1, rep(0,2*p), Sigma0)
    set.seed(k*1000+1)
    X1.bar = apply(mvrnorm(n1, mu[1:p], diag(rep(1,p))),2,mean)
    X2.bar = apply(mvrnorm(n2, mu[(p+1):(2*p)], diag(rep(1,p))),2,mean)
    
    # draw posterior samples
    COV.hat.cor = solve(Sigma0.inv + diag(sample.size)*SigmaX.inv)
    x.hat.cor = COV.hat.cor %*% SigmaX.inv %*% diag(sample.size) %*% c(X1.bar,X2.bar)
    COV.hat.ind = diag(1/(sample.size+1))
    x.hat.ind = COV.hat.ind %*% SigmaX.inv %*% diag(sample.size) %*% c(X1.bar,X2.bar)
    set.seed(k*1000+2)
    MC.cor = mvrnorm(n.mc, x.hat.cor, COV.hat.cor)      # directly from cor model
    MC.ind = mvrnorm(n.mc, x.hat.ind, COV.hat.ind)      # directly from ind model
    # MC samples from the first area
    Theta1.cor = MC.cor[,1:p]
    Theta1.ind = MC.ind[,1:p]
    
    # draw resamples using importance sampling
    B.re = 10^6
    Index1 = sample(n.mc, B.re, replace = TRUE)
    Index2 = sample(n.mc, B.re, replace = TRUE)
    MC.samples = cbind(MC.ind[Index1,1:p], MC.ind[Index2,(p+1):(2*p)])
    
    if (p<=15)   MC.ratio = dmvnorm(MC.samples, x.hat.cor, COV.hat.cor, log = FALSE)/dmvnorm(MC.samples, x.hat.ind, COV.hat.ind, log = FALSE)
    if (p>15){
      MC.ratio = NULL
      for (j in 1:100){
        index.j = 1:10^4 + (j-1)*10^4
        MC.ratio = c(MC.ratio, dmvnorm(MC.samples[index.j,], x.hat.cor, COV.hat.cor, log = FALSE)/dmvnorm(MC.samples[index.j,], x.hat.ind, COV.hat.ind, log = FALSE))
      }
    }
    MC.ratio[MC.ratio>100] = 100
    Reweight = MC.ratio/sum(MC.ratio)
    B.new = 1000
    Index.new = sample(B.re, B.new, replace = TRUE, prob = Reweight)
    # MC samples from the first area
    MC.new = MC.samples[Index.new,]
    Theta1.new = MC.new[,1:p]
    if (p==1)     n.unique = c(n.unique, length(unique(Theta1.new)))
    if (p>1)      n.unique = c(n.unique, length(unique(apply(Theta1.new,1,prod))))
    
    # calculate average MAE over all parameters in the first area
    MAE.cor = c(MAE.cor, mean(abs(apply(MC.cor,2,mean) - mu)[1:p]))
    MAE.ind = c(MAE.ind, mean(abs(apply(MC.ind,2,mean) - mu)[1:p]))
    MAE.new = c(MAE.new, mean(abs(apply(MC.new,2,mean) - mu)[1:p]))
  }
  result = cbind(MAE.cor, MAE.ind, MAE.new, n.unique)
  #result$MAE.cor = MAE.cor
  #result$MAE.ind = MAE.ind
  #result$MAE.new = MAE.new
  #result$n = n.unique
  return(result)
}


Evaluation.K <- function(p, n, rho, variance, n.k){
  result = NULL
  K = length(n)
  # n = c(5,500);  K = length(n);    rho = 0.8;    n.k = 10; variance=1; p=2
  # define the covariance
  # Sigma1 is the covariance of theta_1 across areas
  Sigma1 = matrix(rho, K, K)
  diag(Sigma1) = 1
  Sigma1 = Sigma1*variance
  Sigma1.inv = solve(Sigma1)
  # Sigma0 is the full covariance, K areas each of which has p dimensions
  Sigma0 = diag(rep(1,K*p))
  for (j in 1:p)  Sigma0[((j-1)*K+1):(j*K),((j-1)*K+1):(j*K)] = Sigma1
  
  n.mc = 5000
  # w1 = (n1*(n2+1) - n1*n2*rho^2)/((n1+1)*(n2+1) - n1*n2*rho^2)
  # w2 = (n2*rho)/((n1+1)*(n2+1) - n1*n2*rho^2)
  
  # repeat the analysis n.k times to reduce the noise
  MAE.cor = MAE.ind = MAE.new = n.unique = NULL
  for (k in 1:n.k){
    # generate theta
    set.seed(k*1000)
    theta = mvrnorm(1, rep(0,K*p), Sigma0)
    # sample mean of each area
    X.bar = rnorm(p*K, theta, rep(sqrt(1/n),p))
    
    # draw resamples using importance sampling
    B.re = 10^6;    MC.samples = Index = NULL;    MC.ratio = rep(1,B.re)
    for (i in 1:K)
      Index = rbind(Index, sample(n.mc, B.re, replace = TRUE))
    
    # draw posterior samples
    MC.cor = MC.ind = NULL
    COV.hat.cor = solve(Sigma1.inv + diag(n))
    COV.hat.ind = diag(1/(n+1))
    for (j in 1:p){
      X.bar.j = X.bar[((j-1)*K+1):(j*K)]
      mu.hat.cor = COV.hat.cor %*% diag(n) %*% X.bar.j
      mu.hat.ind = COV.hat.ind %*% diag(n) %*% X.bar.j
      set.seed(k*1000+j)
      MC.cor.j = mvrnorm(n.mc, mu.hat.cor, COV.hat.cor)
      MC.ind.j = mvrnorm(n.mc, mu.hat.ind, COV.hat.ind)
      MC.cor = cbind(MC.cor, MC.cor.j)      # directly from cor model
      MC.ind = cbind(MC.ind, MC.ind.j)      # directly from ind model
      # resamples for proposed model
      MC.samples.j = NULL
      for (i in 1:K)
        MC.samples.j = cbind(MC.samples.j, MC.ind.j[Index[i,],i])
      # MC.samples = cbind(MC.samples, MC.samples.j)
      MC.ratio.j = dmvnorm(MC.samples.j, mu.hat.cor, COV.hat.cor, log = FALSE)/dmvnorm(MC.samples.j, mu.hat.ind, COV.hat.ind, log = FALSE)
      MC.ratio = MC.ratio * MC.ratio.j
    }
    # MC.ratio[MC.ratio>100] = 100
    Reweight = MC.ratio/sum(MC.ratio)
    B.new = 1000
    Index.new = sample(B.re, B.new, replace = TRUE, prob = Reweight)
    
    # MC samples from the first area
    Index1 = seq(1,p*K,K)
    MC.new = MC.ind[Index[1,Index.new],Index1]
    theta1 = theta[Index1]
    MC.cor.mean = apply(MC.cor[,Index1],2,mean)
    MC.ind.mean = apply(MC.ind[,Index1],2,mean)
    MC.new.mean = apply(MC.new,2,mean)
    
    # calculate average MAE over all parameters in the first area
    MAE.cor = c(MAE.cor, mean(abs(MC.cor.mean - theta1)))
    MAE.ind = c(MAE.ind, mean(abs(MC.ind.mean - theta1)))
    MAE.new = c(MAE.new, mean(abs(MC.new.mean - theta1)))
    
    if (p==1)     n.unique = c(n.unique, length(unique(MC.new)))
    if (p>1)      n.unique = c(n.unique, length(unique(apply(MC.new,1,prod))))
  }
  result = cbind(MAE.cor, MAE.ind, MAE.new, n.unique)
  #result$MAE.cor = MAE.cor
  #result$MAE.ind = MAE.ind
  #result$MAE.new = MAE.new
  #result$n = n.unique
  return(result)
}


