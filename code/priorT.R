########################################################################################
# Define prior(x) which calculates the prior density of x
# Define sample.prior(n) which draws n samples from the prior distribution
########################################################################################
t0.low = 1969.5;	t0.high = 1990.5
mu0 = c(1980, 20,  0.42, 0.17, 0.46, -0.68, -0.038, 0.14)
sd0 = c(sqrt((1990-1970)^2/12), 4.5, 0.23, 0.07, 0.12, 0.24, 0.009, 0.045)
r0 = c(0.462, 0.127, 0.800, 0.794, 0.714, 0.206, 0.540, 0.797)
# t0, t1, log r0, beta0 ~ beta4

# calculate the prior density
prior <- function(x){
  if (is.vector(x)){
    value = 0
    for (j in 1:8)
      value = value + dmvt(x[c(2*j-1, 2*j)], delta = c(mu0[j],mu0[j]), df=2, 
                           sigma=matrix(c(1, r0[j], r0[j], 1), 2, 2)*sd0[j]^2)
  }
  if (!is.vector(x)){
    value = rep(0, nrow(x))
    for (j in 1:8)
      value = value + dmvt(x[,c(2*j-1, 2*j)], delta = c(mu0[j],mu0[j]), df=2, 
                           sigma=matrix(c(1, r0[j], r0[j], 1), 2, 2)*sd0[j]^2)
  }
  return(exp(value))
}

# draw samples from the prior distribution
sample.prior <- function(n){
  input = NULL
  j = 1
  theta.j = rtmvt(n, mean = c(mu0[j],mu0[j]), sigma=matrix(c(1, r0[j], r0[j], 1), 2, 2)*sd0[j]^2, df=2,
                  lower=c(1969.51, 1969.51), upper=c(1990.49, 1990.49))
  input = cbind(input, round(theta.j))
  j = 2 
  theta.j = rmvt(n, delta = c(mu0[j],mu0[j]), sigma=matrix(c(1, r0[j], r0[j], 1), 2, 2)*sd0[j]^2, df=2)
  input = cbind(input, round(theta.j))
  for (j in 3:8){
    theta.j = rmvt(n, delta = c(mu0[j],mu0[j]), sigma=matrix(c(1, r0[j], r0[j], 1), 2, 2)*sd0[j]^2, df=2)
    input = cbind(input, theta.j)
  }
  return(input)
}
