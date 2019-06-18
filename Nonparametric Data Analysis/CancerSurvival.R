cancer = read.csv('data/Survive.csv')
hist(cancer$Time,  main = 'Histogram of Cancer Survival Time', xlab = 'Cancer Survival Time')

B = 10000
many.boot.theta = sapply(1:B, function(i){
  boot.sample = sample(cancer$Time, nrow(cancer),replace  = TRUE) # Note that we are sampling with replacement.
  theta.i = c(quantile(boot.sample, .10), quantile(boot.sample, .90))
  return(theta.i)
})

df = data.frame(t(many.boot.theta))

# Bootstrap distributions
hist(df$X10., main = 'Histogram of Bootstrap Sample 10th Percentile', xlab = 'Sampled 10th Percentile of Cancer Survival Time')
hist(df$X90., main = 'Histogram of Bootstrap Sample 90th Percentile', xlab = 'Sampled 90th Percentile of Cancer Survival Time')

#observed values

tenth_obs = quantile(cancer$Time, .10) # looks skewed to the right
nDth_obs = quantile(cancer$Time, .90) # looks approx normal, some noticeable gaps

# boot estimates

# boot expected
tenth_boot = mean(df$X10.)
nineth_boot = mean(df$X90.)

# boot bias
tenth_bias = tenth_boot - tenth_obs
nineth_bias = nineth_boot - nDth_obs

# boot SE
tenth_se =  sqrt( var(df$X10.)*(B-1)/B)
nineth_se =  sqrt( var(df$X90.)*(B-1)/B)

# bias corrected est
tenth_obs - tenth_bias
nDth_obs - nineth_bias
# Empirical CI for 90th percentile

alpha = 0.05
ci.empirical = as.numeric(c(2*nDth_obs - quantile(df$X90.,1-alpha/2),
                            2*nDth_obs -quantile(df$X90.,alpha/2) ))

# BCA CI for 10th percentile
library(boot)
theta.fun = function(the.dataset, random.indices){
  boot.data = the.dataset$Time[random.indices]
  theta.i = quantile(boot.data, 0.10)
  return(theta.i)
}
r.bootstraps = boot(cancer,theta.fun, R = 10000)
boot.ci(r.bootstraps,type = "bca" ,conf = 1-alpha)

theta.fun = function(the.dataset, random.indices){
  boot.data = the.dataset$Time[random.indices]
  theta.i = quantile(boot.data, 0.90)
  return(theta.i)
}
r.bootstraps = boot(cancer,theta.fun, R = 10000)
boot.ci(r.bootstraps,type = "bca" ,conf = 1-alpha)

qqnorm(df$X90., main = 'QQ Plot of Bootstrap 90th Percentiles')