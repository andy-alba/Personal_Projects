data = read.csv("../data/prostate.csv", header = TRUE)


# Exploring Outliers
summary(data)

head(data[order(data$weight, decreasing = TRUE), ])
# Observation 32 may have an outlier or error for weight

head(data[order(data$psa, decreasing = TRUE), ])
# Observation 95, 96, 97 are most likely influential points 

head(data[order(data$c.vol, decreasing = TRUE), ]) 
# Observation 94 may be an influential point

# Check the DFBeta's and STDR plot for psa, weight, c.vol

library(ggplot2)
ggplot(data, aes(group = y, y = psa )) + geom_boxplot()


ggplot(data, aes(group = y, y = c.vol)) + geom_boxplot() + 
  ylim(0, 30) +
  labs(x = 'No Cancer Diagnosis(Left), Cancer Diagnosis (Right)',
       y = 'Serum Prostate-Specific Antigen Level (mg/ml)',
       title = 'Boxplot of PSA Level for Prostate Cancer Diagnosis'
  )
ggsave('psa_boxplot.png', width = 5, height = 10)

ggplot(data, aes(group = y, y = weight)) + geom_boxplot()
# Since the outlier seems to be very 

ggplot(data, aes(group = y, y = weight)) + geom_boxplot() + 
  ylim(0, 100) + 
  labs(x = 'No Cancer Diagnosis(Left), Cancer Diagnosis (Right)',
       y = 'Weight of Prostate (gm)',
       title = 'Boxplot of Prostate Weight for Prostate Cancer Diagnosis'
  )
ggsave('weight_boxplot.png', width = 5, height = 10)

ggplot(data, aes(group = y, y = benign)) + geom_boxplot() +
  labs(x = 'No Cancer Diagnosis(Left), Cancer Diagnosis (Right)',
       y = 'Amount of Benign Prostatic Hyperplasia (cm^2)',
       title = 'Boxplot of Benign Prostatic Hyperplasia\n for Prostate Cancer Diagnosis')
ggsave('benign_boxplot.png', width = 5, height = 10)


ggplot(data, aes(group = y, y = cap)) + geom_boxplot() +
  labs(x = 'No Cancer Diagnosis(Left), Cancer Diagnosis (Right)',
       y = 'Degree of Capsular Penetration (cm)',
       title = 'Boxplot of Capsular Penetration\n for Prostate Cancer Diagnosis')
ggsave('cap_boxplot.png', width = 5, height = 10)

ggplot(data, aes(group= y, y = age)) + geom_boxplot() +
  labs(x = 'No Cancer Diagnosis(Left), Cancer Diagnosis (Right)',
       y = 'Age of Patient (years)',
       title = 'Age for Prostate Cancer Diagnosis')
ggsave('age_boxplot.png', width = 5, height = 10)

full_model = glm(y~., data = data, family = binomial(link = logit))
empty_model = glm(y~1, data = data, family = binomial(link = logit))

# Exploring Outliers

library(LogisticDx)
dxinfo = dx(full_model)
DFBeta = dxinfo$dBhat
std.r = dxinfo$sPr


plot(DFBeta, main = "Index plot of the change in the Betas")
# Many points with dfbeta > .1, a few where dfbeta > 1

dxinfo[DFBeta > 1]

hist(std.r)
std.r[abs(std.r) >= 3]
dxinfo[abs(std.r) >= 3]



data[data$psa == 21.758,] # index 74
data[data$psa == 56.261,] # index 91

# removing outliers
data = data[-c(74, 91),]

# Model Selection

full_model = glm(y~., data = data, family = binomial(link = logit))
empty_model = glm(y~1, data = data, family = binomial(link = logit))

summary(full_model)
f_AIC = step(empty_model,
             scope = list(lower = empty_model, upper = full_model),
             direction = "forward")
f_AIC
# Forward AIC only took c.vol, psa, age

f_BIC = step(empty_model,
             scope = list(lower = empty_model, upper = full_model),
             direction = "forward", k = log(nrow(data)))
f_BIC
# Forward BIC only took c.vol and psa

b_AIC = step(full_model,
             scope = list(lower = empty_model, upper = full_model),
             direction = "backward")
b_AIC
# Backward AIC only took psa, c.vol, age

b_BIC = step(full_model,
             scope = list(lower = empty_model, upper = full_model),
             direion = "backward", k = log(nrow(data)))
b_BIC

step(full_model,
     scope = list(lower = empty_model, upper = full_model),
     direion = "both", k = log(nrow(data)))
step(empty_model,
     scope = list(lower = empty_model, upper = full_model),
     direion = "both", k = log(nrow(data)))

value = c(f_AIC$aic, b_AIC$aic, f_BIC$aic, b_BIC$aic)
models = c(f_AIC$formula, b_AIC$formula, f_BIC$formula, b_BIC$formula)
value
models

temp = data$y
temp
tempdf = data[,-data$y]
tempdf$y = temp
tempdf
best.subset.BIC = bestglm(Xy = tempdf, family = binomial(link=logit),IC = "BIC",method = "exhaustive")
best.subset.BIC  
# Backward BIC only took c.vol and psa

# Conclusion
# Forward and Backward Subset selection selected the same 
# models when using the same criterion. There is probably 
# no need to use FB or BF to select a model. However, I 
# still need to see if age is significant enough to keep.


interaction1_model = 
  glm(y ~ psa + c.vol + age + inv + c.vol*psa, data = data,
      family = binomial(link = logit))
interaction2_model = 
  glm(y ~ psa + c.vol + age + inv + psa*age, data = data,
      family = binomial(link = logit))
interaction3_model = 
  glm(y ~ psa + c.vol + age + inv + psa*inv, data = data,
      family = binomial(link = logit))
interaction4_model = 
  glm(y ~ psa + c.vol + age + inv + c.vol*age, data = data,
      family = binomial(link = logit))
interaction5_model = 
  glm(y ~ psa + c.vol + age + inv + c.vol*inv, data = data,
      family = binomial(link = logit))
interaction6_model = 
  glm(y ~ psa + c.vol + age + inv + age*inv, data = data,
      family = binomial(link = logit))
best_model = glm(y ~ psa + c.vol + age + inv, data = data,
                 family = binomial(link = logit))

best_model$coefficients

interaction_test = function(inter, no_inter, df0, df1) {
  L0 = logLik(no_inter)
  L1 = logLik(inter)
  
  test_stat = -2 * (as.numeric(L0) - as.numeric(L1))
  p.val = pchisq(test_stat, df= df1 - df0, lower.tail=FALSE)
  return(c(test_stat, p.val))
}

# Testing For:
#   No Interaction vs Interaction

interaction_test(interaction1_model, best_model, 4, 5)
interaction_test(interaction2_model, best_model, 4, 5)
interaction_test(interaction3_model, best_model, 4, 5)
interaction_test(interaction4_model, best_model, 4, 5)
interaction_test(interaction5_model, best_model, 4, 5)
interaction_test(interaction6_model, best_model, 4, 5)


# Make Confidence intervals for B

b_AIC

estimates =  summary(best_model)$coefficients[,1] # A vector of only the estimates
estimates
SE =  summary(best_model)$coefficients[,2]

alpha = 0.01
z.a.2 = qnorm(1-alpha/2)
upper.bounds = estimates +z.a.2*SE
lower.bounds = estimates -z.a.2*SE
Wald.CI = cbind(lower.bounds, upper.bounds)
Wald.CI

# Error Rate

truth = data$y #The true values of y
predicted = ifelse(fitted(best_model)>.5,1,0) #The predicted values of y based on pi.0
my.table = table(truth,predicted) 
sens = sum(predicted == 1 & truth == 1 )/sum(truth == 1)
spec = sum(predicted == 0 & truth == 0 )/sum(truth == 0)
error = sum(predicted != truth)/length(predicted)
results = c(sens,spec,error)
names(results) = c("Sensitivity","Specificity","Error-Rate")
results

# AUC
# 
library(pROC)
my.auc = auc(best_model$y,fitted(best_model),
             plot = TRUE,legacy.axes = TRUE)
my.auc
auc.CI = ci(my.auc,level = 1-0.05) #for a 95% confidence interval.
auc.CI

# PRE 

prop.red = 1- sum((best_model$y -best_model$fitted.values)^2)/sum((best_model$y - mean(best_model$y))^2)
prop.red
prop.red = 1- sum((second_best_model$y -second_best_model$fitted.values)^2)/sum((second_best_model$y - mean(second_best_model$y))^2)
prop.red

# Interpretation

exp(0.22971340) # c.vol
exp(0.08372397 ) # psa
exp(0.17706457 ) # age
exp(2.76134535 ) # inv
exp(-19.09174266) # intercept

Wald.CI = cbind(exp(lower.bounds), exp(upper.bounds))
Wald.CI

# Prediction
best_model$coefficients
b = -19.09174266 + 0.08372397 * 5 + 0.22971340 * 10 +  0.17706457 * 67 
exp(b) / (1 + exp(b))