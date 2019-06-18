# STA 135 Final Project

data = read.csv('Wholesale customers data.csv')
data = data[,-2]
channel = data$Channel
data = data[-1]/1000
data$Channel = as.character(channel)
head(data)

table(data$Channel)

library(ggplot2)

# Summary Plots

ggplot(data, aes(x = Channel, y = Fresh)) + geom_boxplot() + ylim(0, 20) +
  labs(title = 'Boxplot of Fresh\nProduct Sales', y = 'Sales (1000 monetary units)')
ggsave('fresh.png', width = 3, height = 4)
ggplot(data, aes(x = Channel, y = Milk)) + geom_boxplot() + ylim(0, 20) +
  labs(title = 'Boxplot of Milk\nProduct Sales', y = 'Sales (1000 monetary units)')
ggsave('milk.png', width = 3, height = 4)
ggplot(data, aes(x = Channel, y = Grocery)) + geom_boxplot() + ylim(0, 40) +
  labs(title = 'Boxplot of Grocery\nProduct Sales', y = 'Sales (1000 monetary units)')
ggsave('grocery.png', width = 3, height = 4)
ggplot(data, aes(x = Channel, y = Frozen)) + geom_boxplot() + ylim(0, 10) +
  labs(title = 'Boxplot of Frozen\nProduct Sales', y = 'Sales (1000 monetary units)')
ggsave('frozen.png', width = 3, height = 4)
ggplot(data, aes(x = Channel, y = Detergents_Paper)) + geom_boxplot() + ylim(0, 15) +
  labs(title = 'Boxplot of Detergents/Paper\nProduct Sales', y = 'Sales (1000 monetary units)')
ggsave('detergents.png', width = 3, height = 4)
ggplot(data, aes(x = Channel, y = Delicassen)) + geom_boxplot() +  ylim(0, 4)+
  labs(title = 'Boxplot of Delicassen\nProduct Sales', y = 'Sales (1000 monetary units)')
ggsave('deli.png', width = 3, height = 4)


# Simultaneous CI of u1 - u2

sample1 = data[data$Channel == 1, 1:6]
sample2 = data[data$Channel == 2, 1:6]

p = 6

n1 = nrow(sample1)
n2 = nrow(sample2)
n = c(n1, n2)

S1 = cov(sample1)
S2 = cov(sample2)

xbar1 = colMeans(sample1)
xbar2 = colMeans(sample2)

Spooled = ((n1 - 1) * S1 + (n2 -1) * S2)/(n1 + n2 - 2)

simul_ci = function(n, p, xbar, S, i, alpha = 0.05, sample_num=1){
  L = xbar[i]-sqrt(((n-1)*p/(n-p))*qf(1 - alpha,p,n-p))*sqrt(S[i,i]/n)
  U = xbar[i]+sqrt(((n-1)*p/(n-p))*qf(1 - alpha,p,n-p))*sqrt(S[i,i]/n)
  ci_t = c(L , U)
  
  L = xbar[i]-qt((1 - alpha)/(2*p),n-1,lower.tail=F)*sqrt(S[i,i]/n)
  U = xbar[i]+qt((1 - alpha)/(2*p),n-1,lower.tail=F)*sqrt(S[i,i]/n)
  ci_b = c(L, U)
  
  ci = rbind(ci_t, ci_b)
  return(ci)
}

simul_ci_2 = function(n, p, d, Sp, alpha = 0.05){
  
  wd<-sqrt(((n[1]+n[2]-2)*p/(n[1]+n[2]-p-1))*qf(1-alpha,p,n[1]+n[2]-p-1))*sqrt(diag(Sp)*sum(1/n))
  wd.b<- qt(1-alpha/(2*p),n[1]+n[2]-2) *sqrt(diag(Sp)*sum(1/n))
  Cis<-cbind(d-wd,d+wd)
  Cis.b<-cbind(d-wd.b,d+wd.b)
  
  return(rbind(Cis, Cis.b))
}

round(simul_ci_2(n = c(n1, n2), p = 6, d = xbar1 - xbar2, Sp = Spooled),2)

# Hotellings T^2

alpha = 0.05
T_2 = t(xbar1 - xbar2) %*% solve(sum(1/n)*Spooled) %*% (xbar1 - xbar2)
c_val =  (sum(n)-2) * p/(sum(n)-p-1) * qf(1 - alpha , p, sum(n)-p-1)
T_2
c_val
pf(T_2, p, sum(n) - p - 1, lower.tail=F)

# PCA 

sample_pca = princomp(data[ , 1:6])
summary(sample_pca, loadings = T) # PCA Summary
(sample_pca$sdev)^2 # Eigen Vlaues

# Scree plot
plot(1:(length(sample_pca$sdev)),  (sample_pca$sdev)^2, type='b', 
     main="Scree Plot", xlab="Number of Components", ylab="Eigenvalue Size")

# PCA Plot
plot(x = sample_pca$scores[,1], y = sample_pca$scores[,2], xlab="PC 1", ylab="PC 2")
biplot(sample_pca, main = 'Effects of Features on Principal Components')

data$PCA1 = sample_pca$scores[,1]
data$PCA2 = sample_pca$scores[,2]

data$Channel = as.character(data$Channel)
ggplot(data, aes(PCA1, PCA2, color = Channel)) + geom_point() + xlim(-10, 25) + ylim(-25, 10)

# LDA
library(MASS)
lda.obj = lda(Channel ~ PCA1 + PCA2, data, prior =c (1,1)/2)
plda<-predict(object = lda.obj , newdata = data)
table(data$Channel, plda$class)

gmean <- lda.obj$prior %*% lda.obj$means
const <- as.numeric(gmean %*%lda.obj$scaling)
slope <- - lda.obj$scaling[1] / lda.obj$scaling[2]
intercept <- const / lda.obj$scaling[2]

#Plot decision boundary
plot(data$PCA1, data$PCA2, pch=rep(c(18,20),each=50), col=rep(c(2,4)))
ggplot(data, aes(PCA1, PCA2, shape = Channel)) + geom_point() + 
  xlim(-10, 25) + ylim(-25, 10) + geom_abline(intercept = intercept, slope = slope) + 
  labs(title = 'Principal Component Plot with Separation Line', x ='Principal Component 1', 'Principal Component 2')
ggsave('pc.png', width = 10, height = 10)
