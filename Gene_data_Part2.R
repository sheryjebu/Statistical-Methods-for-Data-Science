# Load the dataset
library(readxl)
gene_data <- read_excel("gene_data.xlsx",col_names =FALSE)
View(gene_data)

#------------------------------------------------------
#Task 2: Regression â modelling the relationship between gene expression
#------------------------------------------------------

#Estimate model parameters
gene_mat = data.matrix(gene_data)
y <- gene_mat[,1] #sampling time
x1 <- gene_mat[,2]#x1 gene
x2 <- gene_mat[,3]#x2 gene
x3 <- gene_mat[,4]#x3 gene
x4 <- gene_mat[,5]#x4 gene
x5 <- gene_mat[,6]#x5 gene

#define bias
bias <- matrix(1 , length(y),1)

#---------
#Model 1: 
#---------

M1 <- cbind(x4,x3^2,bias)
M1
#Theta is a column vector
thetaHat1 = solve(t(M1) %*% M1) %*% t(M1) %*% x2
print(thetaHat1)

y_Hat1 = M1 %*% thetaHat1

plot(y,x1, col="blue")
#lines(x1,y_Hat1, col="red")

#compute the model residual (error) sum of squared errors (RSS)
#Residual = Observed value - Predicted value
RSS1 = x2 - y_Hat1
RSS1

SSE1 = norm(RSS1, type = "2")^2
SSE1_g1 = sum((x2-y_Hat1)^2)
print(SSE1)
print(SSE1_g1)

#Compute the log-likelihood function
n <- length(y)

sigmahat2_1 = SSE1/(n-1)
LLHFun_1 = -((n/2)*log(2*pi))-((n/2)*log(sigmahat2_1))-((1/(2*sigmahat2_1))*SSE1)
print(LLHFun_1)
print(sigmahat2_1)

#Compute the Akaike information criterion (AIC) 
AIC_1 = (2*length(thetaHat1))-(2*LLHFun_1)
print(AIC_1)

#Bayesian information criterion (BIC)
BIC_1 = (length(thetaHat1)*log(n))-(2*LLHFun_1)
print(BIC_1)


qqnorm(RSS1)

qqplot(RSS1,x2)

density.default(RSS1)
hist.default(RSS1)
shapiro.test(RSS1)
ks.test(RSS1,x2)

#---------
#Model 2
#---------

M2 <- cbind(x4,x3^2,x5,bias)
M2
#Theta is a column vector
thetaHat2 = solve(t(M2) %*% M2) %*% t(M2) %*% x2
print(thetaHat2)

y_Hat2 = M2 %*% thetaHat2

plot(y,x2, col="blue")
#lines(y,x2, col="blue")

#compute the model residual (error) sum of squared errors (RSS)
RSS2 = x2 - y_Hat2
RSS2

SSE2 = norm(RSS2, type = "2")^2
SSE2_g2 = sum((x2-y_Hat2)^2)
print(SSE2)
print(SSE2_g2)

#Compute the log-likelihood function
n <- length(y)

sigmahat2_2 = SSE2/(n-1)
LLHFun_2 = -((n/2)*log(2*pi))-((n/2)*log(sigmahat2_2))-((1/(2*sigmahat2_2))*SSE1)
print(LLHFun_2)
print(sigmahat2_2)

#Compute the Akaike information criterion (AIC) 
AIC_2 = (2*length(thetaHat2))-(2*LLHFun_2)
print(AIC_2)

#Bayesian information criterion (BIC)
BIC_2 = (length(thetaHat2)*log(n))-(2*LLHFun_2)
print(BIC_2)


qqnorm(RSS2)

qqplot(RSS2,x2)

density.default(RSS2)
hist.default(RSS2)
shapiro.test(RSS2)
ks.test(RSS2,x2)

#--------
#Model 3
#--------

M3 <- cbind(x3,x4,x5^3,bias)
M3
#Theta is a column vector
thetaHat3 = solve(t(M3) %*% M3) %*% t(M3) %*% x2
print(thetaHat3)

y_Hat3 = M3 %*% thetaHat3

plot(y,x2, col="blue")
#lines(y,x2, col="blue")

#compute the model residual (error) sum of squared errors (RSS)
RSS3 = x2 - y_Hat3
RSS3

SSE3 = norm(RSS3, type = "2")^2
SSE3_g3 = sum((x2-y_Hat3)^2)
print(SSE3)
print(SSE3_g3)

#Compute the log-likelihood function
n <- length(y)

sigmahat2_3 = SSE3/(n-1)
LLHFun_3 = -((n/2)*log(2*pi))-((n/2)*log(sigmahat2_3))-((1/(2*sigmahat2_3))*SSE1)
print(LLHFun_3)
print(sigmahat2_3)

#Compute the Akaike information criterion (AIC) 
AIC_3 = (2*length(thetaHat3))-(2*LLHFun_3)
print(AIC_3)

#Bayesian information criterion (BIC)
BIC_3 = (length(thetaHat3)*log(n))-(2*LLHFun_3)
print(BIC_3)


qqnorm(RSS3)

qqplot(RSS3,x2)#x2?

density.default(RSS3)
hist.default(RSS3)
shapiro.test(RSS3)
ks.test(RSS3,x2)

#--------
#Model 4 
#--------

M4 <- cbind(x4,x3^2,x5^3,bias)
M4
#Theta is a column vector
thetaHat4 = solve(t(M4) %*% M4) %*% t(M4) %*% x2
print(thetaHat4)

y_Hat4 = M4 %*% thetaHat4

plot(y,x2, col="blue")
#lines(y,x2, col="blue")

#compute the model residual (error) sum of squared errors (RSS)
RSS4 = x2 - y_Hat4
RSS4

SSE4 = norm(RSS4, type = "2")^2
SSE4_g4 = sum((x2-y_Hat4)^2)
print(SSE4)
print(SSE4_g4)

#Compute the log-likelihood function
n <- length(y)

sigmahat2_4 = SSE4/(n-1)
LLHFun_4 = -((n/2)*log(2*pi))-((n/2)*log(sigmahat2_4))-((1/(2*sigmahat2_4))*SSE1)
print(LLHFun_4)
print(sigmahat2_4)

#Compute the Akaike information criterion (AIC) 
AIC_4 = (2*length(thetaHat4))-(2*LLHFun_4)
print(AIC_4)

#Bayesian information criterion (BIC)
BIC_4 = (length(thetaHat4)*log(n))-(2*LLHFun_4)
print(BIC_4)


qqnorm(RSS4)

qqplot(RSS4,x2)#x2?

density.default(RSS4)
hist.default(RSS4)
shapiro.test(RSS4)
ks.test(RSS4,x2)

#--------
#Model 5
#--------

M5 <- cbind(x4,x1^2,x3^3,bias)
M5
#Theta is a column vector
thetaHat5 = solve(t(M5) %*% M5) %*% t(M5) %*% x2# whyx2
print(thetaHat5)

y_Hat5 = M5 %*% thetaHat5

plot(y,x2, col="blue")
#lines(y,x2, col="blue")

#compute the model residual (error) sum of squared errors (RSS)
RSS5 = x2 - y_Hat5
RSS5

SSE5 = norm(RSS5, type = "2")^2
SSE5_g5 = sum((x2-y_Hat5)^2)
print(SSE5)
print(SSE5_g5)

#Compute the log-likelihood function
n <- length(y)

sigmahat2_5 = SSE5/(n-1)
LLHFun_5 = -((n/2)*log(2*pi))-((n/2)*log(sigmahat2_5))-((1/(2*sigmahat2_5))*SSE1)
print(LLHFun_5)
print(sigmahat2_5)

#Compute the Akaike information criterion (AIC) 
AIC_5 = (2*length(thetaHat5))-(2*LLHFun_5)
print(AIC_5)

#Bayesian information criterion (BIC)
BIC_5 = (length(thetaHat5)*log(n))-(2*LLHFun_5)
print(BIC_5)


qqnorm(RSS5)

qqplot(RSS5,x2)#

density.default(RSS5)
hist.default(RSS5)
shapiro.test(RSS5)
ks.test(RSS5,x2)

---------------------------------------
  #Compare Models for Model Selection
  ---------------------------------------
  #AIC for all individual Model at once
  cat("AIC_1 : ",AIC_1)  
cat("AIC_2 : ",AIC_2)
cat("AIC_3 : ",AIC_3)
cat("AIC_4 : ",AIC_4)
cat("AIC_5 : ",AIC_5)
#Negative AIC indicates less information loss than a positive AIC and therefore a better model.


#BIC for all individual Model at once
cat("BIC_1 : ",BIC_1)  
cat("BIC_2 : ",BIC_2)
cat("BIC_3 : ",BIC_3)
cat("BIC_4 : ",BIC_4)
cat("BIC_5 : ",BIC_5)


#-----------------------------------#
# Train and Test data
#-----------------------------------#
#Split the input and output gene dataset (???? and ????) into two parts: one part used to train the model, the
#other used for testing (e.g. 70% for training, 30persentage for testing). For the selected 'best' model, 1) estimate
#model parameters use the training dataset; 2) compute the model's output/prediction on the testing
#data; and 3) also compute the 95% (model prediction) confidence intervals and plot them (with error
#bars) together with the model prediction, as well as the testing data samples.

gen_train <- as.matrix(gene_data[1:210,])
gen_test  <- as.matrix(gene_data[211:301,])

train_x1 = gen_train[,2]
train_x2 = gen_train[,3]
train_x3 = gen_train[,4]
train_x4 = gen_train[,5]
train_x5 = gen_train[,6]
train_y  = gen_train[,1]

test_x1 = gen_test[,2]
test_x2 = gen_test[,3]
test_x3 = gen_test[,4]
test_x4 = gen_test[,5]
test_x5 = gen_test[,6]
test_y  = gen_test[,1]

#compute the model's output/prediction
train_ThetaBias = matrix(1 , length(train_y),1)

test_ThetaBias = matrix(1 , length(test_y),1)

train_M4 = cbind(train_x4, train_x3^2, train_x5^3,train_ThetaBias)

test_M4 = cbind( test_x4, test_x3^2, test_x5^3,test_ThetaBias)

train_thetaHat4 = solve(t(train_M4) %*% train_M4) %*% t(train_M4) %*% train_x2
print(train_thetaHat4)

train_y_Hat4 = test_M4 %*% train_thetaHat4

train_error4 = test_x2 - train_y_Hat4

train_SSE_4 = sum((train_error4)^2)


print(train_SSE_4)

# error variance sigma^2
train_sigma_2 = train_SSE_4/( length(train_y) - 1 ) 
print(train_sigma_2)

train_cov_thetaHat = train_sigma_2 * (solve(t(train_M4) %*% train_M4))

nw_n = length(test_y)
train_var_y_hat = matrix(0 , nw_n , 1)

for( i in 1:nw_n){
  train_y_i = matrix( train_M4[i,] , 1 , length(train_thetaHat4) ) # y[i,] creates a vector. Convert it to matrix
  train_var_y_hat[i,1] = train_y_i %*% train_cov_thetaHat %*% t(train_y_i) 
}

# Confidance interval  95 percent confidence interval
Conf_Inv = 1.96 * sqrt(train_var_y_hat) 


plot(test_y , train_y_Hat4 , type = "b", col = rep(1:3, each = 10))
segments(test_y, train_y_Hat4-Conf_Inv, test_y, train_y_Hat4+Conf_Inv, col = "red", lty = 'solid')


#-----------
# Task 3
#-----------
#Using 'rejection ABC' method to compute the posterior distributions of the 'selected' regression model
#parameters in Task 2.
#1) You only need to compute 2 parameter posterior distributions -- the 2 parameters with largest
#absolute valuesin your least squares estimation (Task 2.1) of the selected model. Fix all the other
#parameters in your model as constant, by using the estimated values from Task 2.1.
#2) Use a Uniform distribution as prior, around the estimated parameter values for those 2
#parameters (from the Task 2.1). You will need to determine the range of the prior distribution.
#3) Draw samples from the above Uniform prior, and perform rejection ABC for those 2 parameters.
#4) Plot the joint and marginal posterior distribution for those 2 parameters.
#5) Explain your results.

#thetaHat4
Theta_v = thetaHat4
Theta_v
N = 100
#taken min and max around the theta values from model4
theta_1 = runif(N,min=0.0,max=0.8) 
theta_2 = runif(N,min=-1.0,max=0.0) 

# To perform rejection

Accept_V = vector(len=N,mode='numeric')
Reject_V  = vector(len=N,mode='numeric')

SSESimulated_V  = vector(len=N,mode='numeric')
Theta1Sampled_V = vector(len=N,mode='numeric')
Theta2Sampled_V = vector(len=N,mode='numeric')

for(i in 1:N)
{
  Theta_v[1] = theta_1[i]
  Theta_v[2] = theta_2[i]
  Theta_v[3] = 0.3982887 #considered const value based on model4
  Theta_v[4] = 0.2009931 #considered const value based on model4
  
  y_simulated = M4 %*% Theta_v
  SSE_simulated = sum((x2-y_simulated)^2)
  if(SSE_simulated < 5)
    Accept_V[i] = 1
  else
    Reject_V[i] = 0
  
  SSESimulated_V[i] = SSE_simulated
  Theta1Sampled_V[i] = theta_1[i]
  Theta2Sampled_V[i] = theta_2[i]
  
}

theta = data.frame(cbind('X1'<-Accept_V,'X2'<-SSESimulated_V,'X3'<-Theta1Sampled_V,'X4'<-Theta2Sampled_V))
theta
# Prior distribution theta

plot(theta_1, theta_2 , col = c("red"), xlab = "theta_1", ylab ="theta_2")

# Posterior distribution theta_1

hist(theta[theta$X1==0,3], col = "coral4",border = "black", probability = TRUE, main = "Posterior Distribution theta_1")


# Posterior distribution theta_2
hist(theta[theta$X1==0,4], col = "coral4",border = "black", probability = TRUE, main = "Posterior Distribution theta_2")


# Rejected parameters plotting

plot(theta[theta$X1==0,3],theta[theta$X1==0,4], col = c("red"), xlab = "theta_1", ylab ="theta_4")
lines(theta[theta$X1==0,3],theta[theta$X1==0,4], col = "blue",type='p')
